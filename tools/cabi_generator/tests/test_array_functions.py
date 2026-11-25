# SPDX-License-Identifier: MIT
"""Independent scalar-result array-function generation and compiled ABI checks."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from test_generator import run
import generate
import generate_python
import pipeline
import discover

SOURCE = """module example_vectors
use iso_c_binding, only: c_double
implicit none
private
public :: vector_score
contains
function vector_score(x, mean, cov, enabled) result(y)
    !! @dimension n = size(x)
    !! @require n >= 1
    implicit none
    real(c_double), intent(in) :: x(:) !! shape (n); point
    real(c_double), intent(in), optional :: mean(:) !! shape (n); center
    real(c_double), intent(in), optional :: cov(:,:) !! shape (n,n); matrix
    logical, intent(in), optional :: enabled !! include matrix contribution
    real(c_double) :: y
    y = sum(x) + sum(mean)
    if (enabled) y = y + sum(cov)
end function
end module
"""


class ArrayFunctionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root/'vectors.f90'
        self.source.write_text(SOURCE, encoding='ascii')
        self.config = dict(schema_version=2, source='vectors.f90', module='example_vectors',
            kind_module='iso_c_binding', kind_name='c_double', wrapper_module='example_array_api',
            header_guard='EXAMPLE_ARRAY_API_H', functions={}, excluded={},
            outputs=dict(fortran='adapter.f90', header='adapter.h', coverage='coverage.json'),
            argument_docs=dict(n='dimension', x='point', mean='center', cov='matrix',
                               enabled='logical flag', value='scalar result', status='call status'),
            array_functions=dict(vector_score=dict(symbol='example_score',
                dimensions=[dict(name='n',minimum=1)], shapes=dict(x=['n'],mean=['n'],cov=['n','n']),
                status_argument='status', result_argument='value', invalid_status=1,
                description='Evaluate a vector score.', optional_policy='required',
                logical_inputs=dict(enabled={'false':0,'true':1}), shape_contract=True)))

    def test_metadata_and_native_call(self):
        outputs = generate.generate(self.config,self.root)
        source = outputs[self.root/'adapter.f90']
        self.assertIn('value=cg_native_0(',source)
        self.assertNotIn('call cg_native_0',source)
        record = json.loads(outputs[self.root/'coverage.json'])['procedures'][0]
        self.assertEqual(record['interface'],'explicit_shape_function')
        self.assertEqual(record['optional_promoted_to_required'],['mean','cov','enabled'])
        self.assertEqual(generate.generate(self.config,self.root),outputs)

    def test_rejects_unsupported_and_conflicting_policies(self):
        original = copy.deepcopy(self.config)
        for update in (dict(result_argument='x'),dict(result_argument='status'),
                       dict(optional_policy='inferred'),
                       dict(logical_inputs=dict(value={'false':0,'true':1}))):
            self.config = copy.deepcopy(original)
            self.config['array_functions']['vector_score'].update(update)
            with self.assertRaises(ValueError):
                generate.generate(self.config,self.root)
        self.config = original
        for source in (SOURCE.replace('real(c_double) :: y','real(c_double) :: y(2)'),
                       SOURCE.replace('intent(in) :: x','intent(inout) :: x'),
                       SOURCE.replace('function vector_score','impure function vector_score',1)):
            self.source.write_text(source,encoding='ascii')
            with self.assertRaises(ValueError):
                generate.generate(self.config,self.root)

    def test_python_defaults_prefix_and_validation(self):
        (self.root/'arrays.json').write_text(json.dumps(self.config))
        # A distribution fixture is required by the existing Python policy schema.
        from test_generator import HERE, fixture_config
        (self.root/'arithmetic.f90').write_bytes((HERE/'fixtures/arithmetic.f90').read_bytes())
        scalar = fixture_config()
        scalar['outputs']['coverage'] = 'scalar_coverage.json'
        (self.root/'scalar.json').write_text(json.dumps(scalar))
        policy = dict(schema_version=1,output_directory='python/example',
            groups=dict(arithmetic=dict(manifest='scalar.json',prefix='',defaults=dict(shift=0.,gain=1.),existing_functions={})),
            array_groups=dict(vectors=dict(manifests=['arrays.json'],prefix='vector_',
                defaults=dict(vector_score=dict(mean={'fill':'zeros'},cov={'fill':'identity'},enabled=False)))))
        path = self.root/'python.json'
        def emit():
            path.write_text(json.dumps(policy))
            return generate_python.artifacts(path,self.root)
        outputs = emit()
        text = outputs[self.root/'python/example/_generated.py']
        self.assertIn('def score(self, x, mean=None, cov=None, enabled=False):',text)
        self.assertIn("'result_argument': 'value'",text)
        compile(text,'<generated>','exec')
        batch = self.root/'batch.json'
        batch.write_text('["scalar.json", "arrays.json"]')
        plan_path = self.root/'pipeline.json'
        plan_path.write_text(json.dumps(dict(schema_version=1,c_manifests='batch.json',
            python_policy='python.json',source_directories=['.'],build_manifest='plan.json',
            coverage_report='interfaces.json',aggregate_header='all.h',aggregate_header_guard='ALL_H')))
        artifacts = pipeline.artifacts(plan_path,self.root)
        report = json.loads(artifacts[self.root/'interfaces.json'])
        row = next(r for r in report['procedures'] if r['procedure']=='vector_score')
        self.assertEqual(row['array_adapter'],'example_score')
        self.assertEqual(row['python'],'vectors.score')
        configured = discover.configured(self.root,batch,path)
        self.assertEqual(configured[('example_vectors','vector_score')]['array_function'],'example_score')
        self.assertEqual(configured[('example_vectors','vector_score')]['python'],'vectors.score')
        for value in ({'fill':'random'}, {'fill':'identity'}, None):
            policy['array_groups']['vectors']['defaults']['vector_score']['mean'] = value
            with self.assertRaises(ValueError):
                emit()
        policy['array_groups']['vectors']['defaults']['vector_score']['mean'] = {'fill':'zeros'}
        policy['array_groups']['vectors']['defaults']['vector_score']['x'] = {'fill':'zeros'}
        with self.assertRaises(ValueError):
            emit()

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc'),'GNU compilers required')
    def test_compiled_checked_optimized_and_c_header(self):
        outputs = generate.generate(self.config,self.root)
        for path,text in outputs.items():
            path.write_text(text,encoding='ascii')
        client = self.root/'client.c'
        client.write_text("""#include "adapter.h"
int main(void) {
 double x[2]={1,2}, mean[2]={3,4}, cov[4]={1,0,0,1}, value=99; int status;
 example_score(2,x,mean,cov,1,&value,&status);
 if (status || value != 12) return 1;
 example_score(2,x,mean,cov,0,&value,&status);
 if (status || value != 10) return 2;
 value=99; example_score(0,x,mean,cov,1,&value,&status);
 if (status != 1 || value != 99) return 3;
 example_score(2,x,mean,cov,7,&value,&status);
 return status != 1 || value != 99;
}
""",encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'],['-O3'])):
            folder = self.root/str(index); folder.mkdir()
            run(['gfortran',*flags,'-std=f2018','-Werror=implicit-interface','-c',self.source,self.root/'adapter.f90'],folder)
            run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',self.root,'-c',client],folder)
            exe = folder/'client.exe'
            run(['gfortran','vectors.o','adapter.o','client.o','-o',exe],folder)
            run([exe],folder)
        if shutil.which('g++'):
            run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',self.root,'-fsyntax-only',client],self.root)
