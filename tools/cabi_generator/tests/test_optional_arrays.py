# SPDX-License-Identifier: MIT
"""Presence-controlled arrays preserve PRESENT semantics without C descriptors."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from test_generator import HERE, fixture_config, run
import generate
import generate_python

SOURCE='''module example_optional
use iso_c_binding, only: c_double
implicit none
private
public :: summarize
contains
pure subroutine summarize(x,total,status,mean,cov)
    !! @dimension n = size(x)
    !! @require n >= 1
    implicit none
    real(c_double), intent(in) :: x(:) !! shape (n); observations
    real(c_double), intent(out) :: total !! presence-sensitive sum
    integer, intent(out) :: status !! zero on success
    real(c_double), intent(in), optional :: mean(:) !! shape (n); optional vector
    real(c_double), intent(in), optional :: cov(:,:) !! shape (n,n); optional matrix
    total=sum(x)
    if (present(mean)) total=total+100+sum(mean)
    if (present(cov)) total=total+200+sum(cov)
    status=0
end subroutine
end module
'''


class OptionalArrayTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.source=self.root/'optional.f90'
        self.source.write_text(SOURCE,encoding='ascii')
        self.config=dict(schema_version=2,source='optional.f90',module='example_optional',
            kind_module='iso_c_binding',kind_name='c_double',wrapper_module='example_optional_api',
            header_guard='EXAMPLE_OPTIONAL_H',functions={},excluded={},
            outputs=dict(fortran='adapter.f90',header='adapter.h',coverage='coverage.json'),
            argument_docs=dict(n='dimension',x='observations',total='sum',status='status',
                               mean='optional vector',cov='optional matrix',use_mean='presence flag',use_cov='presence flag'),
            subroutines=dict(summarize=dict(symbol='example_summary',
                dimensions=[dict(name='n',minimum=1)],shapes=dict(x=['n'],mean=['n'],cov=['n','n']),
                status_argument='status',invalid_status=1,description='Summarize optional arrays.',
                shape_contract=True,optional_arrays=dict(mean='use_mean',cov='use_cov'))))
        (self.root/'arithmetic.f90').write_bytes((HERE/'fixtures/arithmetic.f90').read_bytes())
        scalar=fixture_config();scalar['outputs']['coverage']='scalar_coverage.json'
        (self.root/'scalar.json').write_text(json.dumps(scalar))
        self.policy=dict(schema_version=1,output_directory='python/example',
            groups=dict(arithmetic=dict(manifest='scalar.json',prefix='',defaults=dict(shift=0.,gain=1.),existing_functions={})),
            array_groups=dict(optional=dict(manifests=['arrays.json'],defaults={})))

    def python_outputs(self):
        (self.root/'arrays.json').write_text(json.dumps(self.config))
        path=self.root/'python.json';path.write_text(json.dumps(self.policy))
        return generate_python.artifacts(path,self.root)

    def test_metadata_presence_forwarding_and_python_signature(self):
        outputs=generate.generate(self.config,self.root)
        self.assertEqual(outputs,generate.generate(self.config,self.root))
        source=outputs[self.root/'adapter.f90']
        self.assertIn('allocate(cg_optional_mean, source=mean)',source)
        self.assertIn('mean=cg_optional_mean',source)
        self.assertNotIn('use_mean=use_mean',source)
        record=json.loads(outputs[self.root/'coverage.json'])['procedures'][0]
        self.assertEqual(record['optional_arrays'],dict(mean='use_mean',cov='use_cov'))
        self.assertEqual(record['optional_promoted_to_required'],[])
        self.assertEqual(record['arguments'],['x','total','status','use_mean','mean','use_cov','cov'])
        text=self.python_outputs()[self.root/'python/example/_generated.py']
        self.assertIn('def summarize(self, x, mean=None, cov=None, *, out=None):',text)
        self.assertIn("'inputs': ['x', 'mean', 'cov']",text)

    def test_rejects_unsupported_native_and_flag_policies(self):
        baseline=copy.deepcopy(self.config)
        for mapping in (None,[],{},dict(unknown='use_x'),dict(x='use_x'),dict(status='use_status'),
                        dict(mean='x'),dict(mean='n'),dict(mean='cg_hidden'),
                        dict(mean='use_mean',cov='use_mean'),dict(mean=1)):
            self.config=copy.deepcopy(baseline)
            self.config['subroutines']['summarize']['optional_arrays']=mapping
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        self.config=baseline
        for source in (SOURCE.replace('intent(in), optional :: mean','intent(out), optional :: mean'),
                       SOURCE.replace('intent(in), optional :: mean','intent(in) :: mean'),
                       SOURCE.replace('mean(:) !!','mean !!')):
            self.source.write_text(source,encoding='ascii')
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        self.source.write_text(SOURCE,encoding='ascii')
        self.policy['array_groups']['optional']['defaults']={'summarize':{'mean':{'fill':'zeros'}}}
        with self.assertRaises(ValueError):self.python_outputs()
        self.policy['array_groups']['optional']['defaults']={}
        self.source.write_text(SOURCE.replace('n = size(x)','n = size(mean)'),encoding='ascii')
        with self.assertRaises(ValueError):self.python_outputs()

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc') and shutil.which('python'),
                         'Fortran/C compilers and Python required')
    def test_compiled_c_and_python_all_presence_combinations(self):
        outputs=self.python_outputs();outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        client=self.root/'client.c'
        client.write_text('''#include "adapter.h"
#include <math.h>
int main(void) {
 double x[2]={1,2},mean[2]={4,5},cov[4]={1,0,0,1},total; int status;
 for(int m=0;m<=1;m++) for(int c=0;c<=1;c++) {
   example_summary(2,x,&total,&status,m,mean,c,cov);
   if(status || total!=3+m*109+c*202) return 1;
 }
 mean[0]=NAN;cov[0]=NAN;
 example_summary(2,x,&total,&status,0,mean,0,cov);
 if(status || total!=3) return 2;
 total=99; example_summary(2,x,&total,&status,1,mean,2,cov);
 if(status!=1 || total!=99) return 3;
 example_summary(0,x,&total,&status,0,mean,0,cov);
 return status!=1 || total!=99;
}
''',encoding='ascii')
        script=self.root/'check.py'
        script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Optional
api=Optional(C.CDLL(sys.argv[2]))
for mean in (None,[4.,5.]):
 for cov in (None,np.eye(2)):
  result=api.summarize([1.,2.],mean,cov)
  assert result.status==0
  assert result.outputs['total']==3+109*(mean is not None)+202*(cov is not None)
assert np.isnan(api.summarize([1.,2.],[np.nan,1.]).outputs['total'])
for mean,cov in (([],None),([1.],None),(None,np.eye(3)),(None,[1.,2.])):
 try:api.summarize([1.,2.],mean,cov)
 except ValueError:pass
 else:raise AssertionError((mean,cov))
''',encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all','-fbacktrace'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            library=folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-std=f2018','-Werror=implicit-interface','-shared','-fPIC',
                 self.source,self.root/'adapter.f90','-o',library],folder)
            exe=folder/'client.exe'
            run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',self.root,client,library,'-o',exe],folder)
            run([exe],folder)
            run(['python',script,self.root/'python',library],folder)
        if shutil.which('g++'):
            run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',self.root,'-fsyntax-only',client],self.root)
