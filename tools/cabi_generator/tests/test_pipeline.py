# SPDX-License-Identifier: MIT
"""Independent dependency-closure, regeneration and non-writing workflow tests."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from test_generator import HERE, fixture_config
sys.path.insert(0, str(HERE.parent))
import pipeline


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'source').mkdir()
        native = (HERE/'fixtures/arithmetic.f90').read_text()
        # A dependency used inside a procedure must enter the closure too.
        native = native.replace('module example_arithmetic',
                                'module example_arithmetic\n    use numeric_support, only: identity',1)
        (self.root/'source/unrelated_filename.f90').write_text(native)
        (self.root/'source/dependency.f90').write_text('''module numeric_support
implicit none
private
public :: identity
contains
pure function identity(x) result(y)
implicit none
real,intent(in) :: x
real :: y
y=x
end function
end module
''')
        config = fixture_config()
        config['source']='source/unrelated_filename.f90'
        config['outputs']=dict(fortran='generated/api.f90',header='include/api.h',coverage='coverage.json')
        (self.root/'manifest.json').write_text(json.dumps(config))
        (self.root/'batch.json').write_text('["manifest.json"]')
        (self.root/'python.json').write_text(json.dumps(dict(schema_version=1,
            output_directory='python/example',groups=dict(arithmetic=dict(manifest='manifest.json',
                prefix='',defaults=dict(shift=0.,gain=1.),existing_functions={})))))
        self.config = dict(schema_version=1,c_manifests='batch.json',python_policy='python.json',
            source_directories=['source'],build_manifest='plan.json',coverage_report='interfaces.json',
            aggregate_header='include/all.h',aggregate_header_guard='EXAMPLE_ALL_H')
        self.path=self.root/'pipeline.json'
        self.path.write_text(json.dumps(self.config))

    def collect(self):
        return pipeline.artifacts(self.path,self.root)

    def test_independent_ordered_closure_and_all_outputs_without_writes(self):
        outputs=self.collect()
        self.assertEqual(outputs,self.collect())
        self.assertFalse((self.root/'generated/api.f90').exists())
        plan=json.loads(outputs[self.root/'plan.json'])
        self.assertEqual(plan['sources'][:2],['source/dependency.f90','source/unrelated_filename.f90'])
        self.assertEqual(plan['modules']['example_arithmetic'],'source/unrelated_filename.f90')
        self.assertEqual(plan['families'],['arithmetic'])
        self.assertNotIn('scifort',json.dumps(plan).lower())
        self.assertIn('#include "api.h"',outputs[self.root/'include/all.h'])
        coverage=json.loads(outputs[self.root/'interfaces.json'])
        self.assertEqual(coverage['counts'],dict(exports=2,scalar=1,bulk=1,array_adapters=0,elementwise=0,python=1))
        # No compiler-specific intrinsics need fake source files.
        self.assertNotIn('iso_c_binding',plan['modules'])
        fc=shutil.which('gfortran')
        if fc:
            for path,text in outputs.items():
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(text,encoding='utf-8')
            for flags in (['-O0','-fcheck=all'],['-O3']):
                result=subprocess.run([fc,*flags,'-std=f2018','-fimplicit-none',
                    '-Werror=implicit-interface','-shared','-fPIC',*plan['sources'],'-o','fixture.dll' if sys.platform=='win32' else 'fixture.so'],
                    cwd=self.root,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_missing_ambiguous_and_cyclic_dependencies(self):
        dep=self.root/'source/dependency.f90'
        content=dep.read_text()
        dep.unlink()
        with self.assertRaisesRegex(ValueError,'missing provider for module numeric_support'):
            self.collect()
        dep.write_text(content)
        duplicate=self.root/'source/duplicate.f90'
        duplicate.write_text(content)
        with self.assertRaisesRegex(ValueError,'ambiguous provider for module numeric_support'):
            self.collect()
        duplicate.unlink()
        dep.write_text(content.replace('module numeric_support','module numeric_support\nuse example_arithmetic',1))
        with self.assertRaisesRegex(ValueError,'dependency cycle'):
            self.collect()

    def test_procedure_local_use_and_same_file_module_order(self):
        native=self.root/'source/unrelated_filename.f90'
        text=native.read_text().replace('    use numeric_support, only: identity\n','')
        text=text.replace('pure elemental function affine(t, shift, gain) result(y)',
                          'pure elemental function affine(t, shift, gain) result(y)\nuse numeric_support, only: identity')
        native.write_text(text)
        self.assertIn('source/dependency.f90',json.loads(self.collect()[self.root/'plan.json'])['sources'])
        first='module earlier\nimplicit none\nend module earlier\n'
        second='module later\nuse earlier\nimplicit none\nend module later\n'
        path=self.root/'source/multiple.f90'
        path.write_text(first+second)
        plan,_=pipeline.build_plan([path],['source'],self.root,{})
        self.assertEqual(plan['sources'],['source/multiple.f90'])
        path.write_text(second+first)
        with self.assertRaisesRegex(ValueError,'must precede'):
            pipeline.build_plan([path],['source'],self.root,{})

    def test_intrinsic_resolution_and_unsupported_inputs(self):
        path=self.root/'source/intrinsic.f90'
        path.write_text('module intrinsic_example\nuse, non_intrinsic :: iso_c_binding\nimplicit none\nend module\n')
        with self.assertRaisesRegex(ValueError,'missing provider for module iso_c_binding'):
            pipeline.build_plan([path],['source'],self.root,{})
        path.write_text('module intrinsic_example\nuse, intrinsic :: iso_c_binding\nimplicit none\nend module\n')
        self.assertEqual(pipeline.build_plan([path],['source'],self.root,{})[0]['sources'],['source/intrinsic.f90'])
        path.write_text('#define EXAMPLE 1\nmodule intrinsic_example\nend module\n')
        with self.assertRaisesRegex(ValueError,'preprocessing'):
            pipeline.build_plan([path],['source'],self.root,{})

    def test_cli_checks_missing_stale_plan_python_header_and_never_rewrites_on_error(self):
        command=[sys.executable,HERE.parent/'pipeline.py',self.path,'--root',self.root]
        subprocess.run(command,check=True,capture_output=True)
        subprocess.run([*command,'--check'],check=True,capture_output=True)
        for name in ('plan.json','interfaces.md','include/all.h','python/example/_generated.py','generated/api.f90'):
            path=self.root/name
            original=path.read_text()
            path.write_text('sentinel')
            result=subprocess.run([*command,'--check'],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('stale',result.stderr)
            self.assertEqual(path.read_text(),'sentinel')
            path.write_text(original)
        plan=self.root/'plan.json'
        plan.write_text('sentinel')
        (self.root/'source/dependency.f90').unlink()
        result=subprocess.run(command,capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertEqual(plan.read_text(),'sentinel')

    def test_manual_header_and_unlisted_python_manifest_are_rejected(self):
        header=self.root/'include/all.h'
        header.parent.mkdir()
        header.write_text('/* manually maintained */')
        with self.assertRaisesRegex(ValueError,'manual header'):
            self.collect()
        header.unlink()
        other=json.loads((self.root/'manifest.json').read_text())
        (self.root/'other.json').write_text(json.dumps(other))
        python=json.loads((self.root/'python.json').read_text())
        python['groups']['arithmetic']['manifest']='other.json'
        (self.root/'python.json').write_text(json.dumps(python))
        with self.assertRaisesRegex(ValueError,'same batch|C batch'):
            self.collect()

    def test_pipeline_outputs_cannot_overwrite_original_sources_or_inputs(self):
        self.config['coverage_report']='source/dependency.f90'
        self.path.write_text(json.dumps(self.config))
        with self.assertRaisesRegex(ValueError,'overwrite'):
            self.collect()
        self.config['coverage_report']='interfaces.json'
        self.config['build_manifest']=str(Path(pipeline.__file__).resolve())
        self.path.write_text(json.dumps(self.config))
        with self.assertRaisesRegex(ValueError,'overwrite'):
            self.collect()
