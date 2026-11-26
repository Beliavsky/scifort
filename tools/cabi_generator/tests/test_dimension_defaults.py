# SPDX-License-Identifier: MIT
"""Explicit optional-array dimension fallbacks, independent of SciFort."""
import copy
import json
import shutil
import unittest
import test_array_functions as fixtures
from test_generator import HERE, fixture_config, run
import generate
import generate_python


class DimensionDefaultTests(fixtures.ArrayFunctionTests):
    # Run only these tests here; the base module owns its original tests.
    test_metadata_and_native_call = None
    test_rejects_unsupported_and_conflicting_policies = None
    test_python_defaults_prefix_and_validation = None
    test_compiled_checked_optimized_and_c_header = None

    def setUp(self):
        super().setUp()
        self.source.write_text(fixtures.SOURCE.replace('n = size(x)', 'n = size(mean)'),encoding='ascii')
        (self.root/'arithmetic.f90').write_bytes((HERE/'fixtures/arithmetic.f90').read_bytes())
        scalar=fixture_config();scalar['outputs']['coverage']='scalar_coverage.json'
        (self.root/'scalar.json').write_text(json.dumps(scalar))
        self.policy=dict(schema_version=1,output_directory='python/example',
            groups=dict(arithmetic=dict(manifest='scalar.json',prefix='',defaults=dict(shift=0.,gain=1.),existing_functions={})),
            array_groups=dict(vectors=dict(manifests=['arrays.json'],prefix='vector_',
                defaults=dict(vector_score=dict(mean={'fill':'zeros'},cov={'fill':'identity'},enabled=True)),
                dimension_defaults=dict(vector_score=dict(n=2)))))

    def emit(self):
        (self.root/'arrays.json').write_text(json.dumps(self.config))
        path=self.root/'python.json';path.write_text(json.dumps(self.policy))
        return generate_python.artifacts(path,self.root)

    def test_metadata_and_invalid_fallback_policies(self):
        outputs=self.emit();text=outputs[self.root/'python/example/_generated.py']
        self.assertIn("'dimension_defaults': {'n': 2}",text)
        baseline=copy.deepcopy(self.policy)
        for value in ({},[],{'unknown':2},{'n':0},{'n':True},{'n':2.},{'n':2147483648}):
            self.policy=copy.deepcopy(baseline)
            self.policy['array_groups']['vectors']['dimension_defaults']['vector_score']=value
            with self.assertRaises(ValueError):self.emit()
        self.policy=copy.deepcopy(baseline)
        self.policy['array_groups']['vectors']['dimension_defaults']['unknown']={'n':2}
        with self.assertRaises(ValueError):self.emit()
        self.policy=copy.deepcopy(baseline)
        self.source.write_text(fixtures.SOURCE,encoding='ascii')
        with self.assertRaises(ValueError):self.emit()

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('python'),'compiler and Python required')
    def test_compiled_runtime_omitted_and_supplied_dimensions(self):
        outputs=self.emit()
        outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        script=self.root/'check.py'
        script.write_text("""import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Vectors
api=Vectors(C.CDLL(sys.argv[2]))
assert api.score([1.,2.])==5.
assert api.score([1.,2.,3.],[4.,5.,6.])==24.
assert api.score([1.,2.],None,[[2.,0.],[0.,2.]])==7.
for args in (([1.,2.,3.],),([1.,2.],[],None),([1.,2.],None,np.eye(3))):
    try:api.score(*args)
    except ValueError:pass
    else:raise AssertionError(args)
""",encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            library=folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-shared','-fPIC',self.source,self.root/'adapter.f90','-o',library],folder)
            run(['python',script,self.root/'python',library],folder)
