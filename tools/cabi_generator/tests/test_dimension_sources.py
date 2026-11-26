# SPDX-License-Identifier: MIT
"""Reviewed ordered optional-array sources, with independent compiled fixtures."""
import copy
import shutil
import unittest
import test_dimension_defaults as fixtures
from test_generator import run
import generate


class DimensionSourceTests(fixtures.DimensionDefaultTests):
    test_metadata_and_invalid_fallback_policies = None
    test_compiled_runtime_omitted_and_supplied_dimensions = None

    def setUp(self):
        super().setUp()
        self.group = self.policy['array_groups']['vectors']
        self.group['dimension_sources'] = dict(vector_score=dict(n=[
            dict(argument='mean',axis=None),dict(argument='cov',axis=1)]))

    def test_source_metadata_and_invalid_policies(self):
        text=self.emit()[self.root/'python/example/_generated.py']
        self.assertIn("'dimension_sources': {'n':",text)
        original=copy.deepcopy(self.group['dimension_sources'])
        for rules in (None, [], {'unknown':{}}, {'vector_score':{}},
                      {'vector_score':{'unknown':[]}},
                      {'vector_score':{'n':[]}},
                      {'vector_score':{'n':[dict(argument='mean',axis=None)]}}):
            self.group['dimension_sources']=rules
            with self.assertRaises(ValueError):self.emit()
        for second in (dict(argument='x',axis=None),dict(argument='mean',axis=None),
                       dict(argument='cov',axis=True),dict(argument='cov',axis=3),
                       dict(argument=[],axis=1),dict(argument='cov',axis=None),{}):
            self.group['dimension_sources']=copy.deepcopy(original)
            self.group['dimension_sources']['vector_score']['n'][1]=second
            with self.assertRaises(ValueError):self.emit()
        self.group['dimension_sources']=copy.deepcopy(original)
        self.group['dimension_sources']['vector_score']['n'].reverse()
        with self.assertRaises(ValueError):self.emit()

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('python'),'compiler/Python required')
    def test_compiled_source_precedence_and_shape_guards(self):
        outputs=self.emit();outputs.update(generate.generate(self.config,self.root))
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
        script=self.root/'check.py'
        script.write_text('''import ctypes as C
import sys
import numpy as np
sys.path.insert(0,sys.argv[1])
from example._generated import Vectors
api=Vectors(C.CDLL(sys.argv[2]))
assert api.score([1.,2.])==5.
assert api.score([1.,2.,3.],None,np.eye(3))==9.
assert api.score([1.,2.,3.],[4.,5.,6.])==24.
for args in (([1.,2.,3.],),([1.,2.],[],np.eye(2)),
             ([1.,2.],[0.,0.],np.eye(3)),([1.,2.],None,np.ones((2,3)))):
    try:api.score(*args)
    except ValueError:pass
    else:raise AssertionError(args)
''',encoding='ascii')
        for index,flags in enumerate((['-O0','-g','-fcheck=all'],['-O3'])):
            folder=self.root/str(index);folder.mkdir()
            library=folder/('example.dll' if __import__('os').name=='nt' else 'example.so')
            run(['gfortran',*flags,'-shared','-fPIC',self.source,self.root/'adapter.f90','-o',library],folder)
            run(['python',script,self.root/'python',library],folder)
