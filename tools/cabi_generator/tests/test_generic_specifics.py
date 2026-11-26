# SPDX-License-Identifier: MIT
"""Reviewed public generic selection without exposing private specific names."""
import copy
import json
import shutil
import unittest
import test_array_functions as fixtures
import generate


class GenericSpecificTests(unittest.TestCase):
    def setUp(self):
        fixtures.ArrayFunctionTests.setUp(self)
        source=fixtures.SOURCE.replace('contains',
            'interface vector_score\nmodule procedure vector_score_real, vector_score_int\nend interface\ncontains',1)
        source=source.replace('function vector_score(', 'function vector_score_real(')
        integer=fixtures.SOURCE.split('contains\n',1)[1].rsplit('end module',1)[0]
        integer=integer.replace('function vector_score(', 'function vector_score_int(').replace('real(c_double), intent(in) :: x', 'integer, intent(in) :: x')
        self.source.write_text(source.replace('end module',integer+'end module'),encoding='ascii')
        self.config['generic_specifics']={'vector_score':'vector_score_real'}

    def test_selection_validation_and_coverage(self):
        self.config['outputs']['shape_documentation']='shapes.md'
        outputs=generate.generate(self.config,self.root)
        self.assertIn('reviewed specific `vector_score_real`',outputs[self.root/'shapes.md'])
        text=outputs[self.root/'adapter.f90']
        self.assertIn('cg_native_0 => vector_score',text)
        self.assertNotIn('=> vector_score_real',text)
        record=json.loads(outputs[self.root/'coverage.json'])['procedures'][0]
        self.assertEqual(record['native_specific'],'vector_score_real')
        original=copy.deepcopy(self.config)
        for mapping in ([],{'missing':'vector_score_real'},{'vector_score':'missing'},
                        {'vector_score_real':'vector_score_real'},{'vector_score':'vector_score'}):
            self.config=copy.deepcopy(original);self.config['generic_specifics']=mapping
            with self.assertRaises(ValueError):generate.generate(self.config,self.root)
        self.config=copy.deepcopy(original)
        self.config['array_functions']={};self.config['excluded']={'vector_score':'excluded'}
        with self.assertRaises(ValueError):generate.generate(self.config,self.root)

    def test_python_and_pipeline_generic_metadata(self):
        fixtures.ArrayFunctionTests.test_python_defaults_prefix_and_validation(self)

    @unittest.skipUnless(shutil.which('gfortran') and shutil.which('gcc'),'GNU compilers required')
    def test_compiled_private_specific_via_public_generic(self):
        fixtures.ArrayFunctionTests.test_compiled_checked_optimized_and_c_header(self)
