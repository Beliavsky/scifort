# SPDX-License-Identifier: MIT
"""Project-independent generation of Python array groups from structured contracts."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from test_generator import HERE, fixture_config, array_fixture_config
import generate_python
import pipeline


class PythonArrayGeneratorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'arithmetic.f90').write_bytes((HERE/'fixtures/arithmetic.f90').read_bytes())
        (self.root/'scalar.json').write_text(json.dumps(fixture_config()))
        source = (HERE/'fixtures/array_ops.f90').read_text()
        source = source.replace('        !! Scale row sums;', """        !! @dimension nrows = size(a,1)
        !! @dimension ncols = size(a,2)
        !! @require nrows >= 1
        !! @require ncols >= 1
        !! Scale row sums;""")
        (self.root/'array_ops.f90').write_text(source)
        self.array = array_fixture_config()
        self.array['subroutines']['row_sums']['shape_contract'] = True
        self.policy = dict(schema_version=1,output_directory='python/example',groups=dict(
            arithmetic=dict(manifest='scalar.json',prefix='',defaults=dict(shift=0.,gain=1.),existing_functions={})),
            array_groups=dict(matrices=dict(manifests=['arrays.json'],defaults=dict(row_sums=dict(gain=1.)))))
        self.path = self.root/'python.json'

    def generate(self):
        (self.root/'arrays.json').write_text(json.dumps(self.array))
        self.path.write_text(json.dumps(self.policy))
        return generate_python.artifacts(self.path,self.root)

    def test_generic_group_signatures_and_metadata(self):
        outputs = self.generate()
        self.assertEqual(outputs,self.generate())
        text = outputs[self.root/'python/example/_generated.py']
        self.assertIn('class Matrices(ArrayGroup):',text)
        self.assertIn('def row_sums(self, a, gain=1.0, *, out=None):',text)
        self.assertIn("'dimensions': ['nrows', 'ncols']",text)
        self.assertIn('self.matrices = Matrices(self._library)',text)
        self.assertNotIn('scifort',text.lower())
        compile(text,'<generated>','exec')
        compile(outputs[self.root/'python/example/_runtime.py'],'<runtime>','exec')
        # Array adapters use the primary C source, with no duplicate compatibility module.
        self.assertNotIn(self.root/'python/native/matrices_compat.f90',outputs)

    def test_missing_contract_bad_defaults_and_dimension_sources(self):
        self.array['subroutines']['row_sums'].pop('shape_contract')
        with self.assertRaisesRegex(ValueError,'structured shape contract'):
            self.generate()
        self.array['subroutines']['row_sums']['shape_contract'] = True
        original = copy.deepcopy(self.policy)
        for defaults in (dict(unknown={}),dict(row_sums=dict(a=1)),
                         dict(row_sums=dict(gain='1')),dict(row_sums=dict(gain=float('nan')))):
            self.policy = copy.deepcopy(original)
            self.policy['array_groups']['matrices']['defaults'] = defaults
            with self.assertRaises(ValueError):
                self.generate()
        self.policy = original
        source = self.root/'array_ops.f90'
        source.write_text(source.read_text().replace('nrows = size(a,1)','nrows = size(y)'))
        with self.assertRaisesRegex(ValueError,'input array'):
            self.generate()

    def test_group_collisions_and_invalid_shapes_of_policy(self):
        original = copy.deepcopy(self.policy)
        for group in ('arithmetic','path','library'):
            self.policy = copy.deepcopy(original)
            self.policy['array_groups'][group] = self.policy['array_groups'].pop('matrices')
            with self.assertRaises(ValueError):
                self.generate()
        self.policy = original
        self.policy['array_groups']['matrices']['manifests'] *= 2
        with self.assertRaises(ValueError):
            self.generate()

    def test_pipeline_build_and_coverage_include_array_groups(self):
        self.generate()
        (self.root/'batch.json').write_text('["scalar.json", "arrays.json"]')
        config = dict(schema_version=1,c_manifests='batch.json',python_policy='python.json',
            source_directories=['.'],build_manifest='plan.json',coverage_report='interfaces.json',
            aggregate_header='include/all.h',aggregate_header_guard='ALL_H')
        path = self.root/'pipeline.json'
        path.write_text(json.dumps(config))
        outputs = pipeline.artifacts(path,self.root)
        plan = json.loads(outputs[self.root/'plan.json'])
        self.assertEqual(plan['families'],['arithmetic'])
        self.assertEqual(plan['array_groups'],['matrices'])
        self.assertLess(plan['sources'].index('array_ops.f90'),plan['sources'].index('array_api.f90'))
        report = json.loads(outputs[self.root/'interfaces.json'])
        row = next(r for r in report['procedures'] if r['procedure']=='row_sums')
        self.assertEqual(row['python'],'matrices.row_sums')
        (self.root/'batch.json').write_text('["scalar.json"]')
        with self.assertRaisesRegex(ValueError,'C batch'):
            pipeline.artifacts(path,self.root)
