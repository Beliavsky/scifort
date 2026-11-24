# SPDX-License-Identifier: MIT
"""Project-independent Python generator policy and drift checks."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from test_generator import fixture_config, HERE
sys.path.insert(0, str(HERE.parent))
import generate_python


class PythonGeneratorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'arithmetic.f90').write_bytes((HERE / 'fixtures/arithmetic.f90').read_bytes())
        (self.root / 'manifest.json').write_text(json.dumps(fixture_config()))
        self.config = dict(schema_version=1, output_directory='python/example',
                           groups=dict(arithmetic=dict(manifest='manifest.json', prefix='',
                                      defaults=dict(shift=0., gain=1.), existing_functions={})))
        self.path = self.root / 'python.json'

    def generate(self):
        self.path.write_text(json.dumps(self.config))
        return generate_python.artifacts(self.path, self.root)

    def test_deterministic_generic_signatures_and_no_project_dependency(self):
        outputs = self.generate()
        self.assertEqual(outputs, self.generate())
        generated = outputs[self.root / 'python/example/_generated.py']
        self.assertIn('def affine(self, t, shift=0.0, gain=1.0, *, out=None):', generated)
        self.assertNotIn('scifort', generated.lower())
        compile(generated, '<generated>', 'exec')
        runtime = outputs[self.root / 'python/example/_runtime.py']
        self.assertNotIn('fparser', runtime)
        compile(runtime, '<runtime>', 'exec')

    def test_invalid_defaults_and_unaccounted_existing_exports(self):
        original = copy.deepcopy(self.config)
        for defaults in (dict(unknown=1.), dict(shift=float('nan')), dict(shift=0.)):
            self.config = copy.deepcopy(original)
            self.config['groups']['arithmetic']['defaults'] = defaults
            with self.assertRaises(ValueError):
                self.generate()
        self.config = copy.deepcopy(original)
        self.config['groups']['arithmetic']['existing_functions'] = {'affine': {}}
        with self.assertRaisesRegex(ValueError, 'excluded'):
            self.generate()

    def test_elementwise_external_scalar_and_bulk_mapping(self):
        manifest = fixture_config()
        original = copy.deepcopy(manifest['functions']['affine'])
        manifest['functions']['affine'] = dict(existing_scalar_symbol=original['symbol'],
            elementwise_symbol='example_affine_elements', optional_policy='required')
        (self.root / 'manifest.json').write_text(json.dumps(manifest))
        extras = self.config['groups']['arithmetic']['existing_functions']
        extras['affine'] = original
        outputs = self.generate()
        generated = outputs[self.root / 'python/example/_generated.py']
        self.assertIn("'elementwise': 'example_affine_elements'", generated)
        native = outputs[self.root / 'python/native/arithmetic_compat.f90']
        self.assertNotIn('example_affine_elements', native)
        self.assertIn('function example_affine(', native)
        self.assertIn('subroutine example_affine_array(', native)
        extras['affine']['symbol'] = 'wrong_scalar'
        with self.assertRaisesRegex(ValueError, 'must match'):
            self.generate()
        extras['affine'] = dict(original, elementwise_symbol='wrong_extra')
        with self.assertRaisesRegex(ValueError, 'primary C manifest'):
            self.generate()

    def test_empty_and_reserved_groups_are_rejected(self):
        original = copy.deepcopy(self.config)
        self.config['groups'] = {}
        with self.assertRaisesRegex(ValueError, 'nonempty'):
            self.generate()
        self.config = copy.deepcopy(original)
        self.config['groups']['library'] = self.config['groups'].pop('arithmetic')
        with self.assertRaisesRegex(ValueError, 'reserved'):
            self.generate()

    def test_check_and_invalid_policy_never_rewrite_outputs(self):
        self.generate()
        command = [sys.executable, HERE.parent / 'generate_python.py', self.path, '--root', self.root]
        subprocess.run(command, check=True, capture_output=True)
        subprocess.run([*command, '--check'], check=True, capture_output=True)
        output = self.root / 'python/example/_generated.py'
        output.write_text('sentinel')
        result = subprocess.run([*command, '--check'], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('stale', result.stderr)
        self.assertEqual(output.read_text(), 'sentinel')
        self.config['groups']['arithmetic']['defaults']['gain'] = None
        self.path.write_text(json.dumps(self.config))
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(output.read_text(), 'sentinel')
