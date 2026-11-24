# SPDX-License-Identifier: MIT
"""Project-independent inventory, binding discovery and deliberately unsafe draft checks."""
import json
import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from test_generator import HERE, fixture_config
sys.path.insert(0, str(HERE.parent))
import discover
import generate


class DiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'arithmetic.f90'
        self.source.write_bytes((HERE / 'fixtures/arithmetic.f90').read_bytes())

    def test_types_ranks_optionals_and_explicit_unsupported_reasons(self):
        report = discover.inventory([self.source], self.root, kind='c_double')
        rows = {p['name']: p for p in report['modules'][0]['procedures']}
        self.assertTrue(rows['affine']['scalar_supported'])
        self.assertTrue(rows['affine']['bulk_supported'])
        self.assertEqual(rows['affine']['optional_arguments'], ['shift', 'gain'])
        self.assertIn('c_double', rows['affine']['arguments'][0]['type'])
        self.assertEqual(rows['affine']['arguments'][0]['rank'], 0)
        self.assertFalse(rows['total']['scalar_supported'])
        self.assertEqual(rows['total']['arguments'][0]['rank'], 1)
        self.assertIn('scalar intent', rows['total']['unsupported_reason'])
        self.assertFalse(report['incomplete'])
        self.assertEqual(report, discover.inventory([self.source], self.root, kind='c_double'))

    def test_observed_c_symbols_follow_explicit_import_renames(self):
        config = fixture_config()
        outputs = generate.generate(config, self.root)
        for path, text in outputs.items():
            path.write_text(text)
        report = discover.inventory([self.source], self.root, 'c_double', [self.root/'api.f90'])
        affine = report['modules'][0]['procedures'][0]
        self.assertEqual({s['symbol'] for s in affine['observed_c_symbols']},
                         {'example_affine', 'example_affine_array'})
        self.assertEqual(report['coverage']['observed_c_scalar'], 1)
        self.assertEqual(report['coverage']['observed_c_procedural'], 1)
        self.assertEqual(report['coverage']['configured_scalar'], 0)

    def test_nonexecutable_draft_contains_no_inferred_constraints_or_defaults(self):
        report = discover.inventory([self.source], self.root, 'c_double')
        draft = discover.draft(report['modules'][0], 'c_double', 'iso_c_binding', 'example_')
        self.assertTrue(draft['draft'])
        policy = draft['functions']['affine']
        self.assertEqual(policy['finite'], 'REVIEW_REQUIRED')
        self.assertEqual(policy['positive'], 'REVIEW_REQUIRED')
        self.assertIsNone(policy['array_argument'])
        self.assertNotIn('defaults', draft)
        self.assertIn('total', draft['excluded'])
        with self.assertRaisesRegex(ValueError, 'schema versions'):
            generate.generate(draft, self.root)

    def test_public_types_constants_and_generics_are_not_silently_discarded(self):
        self.source.write_text(self.source.read_text().replace('contains\n', '''integer, public, parameter :: count = 3
    type, public :: box
        real(c_double) :: value
    end type box
    interface evaluate
        module procedure affine
    end interface evaluate
    public :: evaluate
contains
''', 1))
        report = discover.inventory([self.source], self.root, 'c_double')
        rows = {p['name']: p for p in report['modules'][0]['procedures']}
        self.assertEqual(report['coverage']['public_exports'], 5)
        for name in ('box', 'count', 'evaluate'):
            self.assertEqual(rows[name]['kind'], 'unresolved_public_export')
            self.assertIn('explicit resolution', rows[name]['unsupported_reason'])

    def test_dimension_attributes_and_assumed_size_ranks(self):
        self.source.write_text("""module shapes
implicit none
private
public :: run
contains
subroutine run(x,y,z)
implicit none
real, intent(in) :: x(*)
real, dimension(2,*), intent(in) :: y
real, intent(in) :: z(2,3)
end subroutine
end module
""")
        report = discover.inventory([self.source], self.root)
        row = report['modules'][0]['procedures'][0]
        self.assertEqual([a['rank'] for a in row['arguments']], [1, 2, 2])
        self.assertFalse(row['scalar_supported'])

    def test_default_public_module_is_reported_incomplete(self):
        self.source.write_text(self.source.read_text().replace('private', 'public', 1))
        report = discover.inventory([self.source], self.root, 'c_double')
        self.assertTrue(report['incomplete'])
        self.assertIn('private-by-default', report['modules'][0]['error'])

    def test_batch_and_python_policy_coverage_are_separate(self):
        (self.root/'config.json').write_text(json.dumps(fixture_config()))
        batch = self.root/'batch.json'
        batch.write_text('["config.json"]')
        python = self.root/'python.json'
        python.write_text(json.dumps(dict(schema_version=1, output_directory='python/example',
                           groups=dict(arithmetic=dict(manifest='config.json', prefix='', defaults={},
                                                        existing_functions={})))))
        policies = discover.configured(self.root, batch, python)
        report = discover.inventory([self.source], self.root, 'c_double', policies=policies)
        self.assertEqual(report['coverage']['configured_scalar'], 1)
        self.assertEqual(report['coverage']['configured_bulk'], 1)
        self.assertEqual(report['coverage']['python_methods'], 1)

    def test_explicit_size_limit_reports_skipped_source_as_incomplete(self):
        progress = io.StringIO()
        with contextlib.redirect_stderr(progress):
            report = discover.inventory([self.source], self.root, 'c_double', max_source_bytes=1, progress=True)
        self.assertIn('Inspecting 1/1: arithmetic.f90', progress.getvalue())
        self.assertTrue(report['incomplete'])
        self.assertEqual(report['coverage']['public_exports'], 0)
        self.assertIn('max-source-bytes', report['errors'][0]['reason'])
        text = discover.markdown(discover.inventory([self.source], self.root, 'c_double'))
        table = []
        for line in (text+'\n\n').splitlines():
            if line.startswith('|'):
                table.append(line.count('|'))
            elif table:
                self.assertEqual(len(set(table)), 1)
                table = []

    def test_cli_protects_individual_manifests_referenced_by_batch(self):
        config = self.root/'config.json'
        config.write_text(json.dumps(fixture_config()))
        batch = self.root/'batch.json'
        batch.write_text('["config.json"]')
        original = config.read_bytes()
        command = [sys.executable, HERE.parent/'discover.py', self.source, '--root', self.root,
                   '--kind', 'c_double', '--manifest-batch', batch, '--output', config]
        result = subprocess.run(command, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(config.read_bytes(), original)

    def test_cli_never_overwrites_sources_or_existing_drafts(self):
        command = [sys.executable, HERE.parent/'discover.py', self.source, '--root', self.root,
                   '--kind', 'c_double', '--output', self.root/'report.json',
                   '--draft-directory', self.root/'drafts', '--symbol-prefix', 'example_',
                   '--kind-module', 'iso_c_binding']
        subprocess.run(command, check=True, capture_output=True)
        draft = self.root/'drafts/example_arithmetic.json'
        draft.write_text('sentinel')
        result = subprocess.run(command, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(draft.read_text(), 'sentinel')
        original = self.source.read_bytes()
        command = command[:command.index('--draft-directory')]
        command[command.index('--output')+1] = self.source
        result = subprocess.run(command, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.source.read_bytes(), original)
