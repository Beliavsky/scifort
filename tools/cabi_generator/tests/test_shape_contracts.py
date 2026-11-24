# SPDX-License-Identifier: MIT
"""Independent source-contract parsing, manifest agreement and Markdown drift checks."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from test_generator import HERE, array_fixture_config, run
import generate

ANNOTATIONS = """        !! @dimension nrows = size(a,1)
        !! @dimension ncols = size(a,2)
        !! @require nrows >= 1
        !! @require ncols >= 1
        !! @when ncols == 1 : A single column is processed.
"""


class ShapeContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root/'array_ops.f90'
        text = (HERE/'fixtures/array_ops.f90').read_text()
        self.original = text.replace('        !! Scale row sums;', ANNOTATIONS+'        !! Scale row sums;')
        self.source.write_text(self.original)
        self.config = array_fixture_config()
        self.config['subroutines']['row_sums']['shape_contract'] = True
        self.config['outputs']['shape_documentation'] = 'shapes.md'

    def test_deterministic_contract_and_markdown(self):
        before = self.source.read_bytes()
        outputs = generate.generate(self.config,self.root)
        self.assertEqual(outputs,generate.generate(self.config,self.root))
        self.assertEqual(before,self.source.read_bytes())
        report = json.loads(outputs[self.root/'array_coverage.json'])['procedures'][0]
        contract = report['shape_contract']
        self.assertEqual(contract['definitions']['nrows'],dict(argument='a',axis=1))
        self.assertEqual(contract['minimums'],dict(nrows=1,ncols=1))
        self.assertEqual(contract['conditions'][0]['value'],1)
        doc = outputs[self.root/'shapes.md']
        self.assertIn('| `a` | `(:,:)` | `(nrows,ncols)` |',doc)
        self.assertIn('| `y` | `(:)` | `(nrows)` |',doc)
        self.assertIn('When `ncols == 1`:',doc)
        self.assertNotIn('scifort',doc.lower())
        # Reject signed axis typos and unsupported expressions; never eval source text.
        for change in ('@dimension nrows = size(a,1)+1', '@dimension nrows = size(a,2)',
                       '@dimension nrows = size(a)', '@dimension nrows = size(gain,1)',
                       '@dimension nrows = size(a,3)', '@dimension nrows = unknown(a)',
                       '@dimension nrows = __import__(os)', '@dimension other = size(a,1)'):
            self.source.write_text(self.original.replace('@dimension nrows = size(a,1)',change))
            with self.assertRaises(ValueError):
                generate.generate(self.config,self.root)

    def test_incomplete_and_inconsistent_contracts_are_rejected(self):
        for original,replacement in (('@require nrows >= 1','@require nrows >= 0'),
                                     ('@require nrows >= 1','@require unknown >= 1'),
                                     ('@require nrows >= 1','@require nrows > 0'),
                                     ('@require nrows >= 1','plain text'),
                                     ('@when ncols == 1','@when ncols == 0'),
                                     ('@when ncols == 1','@when unknown == 1'),
                                     ('shape (nrows,ncols)','shape (ncols,nrows)')):
            with self.subTest(replacement=replacement):
                self.source.write_text(self.original.replace(original,replacement))
                with self.assertRaises(ValueError):
                    generate.generate(self.config,self.root)
        for duplicate in ('@require nrows >= 1','@dimension nrows = size(a,1)',
                          '@when ncols == 1 : A single column is processed.'):
            self.source.write_text(self.original.replace('        !! @dimension nrows',
                '        !! '+duplicate+'\n        !! @dimension nrows',1))
            with self.assertRaises(ValueError):
                generate.generate(self.config,self.root)
        self.source.write_text(self.original)
        self.config['subroutines']['row_sums']['shape_contract'] = 'true'
        with self.assertRaises(ValueError):
            generate.generate(self.config,self.root)

    def test_cli_checks_docs_and_preserves_outputs_on_invalid_source(self):
        config = self.root/'config.json'
        config.write_text(json.dumps(self.config))
        command = [sys.executable,HERE.parent/'generate.py',config,'--root',self.root]
        run(command,self.root)
        run([*command,'--check'],self.root)
        doc = self.root/'shapes.md'
        doc.write_text('sentinel')
        result = subprocess.run([str(a) for a in [*command,'--check']],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('stale',result.stderr)
        self.assertEqual(doc.read_text(),'sentinel')
        self.source.write_text(self.original.replace('@require nrows >= 1','@require nrows >= 2'))
        result = subprocess.run([str(a) for a in command],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertEqual(doc.read_text(),'sentinel')
        # Contract output cannot overwrite a native source, even through another suffix.
        self.source.write_text(self.original)
        self.config['outputs']['shape_documentation'] = 'array_ops.f90'
        with self.assertRaises(ValueError):
            generate.generate(self.config,self.root)
