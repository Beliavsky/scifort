# SPDX-License-Identifier: MIT
"""Check the public-Python benchmark path and its Markdown summary tables."""
import json
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_python_scipy as benchmark


class BenchmarkTests(unittest.TestCase):
    def test_all_56_cases_and_public_scipy_methods(self):
        self.assertEqual(len(benchmark.CASES), 56)
        class Normal:
            def cdf(self, values):
                return stats.norm.cdf(values)
        class API:
            normal = Normal()
        fort, python, inputs = benchmark.prepare(API(), 'normal_cdf_array', 16)
        for k in range(2):
            np.testing.assert_array_equal(python(k), stats.norm.cdf(inputs[:, k]))
            np.testing.assert_array_equal(fort(k), python(k))
            self.assertFalse(np.shares_memory(fort(k), fort(k)))

    def test_markdown_summary_column_counts_and_check_only(self):
        row = dict(case='normal_cdf_array', category='distribution', n=128,
                   scifort_median_seconds=.001, scipy_median_seconds=.002,
                   scipy_over_scifort=2., scifort_seconds=[.001]*3, scipy_seconds=[.002]*3)
        report = dict(started_utc='test', quick=True, check_only=False, cases=[row])
        text = benchmark.markdown(report)
        self.assertIn('Geometric mean', text)
        self.assertIn('No direct scipy.special kernels', text)
        table = []
        for line in (text+'\n\n').splitlines():
            if line.startswith('|'):
                table.append(line.count('|'))
            elif table:
                self.assertEqual(len(set(table)), 1)
                table = []
        report['check_only'] = True
        text = benchmark.markdown(report)
        self.assertNotIn('Geometric mean', text)
        self.assertIn('Validated 1 workloads', text)


if __name__ == '__main__':
    unittest.main()
