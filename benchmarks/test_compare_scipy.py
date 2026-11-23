# SPDX-License-Identifier: MIT
"""Focused regression checks for benchmark validation and reporting safeguards."""
import subprocess
import sys
import unittest
from pathlib import Path

import compare_scipy as speed
import numpy as np


class BenchmarkSafeguards(unittest.TestCase):
    def test_small_cdf_tail_cannot_pass_on_absolute_error_alone(self):
        calculate, inputs, exact = speed.prepare("normal_cdf_array", 16)
        actual = np.column_stack([calculate(k).copy() for k in range(2)]).astype(complex)
        actual[0, 0] *= 2
        with self.assertRaises(AssertionError):
            speed.validate(actual, calculate, "normal_cdf_array", inputs, exact)

    def test_fft_imaginary_corruption_is_rejected(self):
        calculate, inputs, exact = speed.prepare("ifft_reuse", 64)
        actual = np.column_stack([calculate(k).copy() for k in range(2)])
        actual[17, 0] += .01j
        with self.assertRaises(RuntimeError):
            speed.validate(actual, calculate, "ifft_reuse", inputs, exact)

    def test_solve_checks_known_solution_independently(self):
        calculate, inputs, exact = speed.prepare("cho_reuse", 16)
        wrong = [calculate(k).copy() + .01 for k in range(2)]
        with self.assertRaises(AssertionError):
            speed.validate(np.column_stack(wrong), lambda k: wrong[k],
                           "cho_reuse", inputs, exact)

    def test_timed_checksum_includes_imaginary_output(self):
        expected = np.array([[1+2j, 3+4j], [5+6j, 7+8j]])
        speed.check_checksum(14+22+14, expected, 3)
        with self.assertRaises(RuntimeError):
            speed.check_checksum(6+10+6, expected, 3)

    def test_case_listing_does_not_require_compiler(self):
        result = subprocess.run([sys.executable, str(Path(speed.__file__)),
                                 "--list", "--category", "fft", "--case", "ifft_prime_setup",
                                 "--compiler", "compiler-that-does-not-exist"],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("ifft_prime_setup", result.stdout)
        self.assertNotIn("fft_power2_reuse", result.stdout)

    def test_empty_selection_is_an_error(self):
        result = subprocess.run([sys.executable, str(Path(speed.__file__)),
                                 "--list", "--category", "fft", "--case", "lu_reuse"],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no workloads match", result.stderr)

    def test_ratio_summaries_and_category_selection(self):
        rows = [dict(case=name, category=group, scifort_seconds=[1.], scipy_over_scifort=value)
                for name, group, value in (("normal_cdf_array", "distribution", .25),
                                           ("normal_ppf_array", "distribution", 1.),
                                           ("fft_prime_reuse", "fft", 4.))]
        text = "\n".join(speed.ratio_summary_lines(rows))
        self.assertIn("| All | 3 | 1.000 | 1.000 | 1.750 | 1.984 | 0.250 | 4.000 |", text)
        self.assertIn("| Distributions | 2 | 0.625 | 0.500 | 0.625 | 0.530 | 0.250 | 1.000 |", text)
        self.assertIn("| FFTs | 1 | 4.000 | 4.000 | 4.000 | n/a | 4.000 | 4.000 |", text)
        table = [line for line in text.splitlines() if line.startswith("|")]
        # GFM requires the delimiter row to have the same cell count as the header.
        self.assertEqual(len(table[0].strip("|").split("|")), 8)
        for line in table[1:]:
            self.assertEqual(len(line.strip("|").split("|")), 8, line)
        self.assertNotIn("| Linear solves |", text)
        self.assertNotIn("| Interpolation |", text)
        # Validation-only cases must not enter the timing statistics.
        rows.append(dict(case="barycentric", category="interpolation"))
        self.assertEqual("\n".join(speed.ratio_summary_lines(rows)), text)
        self.assertEqual(speed.ratio_summary_lines([rows[-1]]), [])

    def test_ratio_summaries_reject_invalid_ratios(self):
        for value in (0., -1., float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                speed.ratio_summary_lines([dict(case="normal_cdf_array", scifort_seconds=[1.],
                                               scipy_over_scifort=value)])

    def test_markdown_preserves_timings_and_check_only_results(self):
        base = dict(started_utc="test", platform="test", build=dict(compiler="gfortran", flags=[]),
                    python="3 test", numpy="test", scipy="test", quick=False, check_only=False)
        row = dict(case="test", n=16, accuracy=dict(max_absolute_error=0),
                   scifort_seconds=[.001, .002, .003], scipy_seconds=[.003, .004, .005],
                   scifort_median_seconds=.002, scipy_median_seconds=.004, scipy_over_scifort=2)
        text = speed.markdown_report(dict(base, cases=[row]))
        self.assertIn("| test | 16 | 2 | 4 | 2.000 | 0 |", text)
        self.assertIn("SciFort 1, 2, 3; SciPy 3, 4, 5.", text)
        self.assertIn("| All | 1 | 2.000 | 2.000 | 2.000 | n/a | 2.000 | 2.000 |", text)
        text = speed.markdown_report(dict(base, check_only=True,
                                         cases=[dict(case="test", n=16, accuracy=row["accuracy"])]))
        self.assertIn("| test | 16 | - | - | - | 0 |", text)
        self.assertNotIn("## Summary", text)


if __name__ == "__main__":
    unittest.main()
