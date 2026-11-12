#!/usr/bin/env python3
"""Compare SciFort companion/Fiedler constructors directly with SciPy 1.17+.

Usage:
  python tools/compare_linalg_companion.py /path/to/compiled/linalg_companion_driver

The executable reads a text stream of test cases. This script generates all
cases locally, checks entries against scipy.linalg and verifies eigenvalue
agreement for nontrivial random polynomial companion matrices.
"""

from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass

import numpy as np
from scipy.linalg import companion, fiedler, fiedler_companion
from scipy.optimize import linear_sum_assignment


@dataclass
class Case:
    mode: str
    dtype: str
    a: np.ndarray


def build_cases() -> list[Case]:
    rng = np.random.default_rng(20261109)
    cases: list[Case] = []
    for mode in "CPF":
        for dtype in "RZI":
            for n in range(0, 16):
                for trial in range(10):
                    if dtype == "I":
                        a = rng.integers(-12, 13, n).astype(np.int32)
                    elif dtype == "R":
                        a = rng.normal(size=n)
                    else:
                        a = rng.normal(size=n) + 1j * rng.normal(size=n)
                    if mode != "F" and n >= 2:
                        if trial == 0:
                            a[0] = 0  # error-path coverage
                        elif trial == 1:
                            a[0] = 1
                        elif trial == 2:
                            a[0] = -1
                        elif dtype == "I" and a[0] == 0:
                            a[0] = 2
                        if dtype != "I" and trial in (3, 4):
                            a *= 1e-40 if trial == 3 else 1e40
                    cases.append(Case(mode, dtype, a))
    # Longer polynomials, where the nonzero structure and parity alternate.
    for mode in "CPF":
        for dtype in "RZI":
            for n in (17, 18, 23, 24, 31, 32):
                for trial in range(5):
                    if dtype == "I":
                        a = rng.integers(-8, 9, n).astype(np.int32)
                    elif dtype == "R":
                        a = rng.normal(size=n)
                    else:
                        a = rng.normal(size=n) + 1j * rng.normal(size=n)
                    if mode != "F" and a[0] == 0:
                        a[0] = 1
                    cases.append(Case(mode, dtype, a))
    return cases


def format_row(x: np.ndarray) -> str:
    return " ".join(str(int(t)) if x.dtype.kind in "iu" else f"{float(t):.17e}" for t in x)


def build_input(cases: list[Case]) -> str:
    rows = [str(len(cases))]
    for case in cases:
        rows.append(f"{case.mode} {case.dtype} {len(case.a)}")
        if len(case.a):
            if case.dtype == "Z":
                rows.append(format_row(case.a.real))
                rows.append(format_row(case.a.imag))
            else:
                rows.append(format_row(case.a))
    return "\n".join(rows) + "\n"


def expected(case: Case) -> tuple[int, np.ndarray]:
    n = len(case.a)
    d = n if case.mode == "F" else max(0, n - 1)
    if case.mode == "C" and n < 2:
        return 1, np.zeros((d, d))
    if case.mode in "CP" and n >= 2 and case.a[0] == 0:
        return 1, np.zeros((d, d))
    if case.mode == "C":
        out = companion(case.a)
    elif case.mode == "P":
        out = fiedler_companion(case.a)
    else:
        out = fiedler(case.a)
    return 0, np.asarray(out).reshape((d, d))


def compare(exe: str) -> None:
    cases = build_cases()
    result = subprocess.run(
        [exe], input=build_input(cases), text=True, capture_output=True, check=True
    )
    lines = iter(result.stdout.splitlines())
    max_difference = 0.0
    success_count = 0
    expected_error_count = 0
    spectral_count = 0
    rng = np.random.default_rng(4327)
    for index, case in enumerate(cases):
        status, target = expected(case)
        d = target.shape[0]
        got_status = int(next(lines))
        if got_status != status:
            raise AssertionError((index, case.mode, case.dtype, len(case.a), status, got_status))
        got = np.empty((d, d), dtype=complex if case.dtype == "Z" and case.mode != "F" else float)
        for i in range(d):
            re = np.fromstring(next(lines), sep=" ")
            if len(re) != d:
                raise AssertionError((index, i, d, re))
            if case.dtype == "Z" and case.mode != "F":
                im = np.fromstring(next(lines), sep=" ")
                if len(im) != d:
                    raise AssertionError((index, i, d, im))
                got[i, :] = re + 1j * im
            else:
                got[i, :] = re
        error = np.linalg.norm(got - target)
        scale = max(1.0, float(np.linalg.norm(target)))
        difference = error / scale
        max_difference = max(max_difference, difference)
        if not np.allclose(got, target, rtol=3e-13, atol=5e-13, equal_nan=True):
            raise AssertionError((index, case.mode, case.dtype, case.a, got, target))
        if status:
            expected_error_count += 1
        else:
            success_count += 1
            # Independently check the mathematical invariant: both companion
            # matrices have the same characteristic roots as the polynomial.
            # Restrict this to moderate well-scaled data so that the numerical
            # eigenvalue comparison tests the construction, not conditioning.
            if (
                case.mode in "CP"
                and d >= 2
                and d <= 12
                and case.dtype != "I"
                and rng.random() < 0.25
                and 1e-8 < abs(case.a[0]) < 1e8
                and np.max(abs(case.a)) < 1e8
            ):
                lhs = np.linalg.eigvals(got)
                rhs = np.roots(case.a)
                cost = abs(lhs[:, None] - rhs[None, :])
                row, col = linear_sum_assignment(cost)
                if np.max(cost[row, col]) > 2e-5 * max(1.0, np.max(abs(rhs))):
                    raise AssertionError((index, "eigenvalues", lhs, rhs))
                spectral_count += 1
    try:
        extra = next(lines)
        raise AssertionError(f"Unexpected trailing driver output: {extra}")
    except StopIteration:
        pass
    print(
        f"PASS {len(cases)} SciPy comparisons "
        f"({success_count} valid, {expected_error_count} expected invalid); "
        f"{spectral_count} polynomial-root checks; "
        f"maximum relative Frobenius difference {max_difference:.3e}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("driver", help="Path to compiled Fortran reference driver")
    compare(parser.parse_args().driver)
