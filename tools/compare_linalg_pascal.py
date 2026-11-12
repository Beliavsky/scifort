#!/usr/bin/env python3
"""Compare SciFort Pascal constructors with independent scipy.linalg output.

Compile reference/linalg_pascal_driver.f90 with kinds, linalg, linalg_pascal.
Then run: python tools/compare_linalg_pascal.py PATH_TO_DRIVER

Tests every matrix entry for both Pascal and inverse Pascal, all three
matrix shapes, and floating-point/int64/arbitrary-precision decimal output.
"""

from __future__ import annotations

import argparse
import math
import random
import subprocess
from dataclasses import dataclass

from scipy.linalg import invpascal, pascal


@dataclass(frozen=True)
class Case:
    mode: str
    dtype: str
    kind: str
    n: int


def make_cases() -> list[Case]:
    rng = random.Random(20261109)
    cases = []
    for kind in ("symmetric", "lower", "upper"):
        for mode in ("P", "I"):
            for dtype, sizes in (
                ("R", [0, 1, 2, 3, 4, 5, 7, 9, 13, 20, 28, 35, 45, 64, 80, 120]),
                ("E", [0, 1, 2, 3, 4, 5, 7, 9, 13, 20, 28, 34, 35, 36, 50, 67, 68]),
                ("D", [0, 1, 2, 3, 4, 5, 7, 9, 13, 20, 28, 35, 40, 50, 65, 80]),
            ):
                cases.extend(Case(mode, dtype, kind, n) for n in sizes)
    # Extra randomized small/medium cases catch binomial and parity bugs.
    for _ in range(180):
        kind = rng.choice(("symmetric", "lower", "upper"))
        mode = rng.choice(("P", "I"))
        dtype = rng.choice(("R", "E", "D"))
        n = rng.randint(0, 30)
        cases.append(Case(mode, dtype, kind, n))
    return cases


def expected_overflow(case: Case) -> bool:
    if case.dtype != "E":
        return False
    limit = (35 if case.mode == "P" else 36) if case.kind == "symmetric" else 68
    return case.n >= limit


def evaluate(executable: str) -> None:
    cases = make_cases()
    payload = str(len(cases)) + "\n" + "\n".join(
        f"{c.mode} {c.dtype} {c.kind} {c.n}" for c in cases
    ) + "\n"
    process = subprocess.run(
        [executable], input=payload, text=True, capture_output=True, check=True
    )
    lines = iter(process.stdout.splitlines())
    maximum_scaled_error = 0.0
    floating_entries = 0
    exact_entries = 0
    expected_overflows = 0

    for number, case in enumerate(cases, start=1):
        status = int(next(lines))
        if expected_overflow(case):
            assert status == 4, (number, case, "expected checked overflow", status)
            expected_overflows += 1
            continue
        assert status == 0, (number, case, "unexpected failure", status)
        fun = pascal if case.mode == "P" else invpascal
        # exact=True is the reliable integer mathematical reference even for
        # the floating overload, avoiding roundoff in SciPy's exact=False mode.
        reference = fun(case.n, kind=case.kind, exact=True)
        for i in range(case.n):
            row = next(lines).split()
            assert len(row) == case.n, (number, case, i, "row length", len(row))
            for j, field in enumerate(row):
                target = int(reference[i, j])
                if case.dtype == "R":
                    actual = float(field)
                    assert math.isfinite(actual), (number, case, i, j, actual)
                    scaled_error = abs(actual - float(target)) / max(1.0, abs(float(target)))
                    maximum_scaled_error = max(maximum_scaled_error, scaled_error)
                    assert scaled_error < 5.0e-12, (number, case, i, j, actual, target)
                    floating_entries += 1
                else:
                    actual = int(field)
                    assert actual == target, (number, case, i, j, actual, target)
                    exact_entries += 1

    assert list(lines) == [], "Unexpected trailing output from Fortran driver"
    print(f"SciPy Pascal comparisons: {len(cases)}/{len(cases)} PASS")
    print(f"Floating entries compared: {floating_entries}")
    print(f"Exact integer / decimal entries compared: {exact_entries}")
    print(f"Expected int64 overflow statuses: {expected_overflows}")
    print(f"Maximum scaled floating difference: {maximum_scaled_error:.3e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", help="Path to compiled linalg_pascal_driver")
    evaluate(parser.parse_args().executable)
