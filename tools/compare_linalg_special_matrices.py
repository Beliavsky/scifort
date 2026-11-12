#!/usr/bin/env python3
"""Compare SciFort Helmert/Hilbert constructors with SciPy and exact arithmetic.

Usage:
    python tools/compare_linalg_special_matrices.py path/to/compiled/driver

This checks matrix entries directly, including all exact integer elements.
It uses scipy.linalg as an independent oracle; the installed SciPy package
is not part of SciFort's runtime dependencies.
"""

from __future__ import annotations

import argparse
import dataclasses
import subprocess

import numpy as np
import scipy
from scipy.linalg import helmert, hilbert, invhilbert


@dataclasses.dataclass(frozen=True)
class Case:
    mode: str
    n: int
    bad: int = 0


def cases() -> list[Case]:
    result: list[Case] = []
    for n in range(0, 65):
        for mode in "MFHIED":
            # n > 14 is intentionally rejected by the exact int64 overload.
            result.append(Case(mode, n))
    for n in (72, 80, 96, 128, 250):
        result.append(Case("H", n))
        result.append(Case("I", n))
    for n in (-3, -1, 0, 1, 2, 5, 12, 15, 20):
        for mode in "MFHIED":
            result.append(Case(mode, n, 1))
    for mode in "MFHIED":
        result.append(Case(mode, -1))
    # Repeated independent checks of shape handling with substantial orders.
    for n in (25, 41, 64):
        for mode in "MFHIED":
            result.append(Case(mode, n, 1))
    return result


def expected(c: Case) -> tuple[int, np.ndarray]:
    n = max(0, c.n)
    nr = max(0, n - 1) if c.mode == "M" else n
    nr += c.bad
    target = np.zeros((nr, n), dtype=object if c.mode in "ED" else float)
    if c.bad or c.n < 0 or (c.mode in "MF" and c.n == 0):
        return 1, target
    if c.mode == "E" and n > 14:
        return 1, target
    if c.mode == "M":
        target = helmert(n, full=False)
    elif c.mode == "F":
        target = helmert(n, full=True)
    elif c.mode == "H":
        target = hilbert(n)
    elif c.mode == "I":
        if n == 250:
            return 4, target
        # Compare the floating constructor with SciPy's *exact* integer
        # reference. SciPy's default float path can itself lose ~1e-12
        # relative precision for orders near 128.
        target = np.asarray(invhilbert(n, exact=True), dtype=float)
    else:
        target = invhilbert(n, exact=True)
    return 0, np.asarray(target)


def run(executable: str) -> None:
    all_cases = cases()
    input_data = str(len(all_cases)) + "\n" + "\n".join(
        f"{c.mode} {c.n} {c.bad}" for c in all_cases
    ) + "\n"
    proc = subprocess.run(
        [executable], input=input_data, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
        timeout=180,
    )
    lines = iter(proc.stdout.splitlines())
    max_relative = 0.0
    exact_entries = 0
    floating_entries = 0
    errors = 0
    for ci, case in enumerate(all_cases):
        status, exp = expected(case)
        nr, nc = exp.shape
        try:
            actual_status = int(next(lines))
        except (StopIteration, ValueError) as err:
            raise AssertionError(f"could not read status for case {ci}: {case}") from err
        if actual_status != status:
            raise AssertionError(
                f"case {ci} {case}: status {actual_status}, expected {status}"
            )
        rows = []
        for _ in range(nr):
            try:
                values = next(lines).split()
            except StopIteration as err:
                raise AssertionError(f"missing output row for case {ci}: {case}") from err
            if len(values) != nc and not (case.mode == "D" and status and len(values) == 0):
                raise AssertionError(f"case {ci} {case}: expected {nc} columns; got {values}")
            if case.mode == "D" and status:
                values = [0] * nc
            rows.append(values)
        if case.mode in "ED":
            got = np.array(rows, dtype=object).reshape(nr, nc)
            for i in range(nr):
                for j in range(nc):
                    v = int(got[i, j])
                    w = int(exp[i, j])
                    if v != w:
                        raise AssertionError((ci, case, i, j, v, w))
                    exact_entries += 1
        else:
            got = np.array(rows, dtype=float).reshape(nr, nc)
            if not np.allclose(got, exp, atol=1e-15, rtol=6e-13):
                abs_err = np.abs(got - exp)
                ij = np.unravel_index(np.argmax(abs_err), exp.shape)
                raise AssertionError(
                    (ci, case, ij, got[ij], exp[ij], abs_err[ij])
                )
            floating_entries += nr * nc
            if nr and nc and status == 0:
                rel = np.max(np.abs(got - exp) / np.maximum(1.0, np.abs(exp)))
                max_relative = max(max_relative, float(rel))
        if status:
            errors += 1
    try:
        trailing = next(lines)
    except StopIteration:
        pass
    else:
        raise AssertionError(f"unexpected trailing output: {trailing[:100]}")
    print(
        f"PASS: {len(all_cases)} SciPy matrix comparisons; "
        f"{floating_entries} floating entries, {exact_entries} exact entries; "
        f"{errors} expected error statuses; "
        f"max scaled difference {max_relative:.3e} "
        f"(SciPy {scipy.__version__})"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("driver", help="compiled Fortran reference driver")
    options = parser.parse_args()
    run(options.driver)


if __name__ == "__main__":
    main()
