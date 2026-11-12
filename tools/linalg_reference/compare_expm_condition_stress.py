#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent diagonal divided-difference and scalar-shift stress tests."""
from __future__ import annotations

import os
import subprocess

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy import linalg as la
from compare_expm_condition_scipy import format_cases


def exact_diagonal_condition(values):
    """Exact induced Frobenius norm from the exp divided differences."""
    values = np.asarray(values, dtype=np.complex128)
    shifted = values - max(values.real)
    exponentials = np.exp(shifted)
    differences = values[:, None]-values[None, :]
    exponent_differences = exponentials[:, None]-exponentials[None, :]
    dd = np.empty_like(differences)
    close = np.abs(differences) < 1e-2
    far = ~close
    dd[far] = exponent_differences[far]/differences[far]
    # np.expm1 on complex arguments avoids subtractive cancellation.
    close_diff = differences[close]
    ratio = np.empty_like(close_diff)
    nonzero = np.abs(close_diff) > 1e-14
    ratio[nonzero] = np.expm1(close_diff[nonzero])/close_diff[nonzero]
    ratio[~nonzero] = 1.0 + close_diff[~nonzero]/2
    ei = np.broadcast_to(exponentials[None, :], differences.shape)
    dd[close] = ei[close]*ratio
    return float(np.linalg.norm(values) * np.max(np.abs(dd)) / np.linalg.norm(exponentials))


def make_stress_cases(seed=20261107):
    rng = np.random.default_rng(seed)
    cases = []
    refs = []
    for n in (1, 2, 3, 6, 12, 24, 40, 48):
        for mode in range(8):
            if mode == 0:
                eigs = np.full(n, -300.0)
            elif mode == 1:
                eigs = np.linspace(-1000.0, 1000.0, n)
            elif mode == 2:
                eigs = 1e-10*rng.normal(size=n)
            elif mode == 3:
                eigs = 600.0 + 1j*np.linspace(-30.0, 30.0, n)
            elif mode == 4:
                eigs = rng.uniform(-4.0, 4.0, size=n)
            elif mode == 5:
                eigs = rng.normal(size=n)+1j*rng.normal(size=n)
            elif mode == 6:
                eigs = -700.0 + rng.uniform(-2.0, 2.0, size=n)
            else:
                eigs = -300.0 + 50j + rng.normal(size=n)*0.5 + 1j*rng.normal(size=n)*0.4
            dtype = 2 if np.iscomplexobj(eigs) else 1
            a = np.diag(eigs).astype(np.complex128)
            cases.append((dtype, a))
            refs.append(exact_diagonal_condition(eigs))
    # Dense tests: scalar-shift invariance provides a separate oracle without
    # asking SciPy to exponentiate extreme positive or negative arguments.
    for n in (2, 4, 7, 12):
        for shift in (-600.0, -300.0, 500.0, 650.0):
            is_complex = n in (4, 12)
            a = rng.normal(size=(n,n))*0.12
            if is_complex:
                a = a+1j*rng.normal(size=(n,n))*0.12
            # Include nonnormal upper triangular and general matrices.
            if n < 7:
                a = np.triu(a)
            reference = la.expm_cond(a)*np.linalg.norm(a+shift*np.eye(n))/np.linalg.norm(a)
            cases.append((2 if is_complex else 1,
                          np.asarray(a+shift*np.eye(n), dtype=np.complex128)))
            refs.append(reference)
    return cases, refs


def compare(executable, seed=20261107):
    cases, references = make_stress_cases(seed)
    process = subprocess.run([executable], input=format_cases(cases), text=True,
                             check=True, capture_output=True)
    lines = process.stdout.splitlines()
    assert len(lines) == len(cases), (len(lines), len(cases), process.stderr)
    errors, failures = [], []
    for k, ((dtype, a), reference, line) in enumerate(zip(cases, references, lines)):
        status, steps, value = line.split()
        status, value = int(status), float(value)
        error = abs(value-reference)/max(1e-14, abs(reference))
        errors.append(error)
        if status != 0 or not np.isfinite(value) or error > 1e-6:
            failures.append((k, len(a), dtype, status, steps, value, reference, error))
    print(f"stress_cases={len(cases)} failures={len(failures)} "
          f"worst_relative={max(errors,default=0):.3e}")
    for failed in failures[:12]:
        print("FAIL", failed)
    if failures:
        raise AssertionError("diagonal/shifted expm_cond stress tests failed")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("executable")
    p.add_argument("--seed", type=int, default=20261107)
    args = p.parse_args()
    compare(args.executable, args.seed)
