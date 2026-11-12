#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent SciPy comparison for real/complex exponential condition numbers."""
from __future__ import annotations

import argparse
import os
import subprocess

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy import linalg as la


def make_cases(count: int, seed: int, max_order: int):
    rng = np.random.default_rng(seed)
    cases = []
    for k in range(count):
        complex_input = (k % 2 == 1)
        n = int(rng.integers(1, max_order+1))
        z = rng.standard_normal((n, n))
        if complex_input:
            z = z + 1j*rng.standard_normal((n, n))
        kind = (k // 2) % 13
        if kind == 0:
            a = np.zeros((n,n), dtype=np.complex128 if complex_input else float)
        elif kind == 1:
            alpha = rng.uniform(-5,5) + (1j*rng.uniform(-3,3) if complex_input else 0)
            a = alpha*np.eye(n)
        elif kind == 2:
            a = np.diag(rng.uniform(-4, 4, n)).astype(complex if complex_input else float)
            if complex_input:
                a += 1j*np.diag(rng.uniform(-4,4,n))
        elif kind == 3:
            a = z*0.01
        elif kind == 4:
            a = z*0.5
        elif kind == 5:
            a = np.triu(z*0.8, 1)  # nilpotent
        elif kind == 6:
            a = z*2.0
        elif kind == 7:
            a = (z + z.conj().T)*0.75 if complex_input else (z + z.T)*0.75
        elif kind == 8:
            a = (z-z.conj().T)*0.9 if complex_input else (z-z.T)*0.9
        elif kind == 9:
            a = np.triu(z*0.3)
        elif kind == 10:
            a = np.diag(rng.uniform(-1,1,n)) + np.triu(z*0.05, 1)
        elif kind == 11:
            a = (z-z.conj().T)*0.25 + 1.5*np.eye(n)
        else:
            a = z*1e-9
        cases.append((1+int(complex_input), np.asarray(a, dtype=np.complex128)))
    return cases


def format_cases(cases):
    lines = [str(len(cases))]
    for dtype, a in cases:
        lines.append(f"{dtype} {len(a)}")
        for value in a.flat:
            lines.append(f"{value.real:.17e} {value.imag:.17e}")
    return "\n".join(lines)+"\n"


def compare(executable, count, seed, max_order):
    cases = make_cases(count, seed, max_order)
    output = subprocess.run([executable], input=format_cases(cases), text=True,
                            check=True, capture_output=True)
    lines = output.stdout.splitlines()
    assert len(lines) == len(cases), (len(lines), len(cases), output.stderr)
    failures = []
    errors = []
    iterations = []
    for k, ((dtype,a),line) in enumerate(zip(cases,lines)):
        status, nit, condition = line.split()
        status, nit, condition = int(status), int(nit), float(condition)
        ref = float(la.expm_cond(a if dtype == 2 else a.real))
        if status != 0:
            failures.append((k, dtype, a.shape[0], status, nit, ref, a))
            continue
        err = abs(condition-ref)/max(abs(ref), 1e-14)
        errors.append(err)
        iterations.append(nit)
        if not np.isfinite(condition) or err > 5e-6:
            failures.append((k,dtype,a.shape[0],status,nit,condition,ref,a))
    print(f"cases={len(cases)} successes={len(errors)} failures={len(failures)} "
          f"worst_relative={max(errors,default=0.0):.3e} "
          f"avg_iterations={np.mean(iterations) if iterations else 0:.1f} "
          f"max_iterations={max(iterations,default=0)}")
    if failures:
        for entry in failures[:12]:
            print("FAIL", entry)
        raise AssertionError("expm_cond failed SciPy comparison")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("executable")
    p.add_argument("--cases", type=int, default=300)
    p.add_argument("--seed", type=int, default=20261107)
    p.add_argument("--max-order", type=int, default=6)
    args = p.parse_args()
    compare(args.executable, args.cases, args.seed, args.max_order)
