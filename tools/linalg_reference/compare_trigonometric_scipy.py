#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent SciPy comparison for six real/complex dense matrix functions."""
from __future__ import annotations

import os
import subprocess

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy import linalg as la

REFERENCE = (la.sinm, la.cosm, la.tanm, la.sinhm, la.coshm, la.tanhm)


def make_cases(count: int, seed: int, max_order: int = 6):
    rng = np.random.default_rng(seed)
    cases = []
    for k in range(count):
        op = 1 + (k % 6)
        is_complex = (k // 6) % 2 == 1
        n = int(rng.integers(1, max_order + 1))
        group = (k // 12) % 10
        z = rng.normal(size=(n, n))
        if is_complex:
            z = z + 1j * rng.normal(size=(n, n))
        if group == 0:
            a = np.diag(rng.uniform(-2, 2, size=n))
            if is_complex:
                a = a + 1j * np.diag(rng.uniform(-1.5, 1.5, size=n))
        elif group == 1:
            a = z * 0.15
        elif group == 2:
            a = np.triu(z * 0.5)
            np.fill_diagonal(a, rng.uniform(-0.5, 0.5, n))
        elif group == 3:
            a = np.triu(z * 1.4, 1)  # nilpotent
        elif group == 4:
            a = z * 1.2
        elif group == 5:
            a = z * 3.0
        elif group == 6:
            a = z * (1e-12 if k % 2 == 0 else 1e-70)
        elif group == 7:
            a = (z - z.conj().T) * 0.7 if is_complex else (z - z.T) * 0.7
        elif group == 8:
            a = (z + z.conj().T) * 0.9 if is_complex else (z + z.T) * 0.9
        else:
            a = np.diag(rng.uniform(-0.25, 0.25, size=n)).astype(complex if is_complex else float)
            if op == 3:
                a[0, 0] = np.pi/2 + (1e-3 if k % 2 else -1e-3)
            elif op == 6 and is_complex:
                a[0, 0] = 1j * (np.pi/2 + (1e-3 if k % 2 else -1e-3))
            a += np.triu(z * 0.03, 1)
        cases.append((op, int(is_complex) + 1, np.asarray(a, dtype=np.complex128)))
    return cases


def format_cases(cases):
    lines = [str(len(cases))]
    for op, dtype, a in cases:
        lines.append(f"{op} {dtype} {len(a)}")
        for value in a.flat:
            lines.append(f"{value.real:.17e} {value.imag:.17e}")
    return "\n".join(lines) + "\n"


def compare(executable: str, count: int, seed: int, max_order: int = 6):
    cases = make_cases(count, seed, max_order)
    process = subprocess.run([executable], input=format_cases(cases), check=True,
                             text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    data = iter(process.stdout.splitlines())
    errors = []
    rejected = []
    expected_ill_conditioned = []
    mismatches = []
    per_function = [0] * 6
    for index, (op, dtype, a) in enumerate(cases):
        status = int(next(data))
        n = len(a)
        x = np.array([complex(*map(float, next(data).split()))
                      for _ in range(n * n)]).reshape(n, n)
        reference_input = a if dtype == 2 else a.real
        if status != 0:
            # Tangent is defined only if cos(A), or cosh(A), is invertible;
            # refusing an extremely ill-conditioned solve is intentional.
            if op in (3, 6) and status in (4, 5):
                denom = la.cosm(reference_input) if op == 3 else la.coshm(reference_input)
                condition = np.linalg.cond(denom)
                if not np.isfinite(condition) or condition > 1e10:
                    expected_ill_conditioned.append((index, op, dtype, n, status, condition))
                    continue
            rejected.append((index, op, dtype, n, status))
            continue
        reference = np.asarray(REFERENCE[op-1](reference_input), dtype=np.complex128)
        diff = np.linalg.norm(x-reference, 1) / max(1.0, np.linalg.norm(reference, 1))
        errors.append(diff)
        per_function[op-1] += 1
        if not np.isfinite(diff) or diff > 2e-8:
            mismatches.append((index, op, dtype, n, float(diff), a, x, reference))
        if not np.all(np.isfinite(x)):
            mismatches.append((index, op, dtype, n, np.inf, a, x, reference))
    print(f"count={count} successes={len(errors)} failures={len(rejected)} "
          f"ill_conditioned={len(expected_ill_conditioned)} "
          f"mismatches={len(mismatches)} worst={max(errors, default=0):.3e} "
          f"per_function={per_function}")
    if rejected:
        print("unexpected failures:", rejected[:15])
    if expected_ill_conditioned:
        print("expected ill-conditioned tangent rejections:", expected_ill_conditioned[:8])
    if mismatches:
        for entry in sorted(mismatches, key=lambda v: -v[4])[:8]:
            print("mismatch: ", entry[:5], "input:", entry[5],
                  "got:", entry[6], "reference:", entry[7], sep="\n")
    if rejected or mismatches:
        raise AssertionError("SciFort trigonometric reference regression failed")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("executable")
    p.add_argument("--cases", type=int, default=600)
    p.add_argument("--seed", type=int, default=20261107)
    p.add_argument("--max-order", type=int, default=6)
    args = p.parse_args()
    compare(args.executable, args.cases, args.seed, args.max_order)
