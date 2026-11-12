#!/usr/bin/env python3
"""Reproducible seeded Schur-versus-SciPy shifted-system comparisons."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import scipy.linalg as la


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(14062026)
    cases = []
    lines = []
    for case in range(150):
        ns = int(rng.integers(1, 29))
        ni = int(rng.integers(1, 5))
        no = int(rng.integers(1, 5))
        nw = int(rng.integers(1, 9))
        discrete = case % 2
        # Include complex, real, triangular, and strongly nonnormal cases.
        complex_case = case % 3 != 0
        if case % 7 == 0:
            a = np.diag(-0.3 - 0.03 * np.arange(ns)).astype(complex)
            if ns > 1:
                a += np.diag(np.full(ns - 1, 0.8), k=1)
        else:
            a = 0.06 * (rng.normal(size=(ns, ns)) + 1j * rng.normal(size=(ns, ns)))
            a += np.diag(-0.4 - 0.03 * np.arange(ns))
        if not complex_case:
            a = a.real.astype(complex)
        b = rng.normal(size=(ns, ni)) + (1j * rng.normal(size=(ns, ni)) if complex_case else 0)
        c = rng.normal(size=(no, ns)) + (1j * rng.normal(size=(no, ns)) if complex_case else 0)
        d = 0.03 * (rng.normal(size=(no, ni)) + (1j * rng.normal(size=(no, ni)) if complex_case else 0))
        w = rng.uniform(0.05, 2.8, size=nw)
        cases.append((a, b, c, d, w, discrete))
    lines.append(str(len(cases)))
    for a, b, c, d, w, discrete in cases:
        ns, ni = b.shape
        no = c.shape[0]
        lines.append(f"{ns} {ni} {no} {len(w)} {discrete}")
        for mat in (a, b, c, d):
            for z in mat.flatten(order="F"):
                lines.append(f"{z.real:.17e} {z.imag:.17e}")
        lines.append(" ".join(f"{v:.17e}" for v in w))
    process = subprocess.run([str(args.executable)], input="\n".join(lines) + "\n", text=True,
                             capture_output=True, check=True)
    parsed = iter(process.stdout.splitlines())
    max_error = 0.0
    schur_count = 0
    total = 0
    for a, b, c, d, w, discrete in cases:
        method, status = next(parsed).split()
        assert int(status) == 0, (method, status)
        schur_count += method == "schur"
        no, ni = d.shape
        shape = (len(w), ni, no)
        vals = np.empty(shape, dtype=complex)
        for k in range(len(w)):
            for j in range(ni):
                for i in range(no):
                    real, imag = map(float, next(parsed).split())
                    vals[k, j, i] = real + 1j * imag
        for k, omega in enumerate(w):
            z = np.exp(1j * omega) if discrete else 1j * omega
            reference = c @ la.solve(z * np.eye(len(a)) - a, b) + d
            scaled = np.abs(vals[k].T-reference) / np.maximum(1.0, np.abs(reference))
            max_error = max(max_error, float(np.max(scaled)))
            assert float(np.max(scaled)) < 2e-9, (method, scaled.max())
            total += reference.size
    print(f"PASS: {len(cases)} seeded MIMO cases, {total} channel values, "
          f"{schur_count} Schur-converged plans, maximum scaled error {max_error:.3e}")


if __name__ == "__main__":
    main()
