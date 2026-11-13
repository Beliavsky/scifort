#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Validate the larger AAA arrowhead QZ backend against LAPACK and analytic roots."""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.linalg import eigvals
from scipy.optimize import linear_sum_assignment


def roots(alpha, beta):
    keep = abs(beta) > 64*np.finfo(float).eps*abs(alpha)
    return alpha[keep]/beta[keep]


def matching(a, b):
    assert len(a) == len(b), (len(a), len(b))
    if not len(a):
        return 0.0
    cost = abs(a[:, None]-b[None, :])/(1+abs(b)[None, :])
    i, j = linear_sum_assignment(cost)
    return np.max(cost[i, j])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261117)
    cases = []
    for m in (2, 16, 40, 41, 64, 100, 128, 200, 512):
        s = np.exp(2j*np.pi*np.arange(m)/m)
        u = (1-0.4*s)/1.4
        exact = 0.4**(1/(m-1))*np.exp(2j*np.pi*np.arange(m-1)/(m-1))
        cases.append((s, u, exact))
        if m <= 200:
            for kind in ("real", "complex", "clustered", "degree_loss"):
                if kind == "real":
                    s = np.linspace(-1, 1, m).astype(complex)
                    u = rng.normal(size=m).astype(complex)
                else:
                    s = rng.uniform(-0.7, 0.7, m)+1j*rng.uniform(-0.7, 0.7, m)
                    u = rng.normal(size=m)+1j*rng.normal(size=m)
                    if kind == "clustered":
                        s[:m//2] = 0.2+0.3j+0.001*s[:m//2]
                    if kind == "degree_loss":
                        u[-1] = -sum(u[:-1])
                u = u/max(abs(u))
                cases.append((s, u, None))
    lines = []
    for s, u, _ in cases:
        lines.append(str(len(s)))
        lines += [f"{z.real:.17g} {z.imag:.17g} {w.real:.17g} {w.imag:.17g}"
                  for z, w in zip(s, u)]
    result = subprocess.run([str(args.executable.resolve())], input="\n".join(lines)+"\n",
                            capture_output=True, text=True, check=True)
    rows = iter(result.stdout.splitlines())
    maximum = analytic_maximum = 0.0
    count = analytic_count = 0
    for s, u, exact in cases:
        status = int(next(rows))
        assert status == 0, (len(s), status)
        ab = np.array([np.fromstring(next(rows), sep=" ") for _ in range(len(s)+1)])
        actual = roots(ab[:, 0]+1j*ab[:, 1], ab[:, 2]+1j*ab[:, 3])
        n = len(s)+1
        a = np.zeros((n, n), dtype=complex)
        b = a.copy()
        a[0, 1:] = u
        a[1:, 0] = 1
        a[np.arange(1, n), np.arange(1, n)] = s
        b[np.arange(1, n), np.arange(1, n)] = 1
        ref = eigvals(a, b, homogeneous_eigvals=True)
        expected = roots(*ref)
        error = matching(actual, expected)
        assert error < 2e-8, (len(s), error)
        maximum = max(maximum, error)
        count += len(actual)
        if exact is not None:
            error = matching(actual, exact)
            assert error < 3e-10, (len(s), error)
            analytic_maximum = max(analytic_maximum, error)
            analytic_count += len(actual)
    assert next(rows, None) is None
    print(f"PASS: {len(cases)} pencils, {count} finite SciPy eigenvalues, "
          f"{analytic_count} analytic roots; SciPy {scipy.__version__}")
    print(f"Maximum normalized matching errors: SciPy {maximum:.6g}, "
          f"analytic {analytic_maximum:.6g}")


if __name__ == "__main__":
    main()
