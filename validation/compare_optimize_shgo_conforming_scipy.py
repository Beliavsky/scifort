#!/usr/bin/env python3
"""Compare conforming 2D SHGO samples with SciPy optimized reference values.

The graph triangulation and sampling strategy differ from SciPy's SHGO.
These are quality checks rather than equality of iterates or convergence proofs.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy.optimize import shgo


def quadratic(x: np.ndarray) -> float:
    return float(np.sum((x - [0.625, -0.425]) ** 2))


def himmelblau(x: np.ndarray) -> float:
    return float((x[0] ** 2 + x[1] - 11) ** 2 + (x[0] + x[1] ** 2 - 7) ** 2)


def rastrigin(x: np.ndarray) -> float:
    return float(20 + np.sum(x**2 - 10 * np.cos(2 * np.pi * x)))


def ackley(x: np.ndarray) -> float:
    return float(-20 * np.exp(-0.2 * np.sqrt(np.sum(x**2) / 2))
                 - np.exp(np.sum(np.cos(2 * np.pi * x)) / 2) + 20 + np.e)


def six_hump(x: np.ndarray) -> float:
    return float((4 - 2.1*x[0]**2 + x[0]**4 / 3)*x[0]**2
                 + x[0]*x[1] + (-4 + 4*x[1]**2)*x[1]**2)


CASES = {
    "quadratic": (quadratic, [(-2, 2)] * 2, 0.03, 81),
    "himmelblau": (himmelblau, [(-6, 6)] * 2, 0.5, 1089),
    "rastrigin": (rastrigin, [(-5.12, 5.12)] * 2, 0.05, 289),
    "ackley": (ackley, [(-5, 5)] * 2, 0.05, 81),
    "six_hump": (six_hump, [(-3, 3), (-2, 2)], 0.06, 1089),
}


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_optimize_shgo_conforming_scipy.py DRIVER")
    lines = subprocess.check_output([sys.argv[1]], text=True).splitlines()
    assert len(lines) == len(CASES), (len(lines), len(CASES))
    for line in lines:
        name, status_s, nfev_s, nleaf_s, ntri_s, valid_s, f_s, *coords = line.split()
        f, bounds, tol, full_grid = CASES[name]
        status, nfev, nleaf, ntri, valid = map(
            int, [status_s, nfev_s, nleaf_s, ntri_s, valid_s]
        )
        x = np.array([float(v) for v in coords])
        value = float(f_s)
        assert status == 0 and valid == 1, (name, status, valid)
        assert nfev < full_grid and nleaf >= 4 and ntri >= 4*nleaf
        assert len(x) == 2 and np.all(np.isfinite(x))
        lower, upper = np.array(bounds).T
        assert np.all(x >= lower) and np.all(x <= upper), (name, x)
        assert abs(value - f(x)) < 1e-9*max(1.0, abs(value)), (name, x, value)
        reference = shgo(f, bounds, n=64, iters=2, sampling_method="sobol")
        assert reference.success and np.isfinite(reference.fun), (name, reference.message)
        assert value <= float(reference.fun) + tol, (name, value, reference.fun)
        print(f"{name}: conforming={value:.6g}; SciPy={reference.fun:.6g}; "
              f"vertices={nfev} triangles={ntri} leaves={nleaf}")
    print(f"{len(CASES)} conforming SHGO/SciPy quality comparisons passed")


if __name__ == "__main__":
    main()
