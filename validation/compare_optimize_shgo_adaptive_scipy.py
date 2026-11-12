#!/usr/bin/env python3
"""Validate selective SHGO samples against the known objective and SciPy SHGO.

The adaptive graph and SciPy's algorithm differ; this checks numerical
objective agreement and reasonable proximity to independently refined minima,
not identical evaluation order, outputs, or global-optimality guarantees.
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
    d = x.size
    return float(-20 * np.exp(-0.2 * np.sqrt(np.sum(x**2) / d))
                 - np.exp(np.sum(np.cos(2 * np.pi * x)) / d) + 20 + np.e)


def double_well(x: np.ndarray) -> float:
    return float((x[0] ** 2 - 1) ** 2 + 0.1 * x[0])


def quadratic3(x: np.ndarray) -> float:
    return float(np.sum((x - [0.625, -0.425, 0.775]) ** 2))


# Function, bounds, absolute difference to SciPy-refined objective, maximum
# sample budget expected for the specified adaptive test configuration.
CASES = {
    "quadratic": (quadratic, [(-2, 2)] * 2, 0.03, 81),
    "himmelblau": (himmelblau, [(-6, 6)] * 2, 2.5, 1089),
    "rastrigin": (rastrigin, [(-5.12, 5.12)] * 2, 0.001, 289),
    "ackley": (ackley, [(-5, 5)] * 2, 0.001, 81),
    "double_well": (double_well, [(-2, 2)], 0.003, 9),
    "quadratic3": (quadratic3, [(-2, 2)] * 3, 0.04, 4913),
}


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_optimize_shgo_adaptive_scipy.py DRIVER")
    lines = subprocess.check_output([sys.argv[1]], text=True).splitlines()
    assert len(lines) == len(CASES), (len(lines), len(CASES))
    nchecked = 0
    for line in lines:
        name, raw_nfev, raw_nmesh, raw_nleaf, raw_status, raw_f, *raw_x = line.split()
        fun, bounds, tolerance, complete_grid = CASES[name]
        x = np.array([float(v) for v in raw_x])
        value = float(raw_f)
        nfev, nmesh, nleaf, status = map(
            int, [raw_nfev, raw_nmesh, raw_nleaf, raw_status]
        )
        assert status == 0, (name, status)
        assert nfev == nmesh and nmesh >= 3 and nleaf >= 2, (name, nfev, nmesh, nleaf)
        assert nfev < complete_grid, (name, nfev, complete_grid)
        assert len(x) == len(bounds), (name, x)
        assert np.all(x >= np.array(bounds)[:, 0]) and np.all(x <= np.array(bounds)[:, 1])
        assert abs(value - fun(x)) < 1e-9 * max(1.0, abs(value)), (name, value, x)
        reference = shgo(fun, bounds, n=64, iters=2, sampling_method="sobol")
        assert reference.success and np.isfinite(reference.fun), (name, reference.message)
        assert value <= float(reference.fun) + tolerance, (
            name, value, float(reference.fun), tolerance
        )
        nchecked += 1
        print(f"{name}: adaptive={value:.5g}; SciPy={reference.fun:.5g}; "
              f"adaptive_nfev={nfev} vs complete_mesh={complete_grid}")
    print(f"{nchecked} adaptive SHGO/SciPy quality comparisons passed")


if __name__ == "__main__":
    main()
