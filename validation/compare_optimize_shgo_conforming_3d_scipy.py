#!/usr/bin/env python3
"""Check deterministic 3D SciFort samples against SciPy SHGO references.

SciFort's conforming octree is not the same sampling or local-polishing
algorithm as SciPy SHGO; comparisons are on attainable objective quality.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy.optimize import shgo


def quadratic(x: np.ndarray) -> float:
    return float(np.sum((x - [0.3, -0.2, 0.15]) ** 2))


def rastrigin(x: np.ndarray) -> float:
    return float(30.0 + np.sum(x * x - 10.0 * np.cos(2.0 * np.pi * x)))


def ackley(x: np.ndarray) -> float:
    return float(-20.0 * np.exp(-0.2 * np.sqrt(np.sum(x * x) / 3.0))
                 - np.exp(np.sum(np.cos(2.0 * np.pi * x)) / 3.0)
                 + 20.0 + np.e)


def double_well(x: np.ndarray) -> float:
    return float((x[0] ** 2 - 1.0) ** 2 + (x[1] - 0.35) ** 2
                 + (x[2] + 0.55) ** 2)


def rosenbrock(x: np.ndarray) -> float:
    return float(100.0 * (x[1] - x[0] ** 2) ** 2 + (1.0 - x[0]) ** 2
                 + 100.0 * (x[2] - x[1] ** 2) ** 2 + (1.0 - x[1]) ** 2)


# Tolerances allow for intentionally unpolished grid samples.
CASES = {
    "quadratic": (quadratic, [(-1, 1)] * 3, 0.01),
    "rastrigin": (rastrigin, [(-5.12, 5.12)] * 3, 0.01),
    "ackley": (ackley, [(-5, 5)] * 3, 0.01),
    "double_well": (double_well, [(-2, 2)] * 3, 0.06),
    "rosenbrock": (rosenbrock, [(-3, 3)] * 3, 0.8),
}


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_optimize_shgo_conforming_3d_scipy.py DRIVER")
    lines = subprocess.check_output([sys.argv[1]], text=True).splitlines()
    assert len(lines) == len(CASES), (len(lines), len(CASES))
    for line in lines:
        name, status_s, nfev_s, nleaf_s, ntet_s, valid_s, f_s, *coords = line.split()
        objective, bounds, tolerance = CASES[name]
        status, nfev, nleaf, ntet, valid = map(
            int, [status_s, nfev_s, nleaf_s, ntet_s, valid_s]
        )
        x = np.array([float(v) for v in coords])
        value = float(f_s)
        assert status == 0 and valid == 1, (name, status, valid)
        assert nfev >= 71 and nleaf >= 8 and ntet >= 24 * 8, name
        assert len(x) == 3 and np.all(np.isfinite(x)), name
        low, high = np.array(bounds).T
        assert np.all(x >= low) and np.all(x <= high), (name, x)
        assert abs(value - objective(x)) <= 1e-9 * max(1.0, abs(value)), name
        reference = shgo(objective, bounds, n=64, iters=2, sampling_method="sobol")
        assert reference.success and np.isfinite(reference.fun), (name, reference.message)
        assert value <= float(reference.fun) + tolerance, (name, value, reference.fun)
        print(f"{name}: conforming-3D={value:.6g}; SciPy={reference.fun:.6g}; "
              f"vertices={nfev} tets={ntet} leaves={nleaf}")
    print(f"{len(CASES)} conforming 3D SHGO/SciPy quality comparisons passed")


if __name__ == "__main__":
    main()
