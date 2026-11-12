#!/usr/bin/env python3
"""Independent benchmarks for SciFort's stochastic basin-hopping core.

Build validation/optimize_basinhopping_driver.f90 and supply its executable.
SciFort and SciPy use different random generators and local optimizers;
compare attained objective values, NOT stepwise trajectories or counts.
"""
from __future__ import annotations

import argparse
import math
import subprocess

import numpy as np
from scipy.optimize import basinhopping, minimize, minimize_scalar


def benchmark(case: int, x: np.ndarray) -> float:
    if case == 1:
        return float((x[0] ** 2 - 1.0) ** 2 + 0.25 * x[0])
    if case == 2:
        return float(20.0 + sum(x * x - 10.0 * np.cos(2 * np.pi * x)))
    if case == 3:
        return float((x[0] ** 2 + x[1] - 11.0) ** 2 + (x[0] + x[1] ** 2 - 7.0) ** 2)
    if case == 4:
        return float(-20.0 * np.exp(-0.2 * np.sqrt(np.mean(x * x)))
                     - np.exp(np.mean(np.cos(2.0 * np.pi * x))) + 20.0 + math.e)
    if case == 5:
        return float(np.sin(4.0 * x[0]) ** 2 + 0.03 * x[0] ** 2)
    if case == 6:
        return float((x[0] ** 2 - 1.0) ** 2 + 0.20 * x[0]
                     + 0.5 * (x[1] ** 2 - 1.0) ** 2 + 0.06 * x[1])
    raise ValueError(case)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("executable", help="compiled Fortran reference driver")
    args = p.parse_args()
    output = subprocess.check_output([args.executable], text=True)
    cases = {
        1: ([2.4], 1.5, 0.5),
        2: ([1.8, 1.7], 2.8, 2.0),
        3: ([2.5, 2.1], 3.5, 2.0),
        4: ([2.0, 1.8], 2.3, 1.5),
        5: ([2.0], 1.3, 0.5),
        6: ([2.0, 1.8], 1.8, 0.5),
    }
    root1 = minimize_scalar(lambda z: benchmark(1, np.array([z])),
                            bounds=(-2.0, -0.5), method="bounded",
                            options={"xatol": 1e-13}).fun
    root6 = minimize(lambda z: benchmark(6, z), [-1., -1.],
                     method="BFGS", tol=1e-12).fun
    target = {1: root1, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: root6}
    rows = []
    for line in output.splitlines():
        c, seed, fun, accepted, fev = line.split()
        rows.append((int(c), int(seed), float(fun), int(accepted), int(fev)))
    if len(rows) != 48:
        raise AssertionError(f"expected 48 Fortran reference cases, got {len(rows)}")
    scipy_pass, fortran_pass, worst_gap, worst_scipy_gap = 0, 0, 0.0, 0.0
    for c, seed, f_fortran, accepted, nfev in rows:
        x0, stepsize, temp = cases[c]
        scipy_result = basinhopping(
            lambda z: benchmark(c, z), np.asarray(x0, dtype=float),
            niter=120, T=temp, stepsize=stepsize, rng=119 * seed + 17,
            minimizer_kwargs={"method": "Nelder-Mead",
                              "options": {"maxiter": 240,
                                          "xatol": 1e-8, "fatol": 1e-10}},
        )
        f_ref = float(scipy_result.fun)
        gap = abs(f_fortran - target[c])
        sgap = abs(f_ref - target[c])
        worst_gap = max(worst_gap, gap)
        worst_scipy_gap = max(worst_scipy_gap, sgap)
        if gap < 1e-5:
            fortran_pass += 1
        if sgap < 1e-5:
            scipy_pass += 1
        if not (0 <= accepted <= 120 and nfev >= 121):
            raise AssertionError(f"invalid counts for case {c}, seed {seed}")
        if not math.isfinite(f_fortran) or f_fortran < target[c] - 1e-7:
            raise AssertionError(f"invalid numerical objective for case {c}, seed {seed}")
    print(f"SciFort versus known global objective: {fortran_pass}/48 within 1e-5")
    print(f"SciPy versus known global objective: {scipy_pass}/48 within 1e-5")
    print(f"Worst absolute SciFort objective gap: {worst_gap:.4g}")
    print(f"Worst absolute SciPy objective gap: {worst_scipy_gap:.4g}")
    if fortran_pass != 48:
        raise AssertionError("SciFort missed an objective benchmark")


if __name__ == "__main__":
    main()
