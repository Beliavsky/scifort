#!/usr/bin/env python3
"""Reference comparison for SciFort differential evolution.

Rebuild and run from the project root with
    python test/compare_optimize_differential_evolution_scipy.py

Requires NumPy, SciPy and gfortran. Results are checked by known objective
values, exact constraints, reproducibility and approximate performance;
random streams differ, so solver trajectories must NOT be compared directly.
"""
from __future__ import annotations

import argparse
import math
import pathlib
import subprocess
import tempfile

import numpy as np
from scipy.optimize import differential_evolution

ROOT = pathlib.Path(__file__).resolve().parents[1]
STRATEGIES = [
    "best1bin", "rand1bin", "best2bin", "rand2bin",
    "randtobest1bin", "currenttobest1bin", "best1exp", "rand1exp",
    "best2exp", "rand2exp", "randtobest1exp", "currenttobest1exp",
]


def problem(case: int):
    if case in (1, 2, 3, 4, 5):
        bounds = [(-5.0, 5.0)] * 2
        if case == 2:
            bounds = [(-2.0, 2.0)] * 2
        if case == 3:
            bounds = [(-5.12, 5.12)] * 2
        if case == 5:
            bounds = [(-1.0, 2.0), (-1.0, 1.0)]
    else:
        bounds = [(-5.0, 5.0), (-3.0, 4.0), (2.0, 2.0)]

    def fun(xx):
        x = np.asarray(xx)
        if case in (1, 6):
            f = (x[0] - .75) ** 2 + (x[1] + 1.1) ** 2
            return float(f + ((x[2] - 2.0) ** 2 if case == 6 else 0.0))
        if case == 2:
            return float(100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2)
        if case == 3:
            return float(20 + np.sum(x * x - 10 * np.cos(2 * np.pi * x)))
        if case == 4:
            return float(-20 * np.exp(-.2 * np.sqrt(np.mean(x * x)))
                         - np.exp(np.mean(np.cos(2 * np.pi * x))) + 20 + math.e)
        return float((x[0] - 2) ** 2 + (x[1] + 1) ** 2)

    return bounds, fun


def run(compiler: str, flags: list[str]):
    # About 60 combinations, across all SciPy strategy families and test functions.
    cases = []
    for case in range(1, 7):
        for index in range(1, 13):
            seed = 300 + 19 * case + index
            update = 1 if index % 2 else 2
            init = 1 if index % 3 else 2
            cases.append((case, index, seed, update, init))

    with tempfile.TemporaryDirectory(prefix="scifort_de_compare_") as tmp:
        exe = pathlib.Path(tmp) / "de_reference"
        cmd = [compiler, "-std=f2018", "-J", tmp, "-I", tmp, *flags, "-o", str(exe),
               str(ROOT / "src/kinds.f90"),
               str(ROOT / "src/optimize/differential_evolution.f90"),
               str(ROOT / "tools/reference_differential_evolution_driver.f90")]
        subprocess.run(cmd, check=True)
        text = "".join("%d %d %d %d %d\n" % c for c in cases)
        out = subprocess.run([str(exe)], input=text, text=True, check=True,
                             capture_output=True).stdout.splitlines()
    if len(out) != len(cases):
        raise AssertionError(f"expected {len(cases)} results, got {len(out)}")
    scipy_reference = []
    for case, index, seed, update, init in cases:
        bounds, fun = problem(case)
        result = differential_evolution(
            fun, bounds, strategy=STRATEGIES[index-1], maxiter=170, popsize=9,
            tol=1e-3, atol=1e-8, mutation=(.5, 1.), recombination=.7,
            rng=np.random.default_rng(seed), polish=False,
            updating="immediate" if update == 1 else "deferred",
            init="latinhypercube" if init == 1 else "random")
        scipy_reference.append(result)
    passed = 0
    near_global = 0
    scipy_near_global = 0
    non_worse = 0
    max_validation_error = 0.0
    worst_case = None
    for ix, (params, line, ref) in enumerate(zip(cases, out, scipy_reference)):
        vals = line.split()
        ident = tuple(int(v) for v in vals[:5])
        if ident != params:
            raise AssertionError(f"driver case mismatch at {ix}: {ident} != {params}")
        f_value = float(vals[5])
        x = np.array([float(v) for v in vals[6:9]])
        nit, nfev, status = [int(v) for v in vals[9:12]]
        bounds, fun = problem(params[0])
        x = x[:len(bounds)]
        low, high = np.array(bounds).T
        if not np.all(np.isfinite(x)) or np.any(x < low - 1e-12) or np.any(x > high + 1e-12):
            raise AssertionError(f"out of bounds case {ix}: {x}")
        truth = fun(x)
        error = abs(truth - f_value)
        max_validation_error = max(max_validation_error, error)
        if error > 3e-10 * max(1.0, abs(f_value)):
            raise AssertionError(f"incorrect objective case {ix}: {f_value=} {truth=}")
        if not (0 <= nit <= 170 and 0 < nfev <= (171 * 27)) or status not in (0, 1):
            raise AssertionError(f"invalid metadata case {ix}: {nit=} {nfev=} {status=}")
        if f_value < 1e-3:
            near_global += 1
        if ref.fun < 1e-3:
            scipy_near_global += 1
        if f_value <= ref.fun + 1e-3:
            non_worse += 1
        if worst_case is None or (f_value - ref.fun) > worst_case[0]:
            worst_case = (f_value - ref.fun, params, f_value, float(ref.fun))
        passed += 1
    print(f"VALIDATED {passed} SciFort outputs against independent objective and bounds")
    print(f"SciFort near-global (<1e-3): {near_global}/{passed}")
    print(f"SciPy near-global (<1e-3): {scipy_near_global}/{passed}")
    print(f"SciFort equal or better than SciPy (tolerance 1e-3): {non_worse}/{passed}")
    print(f"Max objective reproduction difference: {max_validation_error:.3g}")
    print(f"Largest SciFort-minus-SciPy objective difference: {worst_case}")
    if near_global < 48:
        raise AssertionError("global optimization quality unexpectedly poor")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fc", default="gfortran")
    parser.add_argument("--mode", choices=["checked", "optimized"], default="checked")
    args = parser.parse_args()
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if args.mode == "checked" else ["-O3", "-march=native"]
    run(args.fc, flags)
