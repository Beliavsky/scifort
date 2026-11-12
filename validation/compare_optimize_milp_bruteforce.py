#!/usr/bin/env python3
"""Exhaustively verify small bounded pure-integer MILPs against enumeration.

This checks optimality and infeasibility without using SciPy's MILP solver.
"""
from __future__ import annotations

from itertools import product
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from compare_optimize_milp_scipy import encode


def main():
    rng = np.random.default_rng(261110)
    models = []
    objectives = []
    for _ in range(220):
        n = int(rng.integers(1, 6))
        m = int(rng.integers(0, 5))
        lower = rng.integers(-2, 2, size=n).astype(float)
        upper = lower + rng.integers(1, 5, size=n)
        c = rng.integers(-6, 7, size=n).astype(float)
        a = rng.integers(-4, 5, size=(m, n)).astype(float)
        b_l = rng.integers(-10, 7, size=m).astype(float)
        b_u = rng.integers(-3, 12, size=m).astype(float)
        # Include both equality and one-sided rows, and some infeasible cases.
        for row in range(m):
            if b_l[row] > b_u[row]:
                b_l[row], b_u[row] = b_u[row], b_l[row]
            kind = int(rng.integers(0, 5))
            if kind == 0:
                b_u[row] = b_l[row]
            elif kind == 1:
                b_l[row] = -np.inf
            elif kind == 2:
                b_u[row] = np.inf
        optimal = np.inf
        for values in product(*(range(int(l), int(u)+1) for l,u in zip(lower,upper))):
            x = np.asarray(values, dtype=float)
            v = a @ x
            if np.all(v >= b_l-1e-10) and np.all(v <= b_u+1e-10):
                optimal = min(optimal, float(c @ x))
        models.append((c, np.ones(n, dtype=np.int32), lower, upper, a, b_l, b_u))
        objectives.append(optimal)

    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix='scifort-milp-exhaust-') as td:
        td = Path(td)
        sources = [root / 'src/kinds.f90', root / 'src/optimize/linprog.f90',
                   root / 'src/optimize/milp.f90', root / 'validation/optimize_milp_driver.f90']
        subprocess.run(['gfortran', '-std=f2018', '-O3', '-fcheck=all',
                        *map(str, sources), '-o', str(td / 'driver')],
                       check=True, cwd=td, capture_output=True, text=True)
        proc = subprocess.run([str(td / 'driver')], input=encode(models),
                              capture_output=True, text=True, cwd=td, check=True)
    results = proc.stdout.strip().splitlines()
    assert len(results) == 2*len(models)
    checks = 0
    for i, optimum in enumerate(objectives):
        split = results[2*i].split()
        status, fun = int(split[1]), float(split[2])
        if np.isfinite(optimum):
            assert status == 0, (i, 'status', status, optimum, models[i])
            assert abs(fun-optimum)<1e-7, (i, 'objective', fun, optimum, models[i])
            checks += 2
        else:
            assert status == 2, (i, 'infeasible status', status, models[i])
            checks += 1
    print(f'MILP exhaustive integer tests passed: {len(models)} problems, {checks} checks')


if __name__ == '__main__':
    main()
