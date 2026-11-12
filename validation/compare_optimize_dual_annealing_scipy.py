#!/usr/bin/env python3
"""Cross-check independent Fortran dual annealing and SciPy 1.17 benchmarks.

Check objective values exactly; compare final objective quality to a separate
SciPy run, not RNG streams or paths (the generators are deliberately different).
"""
import math
import subprocess
import sys

import numpy as np
from scipy.optimize import dual_annealing

ITERATIONS = 360
NCASES = 8
NSEEDS = 6


def setup(case):
    if case == 7:
        return np.tile([-5.0, 5.0], (3, 1))
    if case == 8:
        return np.array([[-4.0, 4.0]])
    bounds = np.tile([-5.0, 5.0], (2, 1))
    if case == 3:
        bounds[:] = [-5.12, 5.12]
    elif case == 4:
        bounds[:] = [-3.0, 3.0]
    elif case == 5:
        bounds[:] = [-8.0, 8.0]
    elif case == 6:
        bounds[:] = [-6.0, 6.0]
    return bounds


def objective(case, x):
    x = np.asarray(x, dtype=float)
    if case == 1:
        return np.sum((x - [1.2, -0.7])**2)
    if case == 2:
        return (-20.0 * np.exp(-0.2 * np.sqrt(np.mean(x*x)))
                - np.exp(np.mean(np.cos(2*np.pi*x))) + 20.0 + np.e)
    if case == 3:
        return 20.0 + np.sum(x*x - 10*np.cos(2*np.pi*x))
    if case == 4:
        return 100*(x[1] - x[0]**2)**2 + (1-x[0])**2
    if case == 5:
        return 1.0 + np.sum(x*x)/4000.0 - np.prod(np.cos(x/np.sqrt(np.arange(1,len(x)+1))))
    if case == 6:
        return (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2
    if case == 7:
        return np.sum((x - [0.5,-1.0,2.0])**2)
    if case == 8:
        return (x[0]-0.25)**2 + 0.3*(1-np.cos(8*x[0]-2))
    raise ValueError(case)


def main():
    output = subprocess.check_output([sys.argv[1]], text=True)
    rows = [line.split() for line in output.splitlines() if line.startswith('CASE')]
    assert len(rows) == NCASES*NSEEDS, len(rows)
    ff = {i:[] for i in range(1,NCASES+1)}
    sf = {i:[] for i in range(1,NCASES+1)}
    max_recalculation = 0.0
    for row in rows:
        _, casetxt, seedtxt, statustxt, nfevtxt, ftxt, *xtxt = row
        case, seedno, status, nfev = map(int, [casetxt,seedtxt,statustxt,nfevtxt])
        n = len(setup(case))
        x = np.array(list(map(float, xtxt)))
        f = float(ftxt)
        bounds = setup(case)
        assert status == 0, (case, seedno, status)
        assert nfev == 1 + 2*n*ITERATIONS, (case, seedno, nfev)
        assert x.shape == (n,)
        assert np.all(np.isfinite(x))
        assert np.all(x >= bounds[:,0] - 1e-12)
        assert np.all(x <= bounds[:,1] + 1e-12)
        assert np.isfinite(f)
        difference = abs(f - objective(case, x))
        max_recalculation = max(max_recalculation, difference)
        assert difference < 2e-11 * (1 + abs(f)), (case, seedno, difference)
        ff[case].append(f)
        reference = dual_annealing(lambda z: objective(case,z), bounds, maxiter=ITERATIONS,
                                   no_local_search=True, rng=1103+seedno*137+case)
        assert np.isfinite(reference.fun)
        assert np.all(reference.x <= bounds[:,1]+1e-12)
        assert np.all(reference.x >= bounds[:,0]-1e-12)
        sf[case].append(float(reference.fun))

    # All tested functions have known global minimum zero. These generously
    # scaled thresholds check actual global-search quality, not just feasibility.
    quality_bounds = {1:0.01, 2:0.2, 3:0.5, 4:0.5, 5:0.05, 6:0.1, 7:0.05, 8:0.01}
    for case in range(1,NCASES+1):
        fval = np.array(ff[case]); sval = np.array(sf[case])
        assert np.max(fval) < quality_bounds[case], (case, fval)
        assert np.max(sval) < quality_bounds[case], (case, sval)
        print(f'case {case} SciFort median {np.median(fval):.6g}, SciPy median {np.median(sval):.6g}, '
              f'SciFort max {np.max(fval):.6g}, SciPy max {np.max(sval):.6g}')
    print(f'DUAL_ANNEALING_REFERENCE_CASES {len(rows)}')
    print(f'DUAL_ANNEALING_MAX_RECALCULATION {max_recalculation:.3e}')

if __name__ == '__main__':
    main()
