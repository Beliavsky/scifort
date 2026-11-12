"""Cross-check SciFort's bounded simplicial optimizer against SciPy SHGO.

The Fortran validation driver deliberately uses a small independent coordinate
pattern search, rather than the standard SciFort Nelder-Mead adapter, to isolate
the global simplicial implementation from SciFort's local solver backend.
Neither evaluation histories nor the ``n`` sampling semantics are expected to
match SciPy's adaptively constructed complex.
"""
from __future__ import annotations

import math
import subprocess
import sys

import numpy as np
from scipy.optimize import shgo as scipy_shgo


def problem(cid):
    dim = 1 if cid in (1, 7, 10) else 3 if cid == 8 else 2
    bound = (-5.12, 5.12) if cid == 4 else (-5.0, 5.0) if cid == 5 else \
        (-6.0, 6.0) if cid == 6 else (-2.0, 2.0)
    return [bound] * dim


def objective(cid, x):
    x = np.asarray(x, dtype=float)
    if cid == 1:
        return (x[0]-0.317)**2
    if cid == 2:
        return (x[0]-0.43)**2 + 3*(x[1]+0.61)**2
    if cid == 3:
        return (x[0]+x[1]-0.3)**2 + 0.1*(x[0]-x[1]-1.2)**2
    if cid == 4:
        return 20 + np.sum(x*x) - 10*np.sum(np.cos(2*np.pi*x))
    if cid == 5:
        return -20*np.exp(-0.2*np.sqrt(np.mean(x*x))) - \
            np.exp(np.mean(np.cos(2*np.pi*x))) + 20 + math.e
    if cid == 6:
        return (x[0]**2+x[1]-11)**2 + (x[0]+x[1]**2-7)**2
    if cid == 7:
        return (x[0]**2-1)**2 + 0.1*x[0]
    if cid == 8:
        return np.sum((x-np.array([0.25,-0.5,0.75]))**2)
    if cid == 9:
        return np.sin(5*x[0])**2+(x[0]-0.71)**2 + \
            np.sin(3*x[1])**2+(x[1]+0.27)**2
    if cid == 10:
        return abs(x[0]-0.123)
    if cid == 11:
        return min((x[0]-1)**2+(x[1]-1)**2,
                   (x[0]+1)**2+(x[1]+1)**2+0.25)
    raise ValueError(cid)


def main(executable):
    lines = subprocess.check_output([executable], text=True).splitlines()
    assert len(lines) == 33, len(lines)
    ncomparisons = 0
    nmatched = 0
    max_objective_error = 0.0
    for line in lines:
        cells = line.split()
        cid, iters, status = map(int, cells[:3])
        reported = float(cells[3])
        nfev, nlfev, nit, nlocal, nminima = map(int, cells[4:9])
        x = np.array([float(c) for c in cells[9:]])
        bounds = problem(cid)
        assert len(x) == len(bounds)
        assert np.all(np.isfinite(x))
        assert all(lo-1e-12 <= xi <= hi+1e-12
                   for xi, (lo, hi) in zip(x,bounds))
        assert status == 0 and nit == iters, (cid, iters, status, nit)
        assert 0 < nfev <= 10000 and 0 <= nlfev <= nfev
        assert nlocal >= nminima >= 1
        value = objective(cid,x)
        error = abs(reported-value)/max(1.0,abs(value))
        max_objective_error = max(max_objective_error,error)
        assert error < 2e-12, (cid,iters,error)
        reference = scipy_shgo(lambda x: objective(cid,x),bounds,iters=iters,
                               sampling_method='simplicial')
        assert np.isfinite(reference.fun), (cid,iters,reference.message)
        ncomparisons += 1
        # True optima of most benchmark cases are near zero. The tilted
        # double-well has a negative known minimum of about -0.100617.
        if cid == 7:
            assert reported < -0.100, (cid,iters,reported)
        elif cid == 9:
            assert reported < 0.10, (cid,iters,reported)
        elif cid == 11 and iters == 1:
            assert reported <= 0.25 + 1e-8, (cid,iters,reported)
        else:
            assert reported <= 1e-6, (cid,iters,reported)
        if iters >= 2:
            # Our independently polished results should be at least comparable
            # to SciPy, allowing a small tolerance for the local solvers.
            assert reported <= reference.fun+1e-4, \
                (cid,iters,reported,reference.fun)
            nmatched += 1
    print(f'{ncomparisons} independent SciPy SHGO comparisons; '
          f'{nmatched} refined-value comparisons; '
          f'max normalized objective recomputation error '
          f'{max_objective_error:.2e}')


if __name__ == '__main__':
    main(sys.argv[1])
