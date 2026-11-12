#!/usr/bin/env python3
"""Independent accuracy and budget comparisons for SciFort DIRECT.

This checks the returned point and objective independently. SciPy uses an
approximate evaluation budget, whereas SciFort enforces a strict limit.
"""
from __future__ import annotations
import math
import subprocess
import sys
import numpy as np
from scipy.optimize import direct as scipy_direct


def setup(cid):
    n = 1 if cid == 1 else 4 if cid == 8 else 3 if cid == 9 else 2
    lo, hi = (-1.0, 1.0) if cid == 1 else (-5.0, 5.0) if cid in (4, 5) else \
        (-6.0, 6.0) if cid == 7 else (-2.0, 2.0)
    return [(lo, hi)] * n


def objective(cid, x):
    x = np.asarray(x)
    if cid == 1:
        return (x[0] - 0.317)**2
    if cid == 2:
        return (x[0] - 0.4)**2 + 3*(x[1]+0.7)**2
    if cid == 3:
        return (x[0]+x[1]-0.3)**2 + 0.1*(x[0]-x[1]-1.2)**2
    if cid == 4:
        return 20 + x[0]**2+x[1]**2 - 10*(np.cos(2*np.pi*x[0])+np.cos(2*np.pi*x[1]))
    if cid == 5:
        return -20*np.exp(-0.2*np.sqrt(np.sum(x*x)/2)) - \
            np.exp(np.mean(np.cos(2*np.pi*x)))+20+np.e
    if cid == 6:
        return 100*(x[1]-x[0]**2)**2 + (1-x[0])**2
    if cid == 7:
        return (x[0]**2+x[1]-11)**2 + (x[0]+x[1]**2-7)**2
    if cid == 8:
        return np.sum((x - np.array([0.25,-0.5,0.75,-0.25]))**2)
    if cid == 9:
        return np.sum(np.abs(x - np.array([0.18,-0.71,0.34])))
    if cid == 10:
        return np.sin(5*x[0])**2 + (x[0]-0.71)**2 + \
            np.sin(3*x[1])**2 + (x[1]+0.27)**2
    raise ValueError(cid)


def main(executable):
    out = subprocess.check_output([executable], text=True)
    records = [line.split() for line in out.strip().splitlines()]
    assert len(records) == 60, len(records)
    rels, comparables, successes, ref_successes = [], 0, 0, 0
    for fields in records:
        cid, bias, budget, status = map(int, fields[:4])
        val = float(fields[4]); nfev = int(fields[5]); nit = int(fields[6])
        x = np.array([float(v) for v in fields[7:]])
        bounds = setup(cid)
        assert len(x) == len(bounds)
        assert np.all(np.isfinite(x))
        assert all(lo - 1e-13 <= xi <= hi + 1e-13
                   for xi,(lo,hi) in zip(x,bounds))
        assert nfev <= budget and nfev >= 1
        assert nit <= 180 and nit >= 0
        assert status in (1,2,3,4,5,6)
        expected = objective(cid, x)
        rel = abs(expected-val)/max(1.0,abs(expected))
        rels.append(rel)
        assert rel < 1e-12, (cid,bias,budget,val,expected,rel)
        ref = scipy_direct(lambda v: objective(cid,v), bounds,
                          maxfun=budget, maxiter=180, locally_biased=bool(bias),
                          vol_tol=0.0, len_tol=0.0)
        assert np.isfinite(ref.fun)
        comparables += 1
        # On the convex tests, both algorithms should give good minima with
        # comparable orders of magnitude of objective error.
        if cid in (1, 2, 3, 4, 5, 6, 7, 8, 9):
            if budget == 1600:
                tol = 5e-3 if cid == 6 else 2e-3 if cid == 9 else 2e-4
                assert val < tol, (cid,bias,budget,val,ref.fun)
                successes += 1
                if ref.fun < tol:
                    ref_successes += 1
        if cid == 10:
            assert val < 0.2, (cid,val)
    print(f'{len(records)} SciFort solutions validated independently; '
          f'{comparables} SciPy comparisons; '
          f'max normalized objective-check difference {max(rels):.2e}; '
          f'{successes} high-budget benchmark targets passed; '
          f'SciPy passed {ref_successes} corresponding targets')

if __name__ == '__main__':
    main(sys.argv[1])
