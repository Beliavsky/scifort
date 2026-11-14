#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Check constructor selection, soft counts, infeasible targets and cyclic equations.

Periodic selection uses the pinned upstream control flow and corrected cyclic
residual evaluator documented in check_generate_knots.py, not an unmodified
SciPy 1.18 binary. Nonperiodic knots use installed SciPy's public generator.
"""
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.interpolate import BSpline, generate_knots
from check_generate_knots import upstream_functions, residuals


def canonical(t, k):
    t = t.copy()
    independent = len(t)-2*k-1
    period = t[-k-1]-t[k]
    for j in range(1, k+1):
        t[k-j] = t[k-j+independent]-period
        t[-k+j-1] = t[k+j]+period
    return t


def main():
    exe = Path(sys.argv[1]).resolve()
    oracle = upstream_functions()
    rng = np.random.default_rng(20261119)
    cases = []
    for k in range(1, 6):
        for periodic in (False, True):
            for scalar in (False, True):
                for mode in range(6):
                    m = 47
                    ndim = 1 if scalar else 3
                    u = np.r_[0., np.cumsum(rng.uniform(.5, 1.5, m-1))]
                    u /= u[-1]
                    x = np.column_stack([np.cos((7+a*4)*u)+.08*rng.normal(size=m)
                                         for a in range(ndim)])
                    if periodic:
                        x[-1] = x[0]
                    w = rng.uniform(.5, 1.5, m)
                    provided = scalar or mode % 2 == 0
                    if not provided:
                        u = np.r_[0., np.cumsum(np.linalg.norm(np.diff(x, axis=0), axis=1))]
                        u /= u[-1]
                    bounds = (0., 1.)
                    capacity = None
                    target = .02*np.sum((w[:, None]*x)**2)
                    if mode == 0:
                        target = 0.
                    elif mode == 1:
                        target = 2*np.sum((w[:, None]*x)**2)
                    elif mode == 3:
                        capacity = 2*(k+1)
                        minimal = np.r_[np.zeros(k+1), np.ones(k+1)]
                        if periodic:
                            minimal = canonical(minimal, k)
                        first = np.insert(minimal, k+1, u[(m+1)//2-1])
                        if periodic:
                            first = canonical(first, k)
                        fp0 = residuals(u, x, minimal, k, w, periodic)[1]
                        fp1 = residuals(u, x, first, k, w, periodic)[1]
                        target = (fp0+fp1)/2
                        assert fp1 < target < fp0
                    elif mode == 4:
                        capacity = 2*(k+1)
                        target = 1e-12
                    elif mode == 5:
                        target = 1e-8
                        if not periodic:
                            bounds = (-.2, 1.3)
                    if target == 0:
                        name = "_periodic_knots" if periodic and k > 1 else "_not_a_knot"
                        selected = oracle[name](u, k)
                    elif periodic:
                        selected = list(oracle["_generate_knots_impl"](
                            u, x, w, *bounds, k, target, capacity, True))[-1]
                        selected = canonical(selected, k)
                    else:
                        selected = list(generate_knots(u, x, w=w, k=k, s=target,
                                                       nest=capacity, xb=bounds[0], xe=bounds[1]))[-1]
                    cases.append((k, periodic, scalar, provided, capacity, target,
                                  bounds, u, x, w, selected, mode))
    rows = []
    q = np.linspace(-.1, 1.1, 29)
    for k, periodic, scalar, provided, capacity, target, bounds, u, x, w, _, _ in cases:
        rows.append(f'{len(u)} {x.shape[1]} {k} {int(periodic)} {int(scalar)} '
                    f'{int(provided)} {capacity or 0} {len(q)}')
        arrays = [[target, *bounds], *x, w]
        if provided:
            arrays.append(u)
        arrays.append(q)
        rows.extend(' '.join(f'{v:.17g}' for v in row) for row in arrays)
    result = subprocess.run([str(exe)], input='\n'.join(rows)+'\n', text=True,
                            capture_output=True, check=True)
    lines = iter(result.stdout.splitlines())
    successes = failures = soft = comparisons = 0
    residual_max = stationarity_max = representation_max = 0.
    for index, case in enumerate(cases):
        k, periodic, scalar, provided, capacity, target, bounds, u, x, w, selected, mode = case
        code, valid, allocated = map(int, next(lines).split())
        if mode == 4:
            assert (code, valid, allocated) == (3, 0, 0), (index, code, valid, allocated)
            failures += 1
            continue
        assert (code, valid, allocated) == (0, 1, 1), (index, case[:7], code)
        successes += 1
        nt = int(next(lines))
        sites = np.fromstring(next(lines), sep=' ')
        t = np.fromstring(next(lines), sep=' ')
        ndim = x.shape[1]
        c = np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(ndim)])
        values = np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(ndim)])
        queries = [np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(ndim)])
                   for _ in range(k+1)]
        np.testing.assert_allclose(sites, u, rtol=2e-14, atol=2e-15)
        assert nt == len(selected), (index, nt, len(selected))
        np.testing.assert_allclose(t, selected, rtol=2e-14, atol=2e-14, err_msg=str(index))
        comparisons += nt
        if mode == 3:
            assert nt > capacity, (index, 'soft count')
            soft += 1
        bs = BSpline(t, c, k)
        np.testing.assert_allclose(values, bs(sites), rtol=3e-10, atol=3e-10)
        for d, actual in enumerate(queries):
            expected = bs(q, nu=d, extrapolate='periodic' if target == 0 and periodic and k > 1 else True)
            error = np.max(abs(actual-expected)/(1+abs(expected)))
            assert error < 3e-8, (index, 'representation', d, error)
            representation_max = max(representation_max, error)
        active = len(u)-int(periodic)
        fp = np.sum((w[:active, None]*(values[:active]-x[:active]))**2)
        if target:
            assert fp <= 1.0011*target, (index, 'target', fp, target)
            if nt > 2*(k+1):
                residual_max = max(residual_max, abs(fp-target)/target)
                assert abs(fp-target) <= .0011*target, (index, 'active target', fp, target)
        else:
            np.testing.assert_allclose(values, x, rtol=2e-11, atol=2e-11)
        if periodic:
            for d in range(k):
                seam = bs([u[0], u[-1]], nu=d)
                scale = max(1., np.max(abs(bs.derivative(d).c)))
                assert np.max(abs(seam[1]-seam[0])) < 4096*np.finfo(float).eps*scale
        if target:
            nc = len(c)
            ni = nc-k if periodic else nc
            wrap = np.eye(ni)[np.arange(nc) % ni]
            basis = BSpline(t, wrap, k)
            a = basis(sites[:active])*w[:active, None]
            rhs = x[:active]*w[:active, None]
            interior = t[k+1:-k-1]
            if len(interior) == 0:
                gradient = a.T@(a@c[:ni]-rhs)
                denom = np.linalg.norm(a)*(np.linalg.norm(a@c[:ni])+np.linalg.norm(rhs))
                error = np.linalg.norm(gradient)/max(denom, 1e-300)
                assert error < 2e-12, (index, 'polynomial/constant normal equations', error)
                stationarity_max = max(stationarity_max, error)
                continue
            jumps = (basis(np.nextafter(interior, np.inf), nu=k)-
                     basis(np.nextafter(interior, -np.inf), nu=k))
            jumps /= np.max(abs(jumps))
            gradient = a.T@(a@c[:ni]-rhs)
            penalty = jumps.T@(jumps@c[:ni])
            lam = max(0., -np.sum(gradient*penalty)/np.sum(penalty**2))
            denom = (np.linalg.norm(a)*np.linalg.norm(a@c[:ni]-rhs)+
                     lam*np.linalg.norm(jumps)**2*np.linalg.norm(c[:ni]))
            error = np.linalg.norm(gradient+lam*penalty)/max(denom, 1e-300)
            assert error < 2e-8, (index, 'shared-multiplier normal equations', error)
            stationarity_max = max(stationarity_max, error)
    assert next(lines, None) is None
    print(f'PASS: {successes} modern fits; {failures} explicit infeasible targets; '
          f'{soft} feasible soft-count overshoots; {comparisons} selected-knot comparisons')
    print(f'Maximum residual relative error {residual_max:.6g}; '
          f'scaled stationarity {stationarity_max:.6g}; representation {representation_max:.6g}')


if __name__ == '__main__':
    main()
