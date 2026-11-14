#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Large native selection against SciPy banded QR and independent sparse LSQR.

Uses pinned upstream selection control flow with canonical periodic exteriors,
as documented in check_generate_knots.py. No dense sample-by-coefficient matrix
is formed by this checker. Peak workspace estimates are reported separately
from elapsed time; they are not measurements of process resident memory.
"""
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.interpolate import BSpline, generate_knots, make_lsq_spline
from scipy.sparse import csr_matrix, diags, eye
from scipy.sparse.linalg import lsqr, norm
from check_generate_knots import upstream_functions, native
from check_modern_selection import canonical


def sparse_design(u, t, k, periodic):
    active = len(u)-int(periodic)
    b = BSpline.design_matrix(u[:active], t, k).tocsr()
    if not periodic:
        return b
    independent = len(t)-2*k-1
    folded = csr_matrix((b.data, b.indices % independent, b.indptr),
                        shape=(active, independent))
    folded.sum_duplicates()
    folded.sort_indices()
    return folded


def sparse_residuals(u, y, t, k, w, periodic=False):
    if not periodic:
        values = make_lsq_spline(u, y, t, k=k, w=w)(u)
    else:
        t = canonical(t, k)
        a = sparse_design(u, t, k, True).multiply(w[:-1, None]).tocsr()
        rhs = y[:-1]*w[:-1, None]
        c = np.empty((a.shape[1], y.shape[1]))
        for j in range(y.shape[1]):
            solution = lsqr(a, rhs[:, j], atol=2e-15, btol=2e-15,
                            iter_lim=max(1000, 4*a.shape[1]))
            assert solution[1] in (0, 1, 2), ("sparse reference did not converge", solution[1])
            c[:, j] = solution[0]
        values = a@c/w[:-1, None]
    active = len(u)-int(periodic)
    errors = np.zeros(len(u))
    errors[:active] = np.sum(((values-y[:active])*w[:active, None])**2, axis=1)
    return errors, errors.sum()


def jumps_matrix(t, k, periodic):
    nc = len(t)-k-1
    derivative = eye(nc, format='csr')
    degree = k
    knots = t.copy()
    while degree:
        rows = derivative.shape[0]-1
        factor = degree/(knots[degree+1:degree+1+rows]-knots[1:1+rows])
        derivative = diags([-factor, factor], [0, 1], shape=(rows, rows+1))@derivative
        knots = knots[1:-1]
        degree -= 1
    jump = (derivative[1:]-derivative[:-1]).tocsr()
    if periodic:
        ni = nc-k
        wrap = csr_matrix((np.ones(nc), (np.arange(nc), np.arange(nc) % ni)), shape=(nc, ni))
        jump = jump@wrap
    if jump.nnz:
        jump /= np.max(abs(jump.data))
    return jump


def main():
    gen_exe, fit_exe = map(lambda s: Path(s).resolve(), sys.argv[1:3])
    oracle = upstream_functions()
    oracle['_get_residuals'] = sparse_residuals
    rng = np.random.default_rng(20261120)
    cases = []
    comparisons = sequence_cases = 0
    native_time = 0.
    for k in range(1, 6):
        for m in (301, 1001, 10000):
            if m == 10000 and k != 3:
                continue
            for periodic in (False, True):
                for ndim in (1, 3):
                    for soft in (False, True):
                        if soft and m != 1001:
                            continue
                        u = np.r_[0., np.cumsum(rng.uniform(.8, 1.2, m-1))]
                        u /= u[-1]
                        x = np.column_stack([np.cos((12+j*3)*u)+.01*rng.normal(size=m)
                                             for j in range(ndim)])
                        if periodic:
                            x[-1] = x[0]
                        w = rng.uniform(.7, 1.3, m)
                        target = .002*np.sum((w[:, None]*x)**2)
                        capacity = None
                        if soft:
                            capacity = 2*(k+1)
                            minimal = np.r_[np.zeros(k+1), np.ones(k+1)]
                            if periodic:
                                minimal = canonical(minimal, k)
                            index = (m+1)//2-1 if periodic else m//2
                            first = np.insert(minimal, k+1, u[index])
                            if periodic:
                                first = canonical(first, k)
                            design = sparse_design(u, first, k, periodic)
                            x = design@rng.normal(size=(design.shape[1], ndim))
                            x += .001*rng.normal(size=x.shape)
                            if periodic:
                                x = np.vstack([x, x[0]])
                            target = .5*(sparse_residuals(u, x, minimal, k, w, periodic)[1]+
                                         sparse_residuals(u, x, first, k, w, periodic)[1])
                        if periodic:
                            expected = list(oracle['_generate_knots_impl'](
                                u, x, w, 0., 1., k, target, capacity, True))
                            expected[-1] = canonical(expected[-1], k)
                        else:
                            expected = list(generate_knots(u, x, w=w, k=k, s=target, nest=capacity))
                        start = time.perf_counter()
                        actual = native(gen_exe, u, x, k, target, w, periodic, capacity, (0., 1.))
                        native_time += time.perf_counter()-start
                        assert len(actual) == len(expected), (k, m, periodic, ndim, soft, 'sequence length')
                        for a, b in zip(actual, expected):
                            np.testing.assert_allclose(a, b, rtol=2e-14, atol=2e-14)
                            comparisons += len(a)
                        sequence_cases += 1
                        cases.append((u, x, w, k, target, periodic, capacity, actual[-1]))
    # Large interpolation vectors need no collocation solve in the generator.
    for k in range(1, 6):
        for periodic in (False, True):
            u = np.linspace(0., 1., 10000)
            x = np.sin(2*np.pi*u)[:, None]
            x[-1] = x[0]
            actual = native(gen_exe, u, x, k, 0., None, periodic, None, (0., 1.))
            name = '_periodic_knots' if periodic and k > 1 else '_not_a_knot'
            assert len(actual) == 1
            np.testing.assert_allclose(actual[0], oracle[name](u, k), rtol=2e-14, atol=2e-14)
            comparisons += len(actual[0])
            sequence_cases += 1
    q = np.linspace(-.05, 1.05, 31)
    rows = []
    for u, x, w, k, target, periodic, capacity, _ in cases:
        rows.append(f'{len(u)} {x.shape[1]} {k} {int(periodic)} {int(x.shape[1] == 1)} '
                    f'1 {capacity or 0} {len(q)}')
        rows.extend(' '.join(f'{v:.17g}' for v in row)
                    for row in [[target, 0., 1.], *x, w, u, q])
    start = time.perf_counter()
    proc = subprocess.run([str(fit_exe)], input='\n'.join(rows)+'\n', text=True,
                          capture_output=True, check=True)
    fit_time = time.perf_counter()-start
    lines = iter(proc.stdout.splitlines())
    residual_max = stationarity_max = seam_max = 0.
    soft_cases = 0
    for index, (u, x, w, k, target, periodic, capacity, selected) in enumerate(cases):
        assert tuple(map(int, next(lines).split())) == (0, 1, 1), index
        nt = int(next(lines))
        sites = np.fromstring(next(lines), sep=' ')
        t = np.fromstring(next(lines), sep=' ')
        nc = nt-k-1
        c = np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(x.shape[1])])
        values = np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(x.shape[1])])
        queries = [np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(x.shape[1])])
                   for _ in range(k+1)]
        np.testing.assert_array_equal(sites, u)
        np.testing.assert_allclose(t, selected, rtol=2e-14, atol=2e-14)
        if capacity:
            assert nt > capacity, (index, k, len(u), periodic, x.shape[1], nt, capacity, target)
            soft_cases += 1
        active = len(u)-int(periodic)
        fp = np.sum((w[:active, None]*(values[:active]-x[:active]))**2)
        error = abs(fp-target)/target
        assert error < .0011, (index, 'residual', error)
        residual_max = max(residual_max, error)
        bs = BSpline(t, c, k)
        np.testing.assert_allclose(values, bs(u), rtol=3e-10, atol=3e-10)
        for d, actual in enumerate(queries):
            np.testing.assert_allclose(actual, bs(q, nu=d), rtol=3e-9, atol=3e-8)
        a = sparse_design(u, t, k, periodic).multiply(w[:active, None]).tocsr()
        rhs = x[:active]*w[:active, None]
        ni = nc-k if periodic else nc
        if periodic:
            np.testing.assert_allclose(c, c[np.arange(nc) % ni], rtol=3e-11, atol=3e-12)
            for d in range(k):
                seam = bs([0., 1.], nu=d)
                scale = max(1., np.max(abs(bs.derivative(d).c)))
                err = np.max(abs(seam[1]-seam[0]))/scale
                assert err < 4096*np.finfo(float).eps, (index, 'seam', d, err)
                seam_max = max(seam_max, err)
        jump = jumps_matrix(t, k, periodic)
        gradient = a.T@(a@c[:ni]-rhs)
        penalty = jump.T@(jump@c[:ni])
        lam = max(0., -np.sum(gradient*penalty)/np.sum(penalty**2))
        denom = norm(a)*np.linalg.norm(a@c[:ni]-rhs)+lam*norm(jump)**2*np.linalg.norm(c[:ni])
        err = np.linalg.norm(gradient+lam*penalty)/max(denom, 1e-300)
        assert err < 2e-8, (index, 'common multiplier', err)
        stationarity_max = max(stationarity_max, err)
    assert next(lines, None) is None
    print(f'PASS: {sequence_cases} large sequences; {comparisons} knots; {len(cases)} fits; '
          f'{soft_cases} successful soft-count overshoots')
    print(f'Maximum residual relative error {residual_max:.6g}; scaled stationarity '
          f'{stationarity_max:.6g}; unwrapped seam {seam_max:.6g}')
    print(f'Native generator process time {native_time:.3f} s; constructor process time {fit_time:.3f} s')
    m, k, channels = 10000, 3, 1
    nt = m+2*k
    workspace = 8*(m*(channels+k+2)+nt*(5*k+8+2*channels)+(nt-k-1)*channels)+4*nt
    dense = 8*m*(nt-k-1)
    print(f'10000-site cyclic kernel declared buffers: {workspace/2**20:.3f} MiB; '
          f'equivalent full dense design alone: {dense/2**20:.3f} MiB')


if __name__ == '__main__':
    main()

