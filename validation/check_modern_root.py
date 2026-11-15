# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Compare feasible supplied-knot smoothing near the least-squares limit.

Execute hash-verified SciPy v1.18.0 rational-root control flow with independent
NumPy augmented least-squares solves. No upstream source is patched. Cached
source and its retained BSD license are documented in CODE_PROVENANCE.md.
"""
import ast
import hashlib
from pathlib import Path
import subprocess
import sys
import warnings

import numpy as np
from scipy.interpolate import BSpline, make_splrep, make_splprep


def root_functions():
    path = Path('build/generate-knots-upstream/_fitpack_repro.py')
    raw = path.read_bytes()
    digest = '82449ae292bb361675439b422cec6c084f9fbfa6519eabec06b1009c51a818bd'
    assert hashlib.sha256(raw).hexdigest() == digest, 'unexpected upstream source'
    namespace = dict(np=np, warnings=warnings, MAXIT=20,
                     _iermesg={i: str(i) for i in (1, 2, 3)})
    for node in ast.parse(raw).body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in {
            'Bunch', 'root_rati', 'fprati', 'disc', 'prodd'
        }:
            exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


def reference(functions, u, y, w, t, k, target, periodic):
    active = len(u)-int(periodic)
    nc = len(t)-k-1
    ni = nc-k if periodic else nc
    wrap = np.eye(ni)[np.arange(nc) % ni]
    a = BSpline(t, wrap, k)(u[:active])*w[:active, None]
    rhs = y[:active]*w[:active, None]
    packed, _, _ = functions['disc'](t, k)
    jumps = np.zeros((len(packed), nc))
    for i, row in enumerate(packed):
        jumps[i, i:i+k+2] = row
    d = jumps@wrap
    ls = np.linalg.lstsq(a, rhs, rcond=None)[0]
    fpinf = np.sum((a@ls-rhs)**2)
    poly = np.ones((active, 1)) if periodic else np.vander(u[:active], k+1, increasing=True)
    poly = poly*w[:active, None]
    cp = np.linalg.lstsq(poly, rhs, rcond=None)[0]
    fp0 = np.sum((poly@cp-rhs)**2)
    # The cyclic/banded triangular factors have the same positive QR diagonal.
    p0 = ni/abs(np.diag(np.linalg.qr(a, mode='reduced')[1])).sum()
    state = {}
    history = []
    def evaluate(p):
        aug = np.vstack((a, d/p))
        brhs = np.vstack((rhs, np.zeros((len(d), y.shape[1]))))
        c = np.linalg.lstsq(aug, brhs, rcond=None)[0]
        fp = np.sum((a@c-rhs)**2)
        state['c'] = wrap@c
        history.append((p, fp))
        return fp-target
    result = functions['root_rati'](evaluate, p0, ((0., fp0-target), (np.inf, fpinf-target)),
                                    .001*target)
    assert result.converged, (len(u), k, periodic, result.ier)
    return state['c'], history, a, rhs, d, wrap


def main():
    exe = Path(sys.argv[1]).resolve()
    functions = root_functions()
    rng = np.random.default_rng(20261123)
    cases, rows = [], []
    q = np.linspace(-.03, 1.03, 13)
    for m in (31, 61, 97):
        for k in range(1, 6):
            for periodic in (False, True):
                for scalar, ndim in ((True, 1), (False, 1), (False, 3)):
                    u = np.linspace(0., 1., m)**1.1
                    y = np.column_stack([np.cos((2+j)*2*np.pi*u)
                                         + .08*np.sin((14+2*j)*np.pi*u) for j in range(ndim)])
                    if m != 61:
                        y += .005*rng.normal(size=y.shape)
                    y[-1] = y[0]
                    w = .8+u
                    if periodic:
                        interior = np.linspace(0., 1., 14)
                        t = np.r_[interior[-k-1:-1]-1, interior, interior[1:k+1]+1]
                    else:
                        t = np.r_[np.zeros(k+1), np.linspace(.08, .92, 11), np.ones(k+1)]
                    nc = len(t)-k-1
                    ni = nc-k if periodic else nc
                    active = m-int(periodic)
                    wrap = np.eye(ni)[np.arange(nc) % ni]
                    a = BSpline(t, wrap, k)(u[:active])*w[:active, None]
                    rhs = y[:active]*w[:active, None]
                    ls = np.linalg.lstsq(a, rhs, rcond=None)[0]
                    fpinf = np.sum((a@ls-rhs)**2)
                    for factor in (1.00001, 1.0005, 1.002):
                        target = fpinf*factor
                        rows.append(f'{m} {ndim} {k} {int(scalar)} 1 {len(q)} 0 '
                                    f'{len(t)} {int(periodic)}')
                        for arr in ([target], *y, w, u, q, t):
                            rows.append(' '.join(f'{v:.17g}' for v in arr))
                        cases.append((u, y, w, t, k, target, periodic, scalar))
    proc = subprocess.run([str(exe)], input='\n'.join(rows)+'\n', text=True,
                          capture_output=True, check=True)
    lines = iter(proc.stdout.splitlines())
    coefficient_max = derivative_max = residual_max = stationarity_max = 0.
    public_max = 0.
    public = steps = 0
    for index, (u, y, w, t, k, target, periodic, scalar) in enumerate(cases):
        code = int(next(lines))
        assert code == 0, (index, k, periodic, scalar, code)
        nt = int(next(lines))
        params = np.fromstring(next(lines), sep=' ')
        knots = np.fromstring(next(lines), sep=' ')
        ndim = y.shape[1]
        c = np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(ndim)])
        values = np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(ndim)])
        queries = [np.column_stack([np.fromstring(next(lines), sep=' ') for _ in range(ndim)])
                   for _ in range(min(k, 3)+1)]
        integral = np.fromstring(next(lines), sep=' ')
        assert nt == len(t)
        np.testing.assert_array_equal(knots, t)
        np.testing.assert_array_equal(params, u)
        cr, history, a, rhs, d, wrap = reference(functions, u, y, w, t, k, target, periodic)
        error = np.max(abs(c-cr)/(1+abs(cr)))
        coefficient_max = max(coefficient_max, error)
        assert error < 3e-11, (index, 'coefficients', error)
        bs, ref = BSpline(t, c, k), BSpline(t, cr, k)
        np.testing.assert_allclose(values, bs(u), rtol=3e-10, atol=3e-10)
        for degree, actual in enumerate(queries):
            np.testing.assert_allclose(actual, bs(q, nu=degree), rtol=3e-9, atol=3e-9)
        for degree in range(k+1):
            expected = ref(q, nu=degree)
            error = np.max(abs(bs(q, nu=degree)-expected)/(1+abs(expected)))
            derivative_max = max(derivative_max, error)
            assert error < 3e-8, (index, 'derivative', degree, error)
        np.testing.assert_allclose(integral, bs.integrate(u[0], u[-1]), rtol=3e-10, atol=3e-10)
        active = len(u)-int(periodic)
        fp = np.sum((w[:active, None]*(values[:active]-y[:active]))**2)
        residual_max = max(residual_max, abs(fp-target)/target)
        assert abs(fp-target) < .0011*target
        ci = c[:wrap.shape[1]]
        penalty = 1/history[-1][0]**2
        gradient = a.T@(a@ci-rhs)+penalty*d.T@(d@ci)
        denom = np.linalg.norm(a)*(np.linalg.norm(a@ci)+np.linalg.norm(rhs))
        error = np.linalg.norm(gradient)/denom
        stationarity_max = max(stationarity_max, error)
        assert error < 2e-11
        if not periodic:
            if scalar:
                ref_public = make_splrep(u, y[:, 0], w=w, k=k, s=target, t=t)
            else:
                ref_public, _ = make_splprep(y.T, u=u, w=w, k=k, s=target, t=t)
            error = np.max(abs(c-ref_public.c.reshape(c.shape))/(1+abs(c)))
            public_max = max(public_max, error)
            assert error < 3e-11, (index, 'public SciPy reference', error)
            public += 1
        steps += len(history)
    assert next(lines, None) is None
    print(f'PASS: {len(cases)} supplied-knot fits; {public} public SciPy checks; '
          f'{steps} independent reference penalty evaluations')
    print(f'Max coefficient error {coefficient_max:.6g}; derivative error {derivative_max:.6g}; '
          f'public SciPy error {public_max:.6g}; relative residual error {residual_max:.6g}; '
          f'stationarity {stationarity_max:.6g}')


if __name__ == '__main__':
    main()
