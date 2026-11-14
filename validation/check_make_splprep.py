#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check native curve construction against SciPy and independent smoothing equations."""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.interpolate import BSpline, make_splprep, make_lsq_spline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261117)
    cases = []
    for k in range(1, 6):
        for j in range(32):
            m = 25+2*j
            ndim = 1+j % 5
            u = np.cumsum(rng.uniform(.5, 1.5, m))
            u = (u-u[0])/(u[-1]-u[0])
            x = np.column_stack([np.sin((4+a)*u)+.08*rng.normal(size=m)
                                 for a in range(ndim)])
            w = np.exp(rng.uniform(-.6, .6, m))
            mode = j % 4
            provided = j % 3 != 0
            bounds = (-.2, 1.3) if j % 7 == 0 and mode else (0., 1.)
            target = .01*np.sum((w[:, None]*x)**2) if mode else 0.
            if j % 7 == 1:
                target = 2*np.sum((w[:, None]*x)**2)
            options = dict(k=k, u=u if provided else None)
            if mode:
                options.update(w=w, s=target, ub=bounds[0], ue=bounds[1])
            ref, ru = make_splprep(x.T, **options)
            t = ref.t.copy() if mode == 2 else np.empty(0)
            if mode == 2:
                ref, ru = make_splprep(x.T, t=t, **options)
            cases.append((k, mode, x, w, u if provided else None, target, bounds, t,
                          np.linspace(*bounds, 37), ref, ru))
        # Independent polynomial represented on a nonminimal supplied knot vector.
        u = np.linspace(0, 1, 61)
        x = np.column_stack([np.sin(8*u)+.15*np.cos(79*u), np.cos(7*u)])
        w = .5+u
        t = make_splprep(x.T, u=u, w=w, k=k, s=.15)[0].t
        ref = make_lsq_spline(u, x, np.r_[np.zeros(k+1), np.ones(k+1)], k=k, w=w)
        for knot in t[k+1:-k-1]:
            ref = ref.insert_knot(knot)
        # Standardize references to SciPy's public curve orientation (channels,sites).
        ref = BSpline(ref.t, ref.c.T, k, axis=1)
        cases.append((k, 2, x, w, u, 100., (0., 1.), t, np.linspace(0, 1, 37), ref, u))
    lines = []
    for k, mode, x, w, u, target, bounds, t, q, _, _ in cases:
        lines.append(f'{len(x)} {x.shape[1]} {k} {int(u is not None)} {len(q)} {mode} {len(t)}')
        arrays = [[target, *bounds], *x, w]
        if u is not None:
            arrays.append(u)
        if len(t):
            arrays.append(t)
        arrays.append(q)
        lines.extend(' '.join(f'{v:.17g}' for v in row) for row in arrays)
    proc = subprocess.run([str(args.executable.resolve())], input='\n'.join(lines)+'\n',
                          capture_output=True, text=True, check=True)
    rows = iter(proc.stdout.splitlines())
    stat_max = residual_max = repr_max = diagnostic_max = fixed_max = interp_max = 0.
    adaptive = different_knots = count = 0
    for index, (k, mode, x, w, u, target, bounds, t, q, ref, ru) in enumerate(cases):
        code = int(next(rows))
        assert code == 0, (index, k, mode, code)
        nt = int(next(rows))
        au = np.fromstring(next(rows), sep=' ')
        knots = np.fromstring(next(rows), sep=' ')
        ndim = x.shape[1]
        coeff = np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
        values = np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
        queries = [np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
                   for _ in range(min(3, k)+1)]
        integral = np.fromstring(next(rows), sep=' ')
        assert len(knots) == nt and coeff.shape == (nt-k-1, ndim)
        np.testing.assert_allclose(au, ru, rtol=2e-14, atol=2e-15)
        if u is not None:
            np.testing.assert_array_equal(au, u)
        if len(t):
            np.testing.assert_array_equal(knots, t)
        same = len(knots) == len(ref.t) and np.allclose(knots, ref.t, atol=3e-13, rtol=3e-13)
        if not same and mode in (1, 3):
            different_knots += 1
        if mode in (0, 2):
            assert same, (index, 'interpolation/supplied knots')
        weights = w if mode else np.ones_like(w)
        fp = np.sum((weights[:, None]*(values-x))**2)
        assert fp <= 1.0011*target+1e-22, (index, 'residual', fp, target)
        if mode:
            basis = BSpline(knots, np.eye(len(coeff)), k)
            design = basis(au)*w[:, None]
            gradient = design.T@(design@coeff-w[:, None]*x)
            energy = max(1., np.linalg.norm(design.T@(w[:, None]*x)))
            interior = knots[k+1:-k-1]
            if len(interior):
                jumps = (basis(np.nextafter(interior, np.inf), nu=k)-
                         basis(np.nextafter(interior, -np.inf), nu=k))
                jumps /= np.max(abs(jumps))
                if fp < .9989*target:
                    assert np.linalg.norm(jumps@coeff) < 1e-9*max(1., np.linalg.norm(coeff))
                    continue_check = False
                else:
                    continue_check = True
                    adaptive += 1
                    residual_max = max(residual_max, abs(fp-target)/target)
                penalty = jumps.T@(jumps@coeff)
                lam = max(0., -np.sum(gradient*penalty)/max(np.sum(penalty**2), 1e-300))
                error = np.linalg.norm(gradient+lam*penalty)/max(
                    energy, np.linalg.norm(gradient), lam*np.linalg.norm(penalty))
                if continue_check:
                    assert error < 2e-9, (index, 'shared-multiplier normal equations', error)
                    stat_max = max(stat_max, error)
            else:
                error = np.linalg.norm(gradient)/energy
                assert error < 2e-9, (index, 'polynomial normal equations', error)
                stat_max = max(stat_max, error)
        native = BSpline(knots, coeff, k)
        independent = [(values, native(au)), (integral, native.integrate(*bounds))]
        independent.extend((v, native(q, nu=d)) for d, v in enumerate(queries))
        for actual, expected in independent:
            error = np.max(abs(actual-expected)/(1+abs(expected)))
            assert error < 3e-8, (index, 'representation', error)
            repr_max = max(repr_max, error)
        comparisons = [(values, ref(au).T), (integral, ref.integrate(*bounds))]
        comparisons.extend((v, ref(q, nu=d).T) for d, v in enumerate(queries))
        if mode in (0, 2):
            comparisons.append((coeff, ref.c))
        for actual, expected in comparisons:
            error = np.max(abs(actual-expected)/(1+abs(expected)))
            diagnostic_max = max(diagnostic_max, error)
            if mode == 0:
                assert error < 3e-8, (index, 'modern interpolation', error)
                interp_max = max(interp_max, error)
                count += actual.size
            if mode == 2:
                assert error < 3e-8, (index, 'modern supplied-knot fit', error)
                fixed_max = max(fixed_max, error)
                count += actual.size
    assert next(rows, None) is None
    print(f'PASS: {len(cases)} fits; {adaptive} active constraints; {count} strict SciPy '
          f'{scipy.__version__} interpolation/supplied-knot comparisons')
    print(f'Maximum interpolation error {interp_max:.6g}; residual relative error '
          f'{residual_max:.6g}; scaled stationarity {stat_max:.6g}; representation {repr_max:.6g}')
    print(f'Modern diagnostics: {different_knots} automatic knot differences; '
          f'fixed-knot maximum {fixed_max:.6g}; overall maximum {diagnostic_max:.6g}')


if __name__ == '__main__':
    main()
