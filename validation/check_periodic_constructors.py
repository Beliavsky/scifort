#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check periodic constructor objects, interpolation and cyclic smoothing equations."""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.interpolate import BSpline, make_interp_spline, splrep, splprep


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261118)
    cases = []
    for k in range(1, 6):
        for j in range(40):
            m = k+1 if j in (0, 4) else 25+2*j
            scalar = j % 8 < 4
            ndim = 1 if scalar else 1+j % 5
            u = np.cumsum(rng.uniform(.5, 1.5, m))
            u = (u-u[0])/(u[-1]-u[0])
            x = np.column_stack([np.cos(2*np.pi*u+a)+.3*np.sin(4*np.pi*u-a)+
                                 .05*rng.normal(size=m) for a in range(ndim)])
            x[-1] = x[0]
            w = np.exp(rng.uniform(-.8, .8, m))
            target = (0., .002, .02, 2.)[j % 4]*np.sum((w[:-1, None]*x[:-1])**2)
            provided = scalar or j % 3 != 0
            capacity = m+2*k if target and j % 7 == 0 else 0
            if not provided:
                ru = np.r_[0., np.cumsum(np.linalg.norm(np.diff(x, axis=0), axis=1))]
                ru /= ru[-1]
            else:
                ru = u
            if target == 0:
                ref = make_interp_spline(ru, x, k=k, bc_type='periodic')
            elif scalar:
                rt, rc, rk = splrep(ru, x[:, 0], w=w, k=k, s=target, per=True)
                ref = BSpline(rt, rc[:len(rt)-k-1, None], rk)
            else:
                tck, _ = splprep(x.T.copy(), u=ru, w=w, k=k, s=target, per=True)
                ref = BSpline(tck[0], np.array(tck[1]).T, k)
            q = np.linspace(-.1, 1.1, 37)
            cases.append((k, scalar, x, w, u if provided else None, target, capacity, q, ref, ru))
    lines = []
    for k, scalar, x, w, u, target, capacity, q, _, _ in cases:
        lines.append(f'{len(x)} {x.shape[1]} {k} {int(scalar)} {int(u is not None)} '
                     f'{len(q)} {capacity}')
        arrays = [[target], *x, w]
        if u is not None:
            arrays.append(u)
        arrays.append(q)
        lines.extend(' '.join(f'{v:.17g}' for v in row) for row in arrays)
    proc = subprocess.run([str(args.executable.resolve())], input='\n'.join(lines)+'\n',
                          capture_output=True, text=True, check=True)
    rows = iter(proc.stdout.splitlines())
    interp_max = repr_max = seam_max = residual_max = stationarity_max = diagnostic_max = 0.
    different_knots = count = adaptive = 0
    for index, (k, scalar, x, w, u, target, capacity, q, ref, ru) in enumerate(cases):
        code = int(next(rows))
        assert code == 0, (index, k, scalar, target, code)
        nt = int(next(rows))
        p = np.fromstring(next(rows), sep=' ')
        t = np.fromstring(next(rows), sep=' ')
        ndim = x.shape[1]
        c = np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
        values = np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
        queries = [np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
                   for _ in range(min(3, k)+1)]
        integral = np.fromstring(next(rows), sep=' ')
        assert nt == len(t) and c.shape == (nt-k-1, ndim)
        if capacity:
            assert nt <= capacity
        np.testing.assert_allclose(p, ru, rtol=2e-14, atol=2e-15)
        same = len(t) == len(ref.t) and np.allclose(t, ref.t, atol=3e-13, rtol=3e-13)
        if not same:
            different_knots += 1
        bs = BSpline(t, c, k)
        for d in range(k):
            seam = bs([p[0], p[-1]], nu=d)
            scale = max(1., np.max(abs(bs.derivative(d).c)))
            error = np.max(abs(seam[0]-seam[1]))/scale
            seam_max = max(seam_max, error)
            assert error < 2048*np.finfo(float).eps, (index, 'unwrapped seam', d, error)
        independent = [(values, bs(p)), (integral, bs.integrate(p[0], p[-1]))]
        for d, v in enumerate(queries):
            independent.append((v, bs(q, nu=d, extrapolate='periodic'
                                      if target == 0 and k > 1 else True)))
        for a, b in independent:
            error = np.max(abs(a-b)/(1+abs(b)))
            repr_max = max(repr_max, error)
            assert error < 3e-8, (index, 'representation', error)
        fp = np.sum((w[:-1, None]*(values[:-1]-x[:-1]))**2)
        assert fp <= 1.0011*target+1e-22, (index, 'residual', fp, target)
        if target and nt > 2*(k+1):
            adaptive += 1
            residual_max = max(residual_max, abs(fp-target)/target)
            nc = len(c)
            ni = nt-2*k-1
            wrap = np.eye(ni)[np.arange(nc) % ni]
            basis = BSpline(t, wrap, k)
            design = basis(p[:-1])*w[:-1, None]
            rhs = x[:-1]*w[:-1, None]
            interior = t[k+1:-k-1]
            jumps = basis(np.nextafter(interior, np.inf), nu=k)-basis(
                np.nextafter(interior, -np.inf), nu=k)
            jumps /= np.max(abs(jumps))
            gradient = design.T@(design@c[:ni]-rhs)
            penalty = jumps.T@(jumps@c[:ni])
            lam = max(0., -np.sum(gradient*penalty)/np.sum(penalty**2))
            denom = (np.linalg.norm(design)*np.linalg.norm(design@c[:ni]-rhs)+
                     lam*np.linalg.norm(jumps)**2*np.linalg.norm(c[:ni]))
            error = np.linalg.norm(gradient+lam*penalty)/max(denom, 1e-300)
            stationarity_max = max(stationarity_max, error)
            assert error < 2e-8, (index, 'cyclic shared-multiplier equations', error)
        comparisons = [(values, ref(p)), (integral, ref.integrate(p[0], p[-1]))]
        comparisons.extend((v, ref(q, nu=d)) for d, v in enumerate(queries))
        if target == 0:
            assert same, (index, 'interpolation knots')
            comparisons.append((c, ref.c))
        for a, b in comparisons:
            error = np.max(abs(a-b)/(1+abs(b)))
            diagnostic_max = max(diagnostic_max, error)
            if target == 0:
                interp_max = max(interp_max, error)
                count += a.size
                assert error < 3e-8, (index, 'SciPy interpolation', error)
    assert next(rows, None) is None
    print(f'PASS: {len(cases)} fits; {adaptive} active constraints; {count} strict '
          f'SciPy {scipy.__version__} interpolation comparisons')
    print(f'Max interpolation {interp_max:.6g}; representation {repr_max:.6g}; '
          f'seam {seam_max:.6g}; residual relative error {residual_max:.6g}; '
          f'stationarity {stationarity_max:.6g}')
    print(f'Legacy smoothing diagnostics: {different_knots} knot differences; '
          f'maximum normalized value/derivative/integral difference {diagnostic_max:.6g}')


if __name__ == '__main__':
    main()
