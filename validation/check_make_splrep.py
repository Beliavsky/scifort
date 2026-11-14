#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare the modern scalar constructor with SciPy and its smoothing equations."""
import argparse
from pathlib import Path
import subprocess
import numpy as np
import scipy
from scipy.interpolate import make_splrep, make_lsq_spline, splrep, BSpline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261116)
    cases = []
    for k in range(1, 6):
        for j in range(32):
            m = 25+2*j
            x = np.cumsum(rng.uniform(.5, 1.5, m))
            x = (x-x[0])/(x[-1]-x[0])
            y = np.sin((4+j % 4)*x)+.08*rng.normal(size=m)
            w = np.exp(rng.uniform(-1, 1, m))
            mode = j % 4
            bounds = (0., 1.)
            if j % 7 == 0 and mode:
                bounds = (-.2, 1.3)
            target = .01*np.sum((w*y)**2)
            if j % 7 == 1:
                target = 2*np.sum((w*y)**2)
            t = np.empty(0)
            if mode == 0:
                target = 0
                reference = make_splrep(x, y, k=k)
            else:
                reference = make_splrep(x, y, w=w, k=k, s=target,
                                        xb=bounds[0], xe=bounds[1])
                if mode == 2:
                    t = reference.t.copy()
                    reference = make_splrep(x, y, w=w, k=k, s=target, t=t,
                                            xb=bounds[0], xe=bounds[1])
            q = np.linspace(*bounds, 47)
            cases.append((k, mode, x, y, w, target, bounds, t, q, reference))
        # A nonminimal supplied knot vector with an inactive residual constraint:
        # verify the polynomial representation independently by exact knot insertion.
        x = np.linspace(0, 1, 61)
        y = np.sin(8*x)+.15*np.cos(79*x)
        w = .5+x
        t = make_splrep(x, y, w=w, k=k, s=.15).t
        reference = make_lsq_spline(x, y, np.r_[np.zeros(k+1), np.ones(k+1)], k=k, w=w)
        for knot in t[k+1:-k-1]:
            reference = reference.insert_knot(knot)
        cases.append((k, 2, x, y, w, 100., (0., 1.), t, np.linspace(0, 1, 47), reference))
    lines = []
    for k, mode, x, y, w, target, bounds, t, q, _ in cases:
        lines.append(f'{len(x)} {k} {mode} {len(t)} {len(x)+k+1} {len(q)}')
        arrays = [[target, *bounds], x, y, w]
        if len(t):
            arrays.append(t)
        arrays.append(q)
        for row in arrays:
            lines.append(' '.join(f'{v:.17g}' for v in row))
    proc = subprocess.run([str(args.executable.resolve())], input='\n'.join(lines)+'\n',
                          capture_output=True, text=True, check=True)
    rows = iter(proc.stdout.splitlines())
    maximum = stationarity_max = residual_max = diagnostic_max = 0.
    count = adaptive = different_knots = 0
    for index, (k, mode, x, y, w, target, bounds, t, q, reference) in enumerate(cases):
        code = int(next(rows))
        assert code == 0, (index, k, mode, code)
        nt = int(next(rows))
        knots = np.fromstring(next(rows), sep=' ')
        coeff = np.fromstring(next(rows), sep=' ')
        values = np.fromstring(next(rows), sep=' ')
        queries = [np.fromstring(next(rows), sep=' ') for _ in range(min(3, k)+1)]
        integral = float(next(rows))
        assert nt == len(knots)
        assert len(coeff) == nt-k-1
        same_knots = (nt == len(reference.t) and
                      np.allclose(knots, reference.t, rtol=3e-13, atol=3e-13))
        if not same_knots:
            different_knots += 1
        if mode in (0, 2):
            assert same_knots, (index, 'specified/interpolation knots')
        if len(t):
            np.testing.assert_array_equal(knots, t)
        weights = w if mode else np.ones_like(w)
        fp = np.sum((weights*(values-y))**2)
        assert fp <= 1.0011*target+1e-22, (index, 'residual', fp, target)
        if mode and nt > 2*(k+1):
            basis = BSpline(knots, np.eye(len(coeff)), k)
            design = basis(x)*w[:, None]
            interior = knots[k+1:-k-1]
            jumps = (basis(np.nextafter(interior, np.inf), nu=k)-
                     basis(np.nextafter(interior, -np.inf), nu=k))
            jumps /= np.max(abs(jumps))
            if fp < .9989*target:
                assert np.linalg.norm(jumps@coeff) < 1e-10*max(1., np.linalg.norm(coeff))
            else:
                adaptive += 1
                residual_max = max(residual_max, abs(fp-target)/target)
            gradient = design.T@(design@coeff-w*y)
            penalty = jumps.T@(jumps@coeff)
            lam = max(0., -np.dot(gradient, penalty)/np.dot(penalty, penalty))
            denom = (np.linalg.norm(design)*np.linalg.norm(design@coeff-w*y)+
                     lam*np.linalg.norm(jumps)**2*np.linalg.norm(coeff))
            error = np.linalg.norm(gradient+lam*penalty)/max(denom, 1e-300)
            if fp >= .9989*target:
                assert error < 2e-8, (index, 'normal equations', error)
                stationarity_max = max(stationarity_max, error)
        bs = BSpline(knots, coeff, k)
        for d, query in enumerate(queries):
            np.testing.assert_allclose(query, bs(q, nu=d), rtol=3e-10, atol=3e-8)
        diagnostic = [(values, reference(x)),
                      (np.array([integral]), np.array([reference.integrate(*bounds)]))]
        diagnostic.extend((v, reference(q, nu=d)) for d, v in enumerate(queries))
        for actual, expected in diagnostic:
            diagnostic_max = max(diagnostic_max,
                                 np.max(abs(actual-expected)/(1+abs(expected))))
        if mode in (1, 3):
            legacy = splrep(x, y, w=w, k=k, s=target, xb=bounds[0], xe=bounds[1])
            reference = BSpline(*legacy)
            np.testing.assert_allclose(knots, reference.t, rtol=3e-13, atol=3e-13)
        comparisons = [(coeff, reference.c[:nt-k-1]), (values, reference(x)),
                       (np.array([integral]), np.array([reference.integrate(*bounds)]))]
        comparisons.extend((v, reference(q, nu=d)) for d, v in enumerate(queries))
        for actual, expected in comparisons:
            error = np.max(abs(actual-expected)/(1+abs(expected)))
            assert error < 3e-8, (index, k, mode, 'SciPy difference', error)
            maximum = max(maximum, error)
            count += actual.size
    assert next(rows, None) is None
    print(f'PASS: {len(cases)} fits; {adaptive} adaptive; {count} SciPy {scipy.__version__} comparisons')
    print(f'Maximum normalized difference {maximum:.6g}; residual relative error '
          f'{residual_max:.6g}; scaled stationarity error {stationarity_max:.6g}')
    print(f'Modern automatic knot differences: {different_knots}; maximum modern '
          f'value/derivative/integral diagnostic difference {diagnostic_max:.6g}')


if __name__ == '__main__':
    main()
