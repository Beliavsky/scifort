#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check native closed parametric FITPACK fits and independent fitting contracts."""
import argparse
from pathlib import Path
import subprocess
import re

import numpy as np
import scipy
from scipy.interpolate import splprep, splev, BSpline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--upstream-executable', type=Path)
    parser.add_argument('--upstream-tests', type=Path,
                        help='Pinned fitpack_curve_tests.f90, for test_closed_fit data')
    args = parser.parse_args()
    rng = np.random.default_rng(20261115)
    cases = []
    for k in range(1, 6):
        for j in range(36):
            m = k+1 if j == 0 else 24+2*j
            ndim = 1+j % 5
            u = np.cumsum(rng.uniform(.5, 1.5, m))
            u = (u-u[0])/(u[-1]-u[0])
            x = np.column_stack([np.cos(2*np.pi*u+c)+.3*np.sin(4*np.pi*u-c)+
                                 .06*rng.normal(size=m) for c in range(ndim)])
            # The final row is ignored and replaced by the first before chord lengths.
            x[-1] = 1234
            w = np.exp(rng.uniform(-1.2, 1.2, m))
            supplied = j % 2 == 0
            target = (0, .01, .1, 2, .00001, .02)[j % 6]*sum(
                (w[:-1, None]*x[:-1])**2).sum()
            cases.append((k, x, u, w, supplied, target, np.linspace(0, 1, 47)))
    if args.upstream_tests:
        source = args.upstream_tests.read_text(encoding='utf-8')
        fixture = source.split('function test_closed_fit(', 1)[1]
        fixture = fixture.split('end function test_closed_fit', 1)[0]
        channels = []
        for channel in (1, 2):
            literals = re.search(r'x\('+str(channel)+r',1:18\)\s*=\s*\[([^]]+)\]',
                                 fixture).group(1)
            channels.append(np.array([float(v) for v in literals.replace('&', '').split(',')],
                                     dtype=np.float32).astype(float))
        x = np.column_stack(channels)
        x = np.vstack([x, x[0]])
        u = np.linspace(0, 1, len(x))
        for k, target in ((3, 900), (3, 10), (3, .1), (3, .5), (5, .5), (5, 0), (3, 0)):
            cases.append((k, x, u, np.ones(len(x)), False, target, np.linspace(0, 1, 47)))
    lines = []
    for k, x, u, w, supplied, target, q in cases:
        lines.append(f'{len(u)} {x.shape[1]} {k} {int(supplied)} {len(q)}')
        values = [[target, 0, 1], *x, w]
        if supplied:
            values.append(u)
        values.append(q)
        for row in values:
            lines.append(' '.join(f'{v:.17g}' for v in row))
    payload = '\n'.join(lines)+'\n'

    def run(exe):
        proc = subprocess.run([str(exe.resolve())], input=payload,
                              text=True, capture_output=True, check=True)
        return iter(proc.stdout.splitlines())

    def read(rows, ndim, m, nq, degree):
        status = int(next(rows))
        assert status == 0, ('fit failure', status)
        nt = int(next(rows))
        parameter = np.fromstring(next(rows), sep=' ')
        knots = np.fromstring(next(rows), sep=' ')
        coeff = np.array([np.fromstring(next(rows), sep=' ') for _ in range(ndim)]).T
        values = np.array([np.fromstring(next(rows), sep=' ') for _ in range(ndim)]).T
        queries = [np.array([np.fromstring(next(rows), sep=' ')
                             for _ in range(ndim)]).T for _ in range(min(3, degree)+1)]
        assert len(knots) == nt and coeff.shape == (nt-degree-1, ndim)
        assert parameter.shape == (m,) and values.shape == (m, ndim)
        assert all(v.shape == (nq, ndim) for v in queries)
        return parameter, knots, coeff, values, queries

    rows = run(args.executable)
    upstream = run(args.upstream_executable) if args.upstream_executable else None
    count = adaptive = different_knots = different_arrays = 0
    maximum = residual_max = stationarity_max = representation_max = 0.0
    for index, (k, x, u, w, supplied, target, q) in enumerate(cases):
        ndim = x.shape[1]
        parameter, knots, coeff, values, queries = read(rows, ndim, len(u), len(q), k)
        if upstream is not None:
            rp, rt, rc, rv, rq = read(upstream, ndim, len(u), len(q), k)
        else:
            tck, rp = splprep(x.T.copy(), u=u if supplied else None, w=w,
                             k=k, s=target, per=True)
            rt, rc = tck[0], np.array(tck[1]).T
            rv = np.array(splev(rp, tck)).T
            rq = [np.array(splev(q, tck, der=d)).T for d in range(min(3, k)+1)]
        np.testing.assert_allclose(parameter, rp, rtol=3e-13, atol=3e-13)
        if supplied:
            np.testing.assert_array_equal(parameter, u)
        else:
            closed = x.copy()
            closed[-1] = closed[0]
            chord = np.r_[0, np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
            np.testing.assert_allclose(parameter, chord/chord[-1], rtol=3e-13, atol=3e-13)
        same = len(knots) == len(rt) and np.allclose(knots, rt, rtol=3e-13, atol=3e-13)
        if not same:
            different_knots += 1
        if upstream is not None:
            assert same, (index, 'upstream knots')
        fp = np.sum((w[:-1, None]*(values[:-1]-x[:-1]))**2)
        ref_fp = np.sum((w[:-1, None]*(rv[:-1]-x[:-1]))**2)
        assert fp <= 1.0011*target+1e-22, (index, 'native target', fp, target)
        assert ref_fp <= 1.0011*target+1e-22, (index, 'reference target')
        if target > 0 and len(knots) > 2*(k+1):
            adaptive += 1
            residual_max = max(residual_max, abs(fp-target)/target)
        bs = BSpline(knots, coeff, k)
        np.testing.assert_allclose(values, bs(parameter), rtol=3e-12, atol=3e-12)
        for d in range(k):
            seam = bs([parameter[0], parameter[-1]], nu=d)
            # Bound cancellation using the derivative coefficient scale.
            scale = max(1., np.max(abs(bs.derivative(d).c)))
            assert np.max(abs(seam[0]-seam[1])) <= 2048*np.finfo(float).eps*scale, (
                index, 'seam', d, seam)
        for d, query in enumerate(queries):
            expected = bs(q, nu=d)
            error = np.max(abs(query-expected)/(1+abs(expected)))
            representation_max = max(representation_max, error)
            assert error < 3e-8, (index, 'representation', d, error)
        if target > 0 and len(knots) > 2*(k+1):
            nc = len(coeff)
            independent = len(knots)-2*k-1
            wrap = np.eye(independent)[np.arange(nc) % independent]
            basis = BSpline(knots, wrap, k)
            design = basis(parameter[:-1])*w[:-1, None]
            rhs = x[:-1]*w[:-1, None]
            interior = knots[k+1:-k-1]
            jumps = (basis(np.nextafter(interior, np.inf), nu=k)-
                     basis(np.nextafter(interior, -np.inf), nu=k))
            jumps /= np.max(abs(jumps))
            c = coeff[:independent]
            gradient = design.T @ (design@c-rhs)
            penalty = jumps.T @ (jumps@c)
            lam = max(0., -np.sum(gradient*penalty)/np.sum(penalty**2))
            denom = (np.linalg.norm(design)*np.linalg.norm(design@c-rhs)+
                     lam*np.linalg.norm(jumps)**2*np.linalg.norm(c))
            error = np.linalg.norm(gradient+lam*penalty)/max(denom, 1e-300)
            assert error < 2e-8, (index, 'shared normal equations', error)
            stationarity_max = max(stationarity_max, error)
        if same and target > 0:
            distance = np.sum((w[:-1, None]*(values[:-1]-rv[:-1]))**2)
            roundoff = 1e-10*max(1., np.sum((w[:-1, None]*x[:-1])**2))
            assert distance <= abs(fp-ref_fp)+roundoff, (index, 'prediction path')
        comparisons = [(values, rv), *zip(queries, rq)]
        if same:
            comparisons.append((coeff, rc))
        for value, expected in comparisons:
            error = np.max(abs(value-expected)/(1+abs(expected)))
            maximum = max(maximum, error)
            count += value.size
            if upstream is not None:
                assert error < 3e-8, (index, k, 'upstream difference', error)
            else:
                if error >= 3e-8:
                    different_arrays += 1
                assert error < 3e-8, (index, k, 'SciPy difference', error)
    assert next(rows, None) is None
    if upstream is not None:
        assert next(upstream, None) is None
    print(f'PASS: {len(cases)} fits, {adaptive} adaptive; {count} compared numbers')
    print(f'Reference: {"original perazz/FITPACK" if upstream else "SciPy "+scipy.__version__}')
    print(f'Maximum normalized difference {maximum:.6g}; '
          f'maximum adaptive residual relative error {residual_max:.6g}')
    print(f'Different knot sets: {different_knots}; diagnostic arrays above 3e-8: '
          f'{different_arrays}; maximum scaled stationarity error {stationarity_max:.6g}; '
          f'maximum independent representation error {representation_max:.6g}')


if __name__ == '__main__':
    main()
