#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check periodic tuple least squares against SciPy and an independent cyclic solve."""
import argparse
from pathlib import Path
import subprocess
import numpy as np
import scipy
from scipy.interpolate import BSpline, splrep, splprep


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261115)
    cases = []
    for k in range(1, 6):
        for j in range(24):
            scalar = j % 2 == 0
            ndim = 1 if scalar else 2+j % 4
            m = 61+j
            u = np.r_[0., np.sort(rng.uniform(.001, .999, m-2)), 1.]
            x = np.column_stack([np.cos(2*np.pi*u+a)+.2*np.sin(6*np.pi*u-a)+
                                 .04*rng.normal(size=m) for a in range(ndim)])
            x[-1] = 1234.
            w = np.exp(rng.uniform(-.8, .8, m))
            provided = scalar or j % 3 != 0
            if provided:
                p = u
            else:
                closed = np.vstack((x[:-1], x[0]))
                p = np.r_[0., np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
                p /= p[-1]
            ni = (2, 3, 7, 15)[j % 4]
            # Most cases are comfortably conditioned; some deliberately lack observation rank.
            core = np.arange(ni)/ni
            if j >= 20:
                core[1:3] = [.00001, .00002][:min(2, ni-1)]
            t = np.array([core[i % ni]+i//ni for i in range(-k, ni+k+1)])
            interior = core[1:]
            wrap = np.eye(ni)[np.arange(ni+k) % ni]
            a = BSpline(t, wrap, k)(p[:-1])*w[:-1, None]
            rhs = x[:-1]*w[:-1, None]
            c = np.linalg.lstsq(a, rhs, rcond=None)[0]
            independent = BSpline(t, wrap@c, k)
            rank = np.linalg.matrix_rank(a)
            reference = None
            if rank == ni:
                if scalar:
                    rt, rc, rk = splrep(p, x[:, 0], w=w, k=k, task=-1,
                                        t=interior, per=True)
                    reference = BSpline(rt, rc[:len(rt)-k-1, None], rk)
                else:
                    template = np.r_[np.zeros(k+1), interior, np.ones(k+1)]
                    tck, _ = splprep(x.T.copy(), u=p, w=w, k=k, task=-1,
                                     t=template, per=True)
                    reference = BSpline(tck[0], np.array(tck[1]).T, k)
            q = np.linspace(-.05, 1.05, 27)
            cases.append((k, scalar, x, w, u if provided else None, p, interior,
                          q, independent, reference, a, rhs, rank))
    lines = []
    for k, scalar, x, w, u, p, interior, q, *_ in cases:
        lines.append(f'{len(x)} {x.shape[1]} {k} {int(scalar)} {int(u is not None)} '
                     f'{len(q)} {len(interior)}')
        arrays = [*x, w]
        if u is not None:
            arrays.append(u)
        arrays.extend((interior, q))
        lines.extend(' '.join(f'{v:.17g}' for v in row) for row in arrays)
    proc = subprocess.run([str(args.executable.resolve())], input='\n'.join(lines)+'\n',
                          capture_output=True, text=True, check=True)
    rows = iter(proc.stdout.splitlines())
    coefficient_max = value_max = stationarity_max = seam_max = 0.
    scipy_comparisons = deficient = 0
    for index, (k, scalar, x, w, u, p, interior, q, independent, ref, a, rhs, rank) in enumerate(cases):
        assert int(next(rows)) == 0, (index, 'status')
        nt = int(next(rows))
        actual_p = np.fromstring(next(rows), sep=' ')
        t = np.fromstring(next(rows), sep=' ')
        ndim = x.shape[1]
        c = np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
        values = np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
        queries = [np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
                   for _ in range(k+1)]
        integral = np.fromstring(next(rows), sep=' ')
        np.testing.assert_allclose(actual_p, p, rtol=2e-14, atol=2e-15)
        np.testing.assert_allclose(t, independent.t, rtol=2e-14, atol=2e-15)
        assert nt == len(t) and c.shape == (nt-k-1, ndim)
        ni = len(interior)+1
        deficient += rank < ni
        error = np.max(abs(c-independent.c)/(1+abs(independent.c)))
        coefficient_max = max(coefficient_max, error)
        assert error < 2e-7, (index, 'independent coefficients', error)
        residual = a@c[:ni]-rhs
        gradient = a.T@residual
        scale = np.linalg.norm(a)*(np.linalg.norm(rhs)+np.linalg.norm(a)*np.linalg.norm(c[:ni]))
        error = np.linalg.norm(gradient)/max(scale, 1e-300)
        stationarity_max = max(stationarity_max, error)
        assert error < 2e-12, (index, 'least-squares normal equations', error)
        bs = BSpline(t, c, k)
        for d in range(k):
            seam = bs([p[0], p[-1]], nu=d)
            scale = max(1., np.max(abs(c))*(k/np.min(np.diff(t)))**d)
            error = np.max(abs(seam[0]-seam[1]))/scale
            seam_max = max(seam_max, error)
            assert error < 2e-12, (index, 'unwrapped seam', error)
        pairs = [(values, independent(p)), (integral, independent.integrate(p[0], p[-1]))]
        pairs.extend((v, independent(q, nu=d)) for d, v in enumerate(queries))
        for actual, expected in pairs:
            error = np.max(abs(actual-expected)/(1+abs(expected)))
            value_max = max(value_max, error)
            assert error < 2e-7, (index, 'independent evaluation', error)
        if ref is not None:
            pairs = [(c, ref.c), (values, ref(p)), (integral, ref.integrate(p[0], p[-1]))]
            pairs.extend((v, ref(q, nu=d)) for d, v in enumerate(queries))
            for actual, expected in pairs:
                np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=2e-7,
                                           err_msg=f'case {index}, SciPy')
                scipy_comparisons += actual.size
    assert next(rows, None) is None
    print(f'PASS: {len(cases)} periodic LSQ fits; {deficient} rank-deficient; '
          f'{scipy_comparisons} SciPy {scipy.__version__} comparisons')
    print(f'Max normalized coefficient error {coefficient_max:.6g}; '
          f'value/derivative/integral error {value_max:.6g}; '
          f'stationarity {stationarity_max:.6g}; seam {seam_max:.6g}')


if __name__ == '__main__':
    main()
