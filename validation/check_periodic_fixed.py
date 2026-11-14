#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check fixed periodic fits against independent cyclic constrained least squares."""
import argparse
from pathlib import Path
import subprocess
import numpy as np
import scipy
from scipy.interpolate import BSpline
from scipy.optimize import brentq


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261114)
    cases = []
    for k in range(1, 6):
        for ni in (1, 2, 3, 4, 7, 13, 23):
            for scalar in (True, False):
                m = 81
                ndim = 1 if scalar else (12 if ni == 7 else 3)
                u = np.r_[0., np.sort(rng.uniform(.001, .999, m-2)), 1.]
                # Clustered knots exercise several knots within a sample interval.
                core = np.r_[0., np.sort(rng.uniform(.01, .99, ni-1))]
                t = np.array([core[i % ni]+i//ni for i in range(-k, ni+k+1)])
                wrap = np.eye(ni)[np.arange(ni+k) % ni]
                basis = BSpline(t, wrap, k)
                x = np.column_stack([np.cos(2*np.pi*u+j)+.13*np.sin(6*np.pi*u-j)
                                     + .03*rng.normal(size=m) for j in range(ndim)])
                x[-1] = x[0]
                w = np.exp(rng.uniform(-.6, .6, m))
                a = basis(u[:-1])*w[:-1, None]
                rhs = x[:-1]*w[:-1, None]
                lsq = np.linalg.lstsq(a, rhs, rcond=None)[0]
                low = np.sum((a@lsq-rhs)**2)
                mean = np.sum(w[:-1, None]**2*x[:-1], axis=0)/np.sum(w[:-1]**2)
                high = np.sum((w[:-1, None]*(x[:-1]-mean))**2)
                q = np.linspace(-.1, 1.1, 27)
                targets = [2*high] if ni == 1 else [low+.45*(high-low), 2*high]
                for target in targets:
                    cases.append((k, scalar, x, w, u, target, q, t, a, rhs, mean,
                                  int(np.linalg.matrix_rank(a) < ni)))
    lines = []
    for k, scalar, x, w, u, target, q, t, *_ in cases:
        lines.append(f'{len(u)} {x.shape[1]} {k} {int(scalar)} 1 {len(q)} 0 {len(t)}')
        lines.extend(' '.join(f'{v:.17g}' for v in row)
                     for row in [[target], *x, w, u, q, t])
    proc = subprocess.run([str(args.executable.resolve())], input='\n'.join(lines)+'\n',
                          capture_output=True, text=True, check=True)
    rows = iter(proc.stdout.splitlines())
    coefficient_max = residual_max = seam_max = representation_max = 0.
    active = rank_deficient = 0
    for index, (k, scalar, x, w, u, target, q, t, a, rhs, mean, deficient) in enumerate(cases):
        code = int(next(rows))
        assert code == 0, (index, k, len(t)-2*k-1, scalar, code)
        rank_deficient += deficient
        assert int(next(rows)) == len(t)
        p = np.fromstring(next(rows), sep=' ')
        actual_t = np.fromstring(next(rows), sep=' ')
        np.testing.assert_array_equal(actual_t, t)
        np.testing.assert_array_equal(p, u)
        ndim = x.shape[1]
        c = np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
        values = np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
        queries = [np.column_stack([np.fromstring(next(rows), sep=' ') for _ in range(ndim)])
                   for _ in range(min(3, k)+1)]
        integral = np.fromstring(next(rows), sep=' ')
        bs = BSpline(t, c, k)
        for d in range(k):
            seam = bs([0., 1.], nu=d)
            scale = max(1., np.max(abs(bs.derivative(d).c)),
                        np.max(abs(c))*(k/np.min(np.diff(t)))**d)
            error = np.max(abs(seam[0]-seam[1]))/scale
            seam_max = max(seam_max, error)
            assert error < 2e-12, (index, 'seam', d, error)
        pairs = [(values, bs(p)), (integral, bs.integrate(0., 1.))]
        pairs.extend((v, bs(q, nu=d)) for d, v in enumerate(queries))
        for actual, ref in pairs:
            error = np.max(abs(actual-ref)/(1+abs(ref)))
            representation_max = max(representation_max, error)
            assert error < 2e-8, (index, 'representation', error)
        fp = np.sum((w[:-1, None]*(values[:-1]-x[:-1]))**2)
        assert fp <= 1.0011*target, (index, 'residual', fp, target)
        ni = len(t)-2*k-1
        high = np.sum((w[:-1, None]*(x[:-1]-mean))**2)
        if target >= high:
            reference = np.tile(mean, (ni, 1))
        else:
            active += 1
            residual_max = max(residual_max, abs(fp-target)/target)
            wrap = np.eye(ni)[np.arange(ni+k) % ni]
            basis = BSpline(t, wrap, k)
            interior = t[k+1:-k-1]
            jumps = basis(np.nextafter(interior, np.inf), nu=k)-basis(
                np.nextafter(interior, -np.inf), nu=k)
            jumps /= np.max(abs(jumps))
            augmented_rhs = np.vstack((rhs, np.zeros((len(jumps), ndim))))
            def solve(lam):
                return np.linalg.lstsq(np.vstack((a, np.sqrt(lam)*jumps)),
                                       augmented_rhs, rcond=None)[0]
            def residual(lam):
                return np.sum((a@solve(lam)-rhs)**2)-fp
            upper = 1.
            while residual(upper) < 0:
                upper *= 10
                assert upper < 1e20, (index, 'multiplier bracket')
            lam = brentq(residual, 0., upper, xtol=1e-28, rtol=1e-14)
            reference = solve(lam)
        error = np.max(abs(c[:ni]-reference)/(1+abs(reference)))
        coefficient_max = max(coefficient_max, error)
        assert error < 2e-7, (index, 'independent constrained coefficients', error)
    assert next(rows, None) is None
    print(f'PASS: {len(cases)} fixed periodic cases; {active} active constraints; '
          f'{rank_deficient} rank-deficient fits; '
          f'SciPy {scipy.__version__} independent cyclic basis and constrained solve')
    print(f'Max normalized coefficient error {coefficient_max:.6g}; '
          f'representation {representation_max:.6g}; seam {seam_max:.6g}; '
          f'relative target error {residual_max:.6g}')


if __name__ == '__main__':
    main()
