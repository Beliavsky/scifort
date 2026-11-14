#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check repeated-knot extension using SciPy piecewise polynomials and constrained LSQ.

SciPy smoothing constructors reject these knots; this is not smoothing API
parity. PPoly leading coefficients give one kth-derivative jump per distinct
interior break. Derivatives here are piecewise, not distributional derivatives."""
import argparse
from pathlib import Path
import subprocess
import math
import numpy as np
import scipy
from scipy.interpolate import BSpline, PPoly, make_lsq_spline
from scipy.optimize import brentq
from scipy.linalg import null_space


def effective_condition(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    retained = singular[singular > max(matrix.shape)*np.finfo(float).eps*singular[0]]
    return singular[0]/retained[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261114)
    cases = []
    for k in range(1, 6):
        for form in range(3):
            for periodic in (False, True):
                for scalar in (False, True):
                    m = 61
                    ndim = 1 if scalar else 3
                    u = np.linspace(0., 1., m)
                    multiplicity = 2+form % k
                    if form == 2:
                        multiplicity = k+1
                    interior = np.repeat([.25, .5, .75], multiplicity)
                    if form == 1:
                        interior = np.repeat([.25, .5, .75], [2, min(k+1, 3), k+1])
                    if form == 2:
                        interior = np.repeat([.001, .002, .5], multiplicity)
                    left = -np.arange(k, 0, -1)*.09
                    right = 1+np.arange(1, k+1)*.13
                    if form != 1:
                        left[:] = 0.
                    if form != 1:
                        right[:] = 1.
                    t = np.r_[left, 0., interior, 1., right]
                    nc = len(t)-k-1
                    basis = BSpline(t, np.eye(nc), k)
                    if periodic:
                        constraint = np.array([basis(0., nu=d)-basis(1., nu=d)
                                               for d in range(k)])
                        constraint /= np.max(abs(constraint), axis=1)[:, None]
                        z = null_space(constraint)
                    else:
                        z = np.eye(nc)
                    x = np.column_stack([np.cos(2*np.pi*u+j)+.1*np.sin(6*np.pi*u-j)
                                         + .02*rng.normal(size=m) for j in range(ndim)])
                    x[-1] = x[0]
                    w = np.exp(rng.uniform(-.5, .5, m))
                    active = m-1 if periodic else m
                    a = (basis(u[:active])*w[:active, None])@z
                    rhs = x[:active]*w[:active, None]
                    # Read polynomial coefficients, independently of native derivative sampling.
                    intervals = (t[:-1] >= 0.) & (t[1:] <= 1.) & (np.diff(t) > 0.)
                    leading = np.column_stack([
                        PPoly.from_spline((t, np.eye(nc)[:, j], k)).c[0, intervals]
                        for j in range(nc)])
                    jumps = np.diff(leading, axis=0)*math.factorial(k)
                    d = jumps@z
                    d /= np.max(abs(d))
                    zero_basis = z@null_space(d)
                    full_a = basis(u[:active])*w[:active, None]
                    zero_coeff = np.linalg.lstsq(full_a@zero_basis, rhs, rcond=None)[0]
                    limit = zero_basis@zero_coeff
                    high = np.sum((full_a@limit-rhs)**2)
                    low = np.sum((a@np.linalg.lstsq(a, rhs, rcond=None)[0]-rhs)**2)
                    breaks = np.unique(t[(t > 0.) & (t < 1.)])
                    q = np.r_[np.linspace(-.05, 1.05, 23), breaks,
                              np.nextafter(breaks, -np.inf), np.nextafter(breaks, np.inf)]
                    assert high > low+1e-10, (k, form, periodic, high, low)
                    # Check ordinary repeated-knot LSQ against the installed public SciPy API.
                    if not periodic and np.linalg.matrix_rank(a) == a.shape[1]:
                        lsq = make_lsq_spline(u, x, t, k=k, w=w)
                        reference = np.linalg.lstsq(full_a, rhs, rcond=None)[0]
                        np.testing.assert_allclose(lsq.c, reference, rtol=2e-8, atol=2e-8)
                    for target in [low+.4*(high-low), 2*high]:
                        cases.append((k, scalar, x, w, u, target, q, t, a, rhs,
                                      z, d, limit, periodic))
    lines = []
    for k, scalar, x, w, u, target, q, t, a, rhs, z, d, limit, periodic in cases:
        lines.append(f'{len(u)} {x.shape[1]} {k} {int(scalar)} 1 {len(q)} 0 {len(t)} {int(periodic)}')
        lines.extend(' '.join(f'{v:.17g}' for v in row)
                     for row in [[target], *x, w, u, q, t])
    proc = subprocess.run([str(args.executable.resolve())], input='\n'.join(lines)+'\n',
                          capture_output=True, text=True, check=True)
    rows = iter(proc.stdout.splitlines())
    coefficient_max = residual_max = seam_max = representation_max = 0.
    active = rank_deficient = 0
    null_max = 0.
    condition_max = stationarity_max = scaled_coefficient_max = 0.
    for index, (k, scalar, x, w, u, target, q, t, a, rhs, z, jumps, limit, periodic) in enumerate(cases):
        code = int(next(rows))
        assert code == 0, (index, k, periodic, scalar, code)
        if np.linalg.matrix_rank(a) < a.shape[1]:
            rank_deficient += 1
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
        for d in range(k if periodic else 0):
            seam = bs([0., 1.], nu=d)
            scale = max(1., 
                        np.max(abs(c))*(k/np.min(np.diff(np.unique(t))))**d)
            error = np.max(abs(seam[0]-seam[1]))/scale
            seam_max = max(seam_max, error)
            assert error < 2e-12, (index, 'seam', d, error)
        identity_basis = BSpline(t, np.eye(len(c)), k)
        pairs = [(values, bs(p), abs(identity_basis(p))@abs(c)),
                 (integral, bs.integrate(0., 1.), abs(identity_basis.integrate(0., 1.))@abs(c))]
        pairs.extend((v, bs(q, nu=d), abs(identity_basis(q, nu=d))@abs(c))
                     for d, v in enumerate(queries))
        for actual, ref, roundoff in pairs:
            # Derivative cancellation at clustered knots requires an absolute rounding bound.
            error = np.max(abs(actual-ref)/(1+abs(ref)+64*np.finfo(float).eps*roundoff/2e-8))
            representation_max = max(representation_max, error)
            assert error < 2e-8, (index, 'representation', error)
        nactive = len(u)-int(periodic)
        fp = np.sum((w[:nactive, None]*(values[:nactive]-x[:nactive]))**2)
        assert fp <= 1.0011*target, (index, 'residual', fp, target)
        high = np.sum(((BSpline(t, np.eye(len(c)), k)(u[:nactive])*
                            w[:nactive, None])@limit-rhs)**2)
        if target >= high:
            reference = limit
        else:
            active += 1
            residual_max = max(residual_max, abs(fp-target)/target)
            assert abs(fp-target) < .0011*target, (index, 'active target', fp, target, high)
            augmented_rhs = np.vstack((rhs, np.zeros((len(jumps), ndim))))
            def solve(mu):
                return np.linalg.lstsq(np.vstack((a, mu*jumps)),
                                       augmented_rhs, rcond=None)[0]
            def residual(mu):
                return np.sum((a@solve(mu)-rhs)**2)-fp
            upper = 1.
            while residual(upper) < 0:
                upper *= 10
                assert upper < 1e20, (index, 'multiplier bracket')
            mu = brentq(residual, 0., upper, xtol=1e-28, rtol=1e-14)
            reference = z@solve(mu)
        error = np.max(abs(c-reference)/(1+abs(reference)))
        coefficient_max = max(coefficient_max, error)
        if target < high:
            condition = effective_condition(np.vstack((a,mu*jumps)))
            native = z.T@c
            residual_vector = a@native-rhs
            gradient = a.T@residual_vector
            penalty_vector = jumps.T@(jumps@native)
            lam = mu**2
            denominator = (np.linalg.norm(a)*np.linalg.norm(residual_vector)+
                           lam*np.linalg.norm(jumps)**2*np.linalg.norm(native))
            stationarity = np.linalg.norm(gradient+lam*penalty_vector)/max(denominator, 1e-300)
            stationarity_max = max(stationarity_max, stationarity)
            assert stationarity < 2e-8, (index, 'stationarity', stationarity)
        else:
            full_a = identity_basis(u[:nactive])*w[:nactive, None]
            zero_basis = z@null_space(jumps)
            condition = effective_condition(full_a@zero_basis)
        joint = null_space(np.vstack((a, jumps)))
        null_error = np.linalg.norm(joint.T@(z.T@c))/max(1., np.linalg.norm(c))
        null_max = max(null_max, null_error)
        assert null_error < 2e-7, (index, 'minimum-norm null component', null_error)
        condition_max = max(condition_max, condition)
        tolerance = 2e-7+16*np.finfo(float).eps*condition
        scaled_coefficient_max = max(scaled_coefficient_max, error/tolerance)
        assert error < tolerance, (index, 'conditioned constrained coefficients', error, condition)
    assert next(rows, None) is None
    print(f'PASS: {len(cases)} repeated-knot cases; {active} active constraints; '
          f'{rank_deficient} rank-deficient cases; '
          f'SciPy {scipy.__version__} independent cyclic basis and constrained solve')
    print(f'Max normalized coefficient error {coefficient_max:.6g}; '
          f'representation {representation_max:.6g}; seam {seam_max:.6g}; '
          f'relative target error {residual_max:.6g}')
    print(f'Max minimum-norm null component {null_max:.6g}')
    print(f'Max retained-subspace condition {condition_max:.6g}; stationarity {stationarity_max:.6g}; '
          f'coefficient error / conditioned tolerance {scaled_coefficient_max:.6g}')


if __name__ == '__main__':
    main()
