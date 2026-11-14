#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Validate periodic scalar fits against SciPy or the unchanged pinned FITPACK core."""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.interpolate import splrep, splev, splint, BSpline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--upstream-executable', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261114)
    cases = []
    for k in range(1, 6):
        for j in range(36):
            n = k+1 if j == 0 else 24 + 2*j
            x = np.cumsum(rng.uniform(.5, 1.5, n))
            x = (x-x[0])/(x[-1]-x[0])
            y = np.cos(2*np.pi*x)+np.sin(4*np.pi*x)+.08*rng.normal(size=n)
            w = np.exp(rng.uniform(-1.2, 1.2, n))
            # The closing ordinate deliberately does not equal the first.
            y[-1] = 100
            target = (0, .01, .1, 2, .00001, .02)[j % 6]*sum((w[:-1]*y[:-1])**2)
            q = np.linspace(0, 1, 47)
            cases.append((k, x, y, w, target, q))
    # Reproduce the sine/cosine interpolation geometry of upstream test_periodic_fit.
    x = np.linspace(0, 2*np.pi, 200)
    cases.append((3, x, np.cos(x)+np.sin(2*x), np.ones(200), 0,
                  np.r_[x[0], (x[:-2]+x[1:-1])/2, x[-1]]))
    lines = []
    for k, x, y, w, target, q in cases:
        lines.append(f'{len(x)} {k} {len(q)}')
        for values in ([target, x[0], x[-1]], x, y, w, q):
            lines.append(' '.join(f'{v:.17g}' for v in values))
    payload = '\n'.join(lines)+'\n'

    def run(exe):
        proc = subprocess.run([str(exe.resolve())], input=payload,
                              text=True, capture_output=True, check=True)
        return iter(proc.stdout.splitlines())

    rows = run(args.executable)
    upstream = run(args.upstream_executable) if args.upstream_executable else None
    maximum = residual_max = stationarity_max = 0.0
    count = adaptive = differing_knots = differing_results = 0
    for index, (k, x, y, w, target, q) in enumerate(cases):
        assert int(next(rows)) == 0, ('native failure', index, k)
        nt = int(next(rows))
        knots = np.fromstring(next(rows), sep=' ')
        coeff = np.fromstring(next(rows), sep=' ')
        actual = np.fromstring(next(rows), sep=' ')
        values = [np.fromstring(next(rows), sep=' ') for _ in range(min(3, k)+1)]
        integral = float(next(rows))
        if upstream is None:
            ref = splrep(x, y, w=w, k=k, s=target, per=True)
            ref_actual = splev(x, ref)
            ref_values = [splev(q, ref, der=d) for d in range(min(3, k)+1)]
            ref_integral = splint(x[0], x[-1], ref)
        else:
            assert int(next(upstream)) == 0, ('upstream failure', index)
            assert int(next(upstream)) == nt, ('upstream knot count', index)
            ref = (np.fromstring(next(upstream), sep=' '),
                   np.fromstring(next(upstream), sep=' '), k)
            ref_actual = np.fromstring(next(upstream), sep=' ')
            ref_values = [np.fromstring(next(upstream), sep=' ')
                          for _ in range(min(3, k)+1)]
            ref_integral = float(next(upstream))
        assert nt == len(knots) == len(coeff), (index, 'tuple shape')
        same_knots = (nt == len(ref[0]) and
                      np.allclose(knots, ref[0], rtol=3e-13, atol=3e-13))
        if not same_knots:
            differing_knots += 1
        if upstream is not None:
            assert same_knots, (index, 'upstream knots')
        fp = sum((w[:-1]*(actual[:-1]-y[:-1]))**2)
        if target == 0:
            assert fp < 1e-22, (index, 'interpolation', fp)
        else:
            assert fp <= 1.0011*target, (index, 'residual', fp, target)
            if nt > 2*(k+1):
                adaptive += 1
                residual_max = max(residual_max, abs(fp-target)/target)
        # Independently verify seam derivatives and native representation evaluation.
        bs = BSpline(knots, coeff[:nt-k-1], k)
        if target > 0 and nt > 2*(k+1):
            nc = nt-k-1
            independent = nt-2*k-1
            wrap = np.eye(independent)[np.arange(nc) % independent]
            basis = BSpline(knots, wrap, k)
            design = basis(x[:-1])*w[:-1, None]
            rhs = w[:-1]*y[:-1]
            interior = knots[k+1:nt-k-1]
            jumps = (basis(np.nextafter(interior, np.inf), nu=k)-
                     basis(np.nextafter(interior, -np.inf), nu=k))
            jumps /= np.max(abs(jumps))
            c = coeff[:independent]
            gradient = design.T @ (design @ c-rhs)
            penalty = jumps.T @ (jumps @ c)
            lam = max(0., -np.dot(gradient, penalty)/np.dot(penalty, penalty))
            denom = (np.linalg.norm(design)*np.linalg.norm(design@c-rhs)+
                     lam*np.linalg.norm(jumps)**2*np.linalg.norm(c))
            error = np.linalg.norm(gradient+lam*penalty)/max(denom, 1e-300)
            assert error < 2e-8, (index, 'normal equations', error)
            stationarity_max = max(stationarity_max, error)
        for d in range(k):
            seam = bs([x[0], x[-1]], nu=d)
            assert abs(seam[0]-seam[1]) < 2e-8*(1+np.max(abs(seam))), (index, d, seam)
        for d, value in enumerate(values):
            np.testing.assert_allclose(value, bs(q, nu=d), rtol=3e-10, atol=3e-8)
        ref_fp = sum((w[:-1]*(ref_actual[:-1]-y[:-1]))**2)
        assert ref_fp <= 1.0011*target+1e-22, (index, 'reference target')
        if same_knots and target > 0:
            # For this common penalty family, monotonic residual eigencomponents
            # bound squared prediction distance by the residual difference.
            distance = sum((w[:-1]*(actual[:-1]-ref_actual[:-1]))**2)
            roundoff = 1e-10*max(1., sum((w[:-1]*y[:-1])**2))
            assert distance <= abs(fp-ref_fp)+roundoff, (index, 'prediction path')
        comparisons = [(actual, ref_actual),
                       (np.array([integral]), np.array([ref_integral]))]
        if same_knots:
            comparisons.append((coeff, ref[1]))
        comparisons.extend(zip(values, ref_values))
        for value, expected in comparisons:
            err = np.max(abs(value-expected)/(1+abs(expected)))
            if upstream is not None:
                assert err < 3e-8, (index, k, 'comparison', err)
            elif err >= 3e-8:
                differing_results += 1
            maximum = max(maximum, err)
            count += value.size
        if index == len(cases)-1:
            for d, value in enumerate(values):
                exact = (np.cos(q+d*np.pi/2)+2**d*np.sin(2*q+d*np.pi/2))
                # Upstream uses RTOL=.01 and ATOL=.0001 for this fixture.
                assert np.max(abs(value-exact)/(.0001+.01*abs(exact))) <= 1
    assert next(rows, None) is None
    if upstream is not None:
        assert next(upstream, None) is None
    print(f'PASS: {len(cases)} fits, {adaptive} adaptive; {count} comparisons')
    print(f'Reference: {"original perazz/FITPACK" if upstream else "SciPy "+scipy.__version__}')
    print(f'Maximum normalized difference {maximum:.6g}; '
          f'maximum adaptive residual relative error {residual_max:.6g}')
    print(f'Different knot sets: {differing_knots}; comparison arrays above 3e-8: '
          f'{differing_results}; maximum scaled stationarity error {stationarity_max:.6g}')


if __name__ == '__main__':
    main()
