#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Compare explicit native continuation with exact cached upstream FITPACK state.

The Fortran driver links the unchanged, separately cached attributed source.
Periodic scalar references use one-channel clocur: upstream percur has an
overlapping coefficient assignment bug for short cyclic bases.
Upstream closed-curve restart failures after interpolation are recorded before
an independent task-0 restart; no upstream source patch is made.
"""
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.interpolate import BSpline, splrep


def main():
    exe = Path(sys.argv[1]).resolve()
    rng = np.random.default_rng(20261121)
    cases = []
    rows = []
    for k in range(1, 6):
        for periodic in (False, True):
            for ndim in (1, 3):
                for m in (47, 91):
                    u = np.r_[0., np.cumsum(rng.uniform(.7, 1.3, m-1))]
                    u /= u[-1]
                    x = np.column_stack([np.cos(2*np.pi*u+j)+.12*np.sin((26+j*4)*np.pi*u)
                                         + .01*rng.normal(size=m) for j in range(ndim)])
                    x[-1] = x[0]
                    w = rng.uniform(.6, 1.4, m)
                    active = m-int(periodic)
                    energy = np.sum((x[:active]*w[:active, None])**2)
                    targets = energy*np.array([10., .0005, .001, .002, .0002, 10., 0., .0005])
                    cases.append((u, x, w, k, periodic, targets))
                    rows.append(f'{m} {ndim} {k} {int(periodic)} {len(targets)}')
                    rows.append(' '.join(f'{v:.17g}' for v in u))
                    rows.extend(' '.join(f'{v:.17g}' for v in row) for row in x)
                    rows.extend(' '.join(f'{v:.17g}' for v in row) for row in (w, targets))
    data = '\n'.join(rows)+'\n'
    def run(reference):
        proc = subprocess.run([str(exe), reference], input=data, text=True,
                              capture_output=True, check=True)
        return iter(proc.stdout.splitlines())
    native, upstream = run('native'), run('reference-cyclic')
    fits = fixes = history = 0
    coefficient_max = residual_max = derivative_max = seam_max = 0.
    for index, (u, x, w, k, periodic, targets) in enumerate(cases):
        for step, target in enumerate(targets):
            code, _ = map(int, next(native).split())
            refcode, original = map(int, next(upstream).split())
            assert code == refcode == 0, (index, step, k, periodic, code, refcode, original)
            fixes += int(original > 0)
            n = int(next(native)); nr = int(next(upstream))
            assert n == nr, (index, step, 'knot count', n, nr)
            t = np.fromstring(next(native), sep=' ')
            tr = np.fromstring(next(upstream), sep=' ')
            np.testing.assert_allclose(t, tr, rtol=2e-13, atol=2e-13)
            c = np.column_stack([np.fromstring(next(native), sep=' ') for _ in range(x.shape[1])])
            cr = np.column_stack([np.fromstring(next(upstream), sep=' ') for _ in range(x.shape[1])])
            error = np.max(abs(c-cr)/(1+abs(cr)))
            coefficient_max = max(coefficient_max, error)
            assert error < 2e-7, (index, step, 'coefficients', error)
            bs = BSpline(t, c, k); ref = BSpline(tr, cr, k)
            active = len(u)-int(periodic)
            fp = np.sum((w[:active, None]*(bs(u[:active])-x[:active]))**2)
            assert fp <= 1.0011*target+1e-22, (index, step, 'residual', fp, target)
            if target > 0. and step not in (0, 5):
                err = abs(fp-target)/target
                residual_max = max(residual_max, err)
                assert err < .0011, (index, step, 'active residual', err)
            q = np.linspace(-.05, 1.05, 31)
            for d in range(k+1):
                err = np.max(abs(bs(q, nu=d)-ref(q, nu=d))/(1+abs(ref(q, nu=d))))
                derivative_max = max(derivative_max, err)
                assert err < 2e-6, (index, step, 'derivative', d, err)
                if periodic and d < k:
                    ends = bs([0., 1.], nu=d)
                    scale = max(1., np.max(abs(bs.derivative(d).c)))
                    err = np.max(abs(ends[1]-ends[0]))/scale
                    seam_max = max(seam_max, err)
                    assert err < 4096*np.finfo(float).eps, (index, step, 'seam', d, err)
            if x.shape[1] == 1 and target > 0.:
                fresh = splrep(u, x[:, 0], w=w, k=k, s=target, per=periodic)
                history += int(len(fresh[0]) != len(t) or
                               not np.allclose(fresh[0], t, rtol=2e-13, atol=2e-13))
            fits += 1
    assert next(native, None) is None and next(upstream, None) is None
    print(f'PASS: {len(cases)} sequences; {fits} continued/fresh fits; '
          f'{fixes} recorded upstream closed restart failures recovered independently')
    print(f'Max normalized coefficient difference {coefficient_max:.6g}; '
          f'derivative difference {derivative_max:.6g}; relative residual error {residual_max:.6g}; '
          f'seam {seam_max:.6g}')
    print(f'{history} scalar fits have different knots from independent fresh SciPy fits')


if __name__ == '__main__':
    main()
