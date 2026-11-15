# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Check short cyclic scalar bases against SciPy and independent B-spline evaluation."""
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.interpolate import BSpline, splrep


def main():
    exe = Path(sys.argv[1]).resolve()
    rng = np.random.default_rng(20261122)
    cases, rows = [], []
    for m in (31, 61, 97):
        for k in range(1, 6):
            u = np.linspace(0, 1, m)**1.1
            y = np.cos(2*np.pi*u)+.08*np.sin(14*np.pi*u)
            if m != 61:
                y += .005*rng.normal(size=m)
            y[-1] = y[0]
            w = .8+u
            energy = np.sum((w[:-1]*y[:-1])**2)
            for factor in (.01, .03, .1, .2, .6):
                target = factor*energy
                targets = target*np.array([1., .75, 1.25])
                cases.append((u, y, w, k, targets))
                rows.append(f'{m} 1 {k} 1 3')
                rows.append(' '.join(map(repr, u.tolist())))
                rows.extend(map(repr, y.tolist()))
                rows.extend(' '.join(map(repr, a.tolist())) for a in (w, targets))
    proc = subprocess.run([str(exe), 'native'], input='\n'.join(rows)+'\n',
                          text=True, capture_output=True, check=True)
    lines = iter(proc.stdout.splitlines())
    fits = short = 0
    coefficient_max = residual_max = seam_max = 0.
    for u, y, w, k, targets in cases:
        for step, target in enumerate(targets):
            code, _ = map(int, next(lines).split())
            assert code == 0, (len(u), k, target, step, code)
            n = int(next(lines))
            t = np.fromstring(next(lines), sep=' ')
            c = np.fromstring(next(lines), sep=' ')
            n7 = n-2*k-1
            short += int(n7 < k)
            bs = BSpline(t, c, k)
            fp = np.sum((w[:-1]*(bs(u[:-1])-y[:-1]))**2)
            error = abs(fp-target)/target
            residual_max = max(residual_max, error)
            assert error < .0011, (len(u), k, step, fp, target)
            for i in range(k):
                assert c[n7+i] == c[i]
            for d in range(k):
                ends = bs([0., 1.], nu=d)
                scale = max(1., np.max(abs(bs.derivative(d).c)))
                error = abs(ends[1]-ends[0])/scale
                seam_max = max(seam_max, error)
                assert error < 4096*np.finfo(float).eps
            if step == 0:
                # Match the wrapper's scaling to avoid knot/root branch roundoff.
                ys = np.max(abs(y[:-1])); ws = np.max(w[:-1])
                reference = splrep(u, y/ys, w=w/ws, k=k, per=True,
                                   s=np.exp(np.log(target)-2*np.log(ys)-2*np.log(ws)),
                                   full_output=True)
                (tr, cr, _), _, ier, _ = reference
                assert ier <= 0, (len(u), k, ier)
                np.testing.assert_allclose(t, tr, rtol=2e-13, atol=2e-13)
                error = np.max(abs(c-cr[:len(c)]*ys)/(1+abs(c)))
                coefficient_max = max(coefficient_max, error)
                assert error < 2e-7
            fits += 1
    assert next(lines, None) is None
    print(f'PASS: {fits} fits, {short} short cyclic bases; 75 fresh SciPy comparisons')
    print(f'Max coefficient difference {coefficient_max:.6g}; '
          f'relative residual error {residual_max:.6g}; seam {seam_max:.6g}')


if __name__ == '__main__':
    main()
