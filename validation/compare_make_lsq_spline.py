#!/usr/bin/env python3
"""Compare SciFort's weighted complete-knot LSQ B-splines to SciPy 1.17+.

Run from anywhere: python validation/compare_make_lsq_spline.py
Requires NumPy, SciPy, and gfortran; only files in a temporary directory
are compiled.  SciPy is an independent numerical reference, not a dependency
of SciFort's Fortran library.
"""
from pathlib import Path
import argparse
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import BSpline, make_lsq_spline

ROOT = Path(__file__).resolve().parents[1]


def case(rng, degree):
    m = max(25, 2 * degree + 18) + int(rng.integers(0, 29))
    x = np.linspace(-1.0, 1.0, m)
    x[1:-1] += rng.uniform(-0.002, 0.002, size=m-2)
    x.sort()
    if rng.random() < 0.2:
        x[m//3] = x[m//3-1]
    ninterior = int(rng.integers(0, 5))
    interior = np.sort(rng.uniform(-0.65, 0.65, ninterior))
    if degree > 0 and interior.size > 1 and rng.random() < 0.3:
        interior[1] = interior[0]
    t = np.r_[np.repeat(-1.0, degree + 1), interior, np.repeat(1.0, degree + 1)]
    w = rng.uniform(0.4, 1.5, size=m)
    if rng.random() < 0.5:
        w[rng.integers(1, m-1)] = 0.0
    nch = int(rng.integers(1, 5))
    y = np.empty((m, nch), dtype=float, order='F')
    for ch in range(nch):
        y[:, ch] = (0.8 * np.sin((ch+2.5)*x) + 0.15*np.cos((ch+1)*x)
                    + (ch+1)*0.02*x**2 + 0.025*rng.normal(size=m))
    return x, y, t, degree, w


def run(seed, count):
    rng = np.random.default_rng(seed)
    cases = []
    for degree in rng.permutation(np.tile(np.arange(16), (count+15)//16))[:count]:
        # SciPy's QR solve may return inf without reporting a singular
        # matrix when a sample-free knot interval causes rank deficiency.
        # Compare only valid full-rank fits; rank rejection is tested natively.
        for attempt in range(200):
            item = case(rng, int(degree))
            x, y, t, k, w = item
            a = BSpline.design_matrix(x,t,k).toarray()*w[:,None]
            if np.linalg.matrix_rank(a, tol=2e-12*np.linalg.norm(a,2)) == a.shape[1]:
                cases.append(item)
                break
        else:
            raise RuntimeError(f'Unable to generate a full-rank case for k={degree}')
    references = []
    data = [str(count)]
    for x, y, t, degree, w in cases:
        r = make_lsq_spline(x, y, t, k=degree, w=w, method='qr')
        references.append(r.c)
        data.append(f'{len(x)} {degree} {y.shape[1]} {len(t)}')
        for a in (x,t,w,y.ravel(order='F')):
            data.append(' '.join(f'{v:.17e}' for v in a))
    with tempfile.TemporaryDirectory(prefix='scifort-lsq-') as tmp:
        out = Path(tmp)/'lsq-compare'
        cmd = ['gfortran','-std=f2018','-O0','-g','-fcheck=all',
               '-ffree-line-length-none','-J',tmp,'-I',tmp,
               str(ROOT/'src/kinds.f90'),str(ROOT/'src/interpolate_bspline.f90'),
               str(ROOT/'validation/make_lsq_driver.f90'),'-o',str(out)]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        proc = subprocess.run([str(out)], input='\n'.join(data)+'\n',
                              capture_output=True,text=True,check=True)
    lines = iter(proc.stdout.strip().splitlines())
    max_abs, max_rel = 0.0, 0.0
    for i, ref in enumerate(references):
        code = int(next(lines))
        assert code == 0, f'Case {i} (k={cases[i][3]}): Fortran status {code}'
        got = np.array([[float(v) for v in next(lines).split()]
                        for _ in range(ref.shape[1])]).T
        assert got.shape == ref.shape
        err = float(np.max(np.abs(got-ref)))
        rel = err/max(1.0, float(np.max(np.abs(ref))))
        max_abs = max(max_abs,err)
        max_rel = max(max_rel,rel)
        assert rel < 3e-7, (f'Case {i} (k={cases[i][3]}): coefficients disagree '
                                 f'abs={err:.4e}, relative={rel:.4e}')
    print(f'PASS: {count} SciPy QR reference fits, seed={seed}, '
          f'max absolute={max_abs:.3e}, max relative={max_rel:.3e}')


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=20261101)
    parser.add_argument('--count',type=int,default=160)
    a=parser.parse_args()
    run(a.seed,a.count)
