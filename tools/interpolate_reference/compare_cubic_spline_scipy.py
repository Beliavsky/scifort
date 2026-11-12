#!/usr/bin/env python3
"""Compare independent Fortran CubicSpline with SciPy 1.17.0.

Compiles Fortran with gfortran in temporary directories, without editing source.
Covers periodic, not-a-knot, natural, clamped and asymmetric first/second BCs.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[2]


def cases(seed: int):
    rng = np.random.default_rng(seed)
    for flavor in range(11):
        for n in (2, 3, 4, 5, 8, 17, 49):
            for pattern in range(11):
                if pattern in (0, 3, 7):
                    h = 10.0**rng.uniform(-1.6, 1.3, n-1)
                else:
                    h = rng.uniform(0.15, 2.7, n-1)
                x = np.r_[0.0, np.cumsum(h)] - rng.uniform(-4.0, 4.0)
                if pattern == 0:
                    y = np.full(n, rng.uniform(-3, 3))
                elif pattern == 1:
                    y = rng.normal() + rng.normal()*x
                elif pattern == 2:
                    y = (x-np.mean(x))**2
                elif pattern == 3:
                    z = (x-x.mean()) / max(1.0, 0.25*(x[-1]-x[0]))
                    y = 0.2*z**3-0.3*z**2+0.1*z
                elif pattern == 4:
                    y = np.sin(1.4*x)
                elif pattern == 5:
                    y = np.cumsum(rng.uniform(-1.0, 1.0, n))
                elif pattern == 6:
                    y = rng.integers(-3, 4, n).astype(float)
                elif pattern == 7:
                    y = rng.normal(0, 0.003, n)
                elif pattern == 8:
                    y = np.zeros(n)
                    y[n//2:] = 1
                elif pattern == 9:
                    y = rng.normal(size=n)
                else:
                    y = np.cos(0.4*x)+0.1*np.sin(2.0*x)
                if flavor == 3:
                    y[-1] = y[0]
                span = x[-1]-x[0]
                q = np.r_[x, (x[:-1]+x[1:])*0.5, x[0]-0.9*span, x[-1]+0.7*span,
                          rng.uniform(x[0]-1.4*span, x[-1]+1.7*span, 11)]
                a = float(rng.uniform(x[0]-2.2*span, x[-1]+1.8*span))
                b = float(rng.uniform(x[0]-1.8*span, x[-1]+2.2*span))
                lv = float(rng.normal())
                rv = float(rng.normal())
                ext = (-1, 0, 1)[(pattern+flavor+n)%3]
                yield flavor, n, x, y, q, ext, a, b, lv, rv


def scipy_bc(flavor: int, lv: float, rv: float):
    return [
        'not-a-knot', 'natural', 'clamped', 'periodic',
        ((1, lv), 'natural'), ('not-a-knot', (1, rv)),
        ((2, lv), (2, rv)), ('natural', 'clamped'),
        ('clamped', 'natural'), ((1, lv), (1, rv)),
        ((2, lv), 'not-a-knot'),
    ][flavor]


def run(flags: list[str], seed: int) -> None:
    trials = list(cases(seed))
    with tempfile.TemporaryDirectory(prefix='scifort_cubic_spline_') as name:
        path = Path(name)
        sources = [ROOT/'src/kinds.f90', ROOT/'src/interpolate_bspline.f90',ROOT/'src/interpolate_regular_grid.f90',ROOT/'src/interpolate_piecewise.f90', ROOT/'src/interpolate_barycentric.f90', ROOT/'src/interpolate_global.f90', ROOT/'src/interpolate_kdtree.f90', ROOT/'src/interpolate_rbf.f90', ROOT/'src/interpolate_legacy.f90', ROOT/'src/interpolate_scattered.f90', ROOT/'src/interpolate_clough_tocher.f90', ROOT/'src/interpolate_fitpack_bivariate.f90', ROOT/'src/interpolate_fitpack_compat.f90', ROOT/'src/interpolate.f90',
                   ROOT/'tools/interpolate_reference/cubic_spline_driver.f90']
        compiler = os.environ.get('FC', 'gfortran')
        proc = subprocess.run([compiler,'-std=f2018',*flags,'-J'+str(path),'-I'+str(path),
                               *map(str,sources),'-o',str(path/'driver')],cwd=path,
                              capture_output=True,text=True)
        if proc.returncode:
            raise RuntimeError(proc.stderr)
        payload=[str(len(trials))]
        for flavor, n, x, y, q, ext, a, b, lv, rv in trials:
            payload.extend((f'{n} {len(q)} {flavor} {ext}',
                           ' '.join(format(v,'.17g') for v in x),
                           ' '.join(format(v,'.17g') for v in y),
                           ' '.join(format(v,'.17g') for v in q),
                           f'{a:.17g} {b:.17g} {lv:.17g} {rv:.17g}'))
        proc = subprocess.run([str(path/'driver')],input='\n'.join(payload)+'\n',
                              capture_output=True,text=True)
        if proc.returncode:
            raise RuntimeError(f'Fortran driver failed:\n{proc.stderr}\n{proc.stdout[-2000:]}')
        lines=iter(proc.stdout.splitlines())
        errs=[]
        maximum=0.0
        max_scaled=0.0
        comparisons=0
        for number, (flavor,n,x,y,q,ext,a,b,lv,rv) in enumerate(trials,1):
            header=next(lines).split()
            if header != ['CASE',str(number),'0']:
                raise AssertionError((number,header))
            actual=np.array([[float(z) for z in next(lines).split()] for _ in q])
            integral=float(next(lines))
            coeff=np.array([[float(z) for z in next(lines).split()]
                            for _ in range(n-1)]).T
            extr=None if ext==-1 else bool(ext)
            ref=CubicSpline(x,y,bc_type=scipy_bc(flavor,lv,rv),extrapolate=extr)
            expected=np.array([ref(q,nu=nu) for nu in range(5)]).T
            for label,got,want,rtol,atol in [
                ('values',actual,expected,8e-10,1e-8),
                ('coeff',coeff,ref.c[::-1],8e-10,1e-8),
                ('integral',np.array([integral]),np.array([ref.integrate(a,b)]),5e-9,1e-7),
            ]:
                ok=np.isclose(got,want,rtol=rtol,atol=atol,equal_nan=True)
                if not np.all(ok):
                    ix=tuple(np.argwhere(~ok)[0])
                    errs.append((number, flavor, n, label, ix, got[ix],want[ix],
                                 x.tolist()[:4], y.tolist()[:4], ext))
                finite=np.isfinite(got)&np.isfinite(want)
                if np.any(finite):
                    maximum=max(maximum,float(np.max(abs(got[finite]-want[finite]))))
                    max_scaled=max(max_scaled,float(np.max(abs(got[finite]-want[finite])/(1.0+abs(want[finite])))))
                comparisons+=got.size
        print(f'{len(trials)} SciPy CubicSpline cases, {comparisons} comparison values; '
              f'maximum absolute discrepancy {maximum:.4g}; max |error|/(1+|ref|) {max_scaled:.4g}')
        if errs:
            print('MISMATCHES:',*errs[:12],sep='\n  ')
            raise AssertionError(f'{len(errs)} mismatches')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=20261031)
    parser.add_argument('--checked-only',action='store_true')
    args=parser.parse_args()
    run(['-O0','-g','-fcheck=all','-Wall','-Wextra'],args.seed)
    if not args.checked_only:
        run(['-O3','-march=native'],args.seed)
