#!/usr/bin/env python3
"""Independent SciPy CubicSpline comparison for matrix ordinates (axis=0)."""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import CubicSpline

ROOT=Path(__file__).resolve().parents[2]


def cases(seed):
    rng=np.random.default_rng(seed)
    for flavor in range(6):
        for n in (2,3,4,7,20):
            for m in (1,2,4,7):
                for pattern in range(4):
                    h=10.0**rng.uniform(-0.9,0.7,size=n-1)
                    x=np.r_[0,np.cumsum(h)]-rng.uniform(-2,2)
                    if pattern == 0:
                        y=rng.normal(size=(n,m))
                    elif pattern == 1:
                        y=rng.normal(size=(1,m))+np.outer(x,rng.normal(size=m))
                    elif pattern == 2:
                        y=np.sin(np.outer(x,rng.uniform(0.2,1.5,m)))
                    else:
                        z=x/max(1.0,x[-1]-x[0])
                        y=rng.normal(size=(1,m))*z[:,None]**3+np.outer(z,rng.normal(size=m))
                    if flavor==3:
                        y[-1,:]=y[0,:]
                    span=x[-1]-x[0]
                    q=np.r_[x, (x[:-1]+x[1:])/2, rng.uniform(x[0]-1.3*span,x[-1]+1.1*span,7)]
                    a=float(rng.uniform(x[0]-1.5*span,x[-1]+1.5*span))
                    b=float(rng.uniform(x[0]-1.5*span,x[-1]+1.5*span))
                    left=rng.normal(size=m)
                    right=rng.normal(size=m)
                    ext=(-1,0,1)[(n+m+flavor+pattern)%3]
                    yield n,m,q,flavor,ext,x,y,a,b,left,right


def bc_for(flavor,left,right):
    return ['not-a-knot','natural','clamped','periodic',
            ((1,left),(2,right)),((2,left),(1,right))][flavor]


def run(flags,seed):
    trials=list(cases(seed))
    with tempfile.TemporaryDirectory(prefix='scifort_cubic_spline_multi_') as tmp:
        tmp=Path(tmp)
        sources=[ROOT/'src/kinds.f90',ROOT/'src/interpolate_bspline.f90',ROOT/'src/interpolate_regular_grid.f90',ROOT/'src/interpolate_piecewise.f90', ROOT/'src/interpolate_barycentric.f90', ROOT/'src/interpolate_global.f90', ROOT/'src/interpolate_kdtree.f90', ROOT/'src/interpolate_rbf.f90', ROOT/'src/interpolate_legacy.f90', ROOT/'src/interpolate_scattered.f90', ROOT/'src/interpolate_clough_tocher.f90', ROOT/'src/interpolate_fitpack_bivariate.f90', ROOT/'src/interpolate_fitpack_compat.f90', ROOT/'src/interpolate.f90',
                 ROOT/'tools/interpolate_reference/cubic_spline_multi_driver.f90']
        cmd=[os.environ.get('FC','gfortran'),'-std=f2018',*flags,'-J'+str(tmp),'-I'+str(tmp),
             *map(str,sources),'-o',str(tmp/'driver')]
        proc=subprocess.run(cmd,cwd=tmp,capture_output=True,text=True)
        if proc.returncode: raise RuntimeError(proc.stderr)
        inp=[str(len(trials))]
        for n,m,q,flavor,ext,x,y,a,b,left,right in trials:
            inp.append(f'{n} {m} {len(q)} {flavor} {ext}')
            inp.append(' '.join(format(v,'.17g') for v in x))
            for col in y.T:
                inp.append(' '.join(format(v,'.17g') for v in col))
            inp.append(' '.join(format(v,'.17g') for v in q))
            inp.append(f'{a:.17g} {b:.17g}')
            inp.append(' '.join(format(v,'.17g') for v in left))
            inp.append(' '.join(format(v,'.17g') for v in right))
        proc=subprocess.run([str(tmp/'driver')],input='\n'.join(inp)+'\n',
                            capture_output=True,text=True)
        if proc.returncode: raise RuntimeError(proc.stderr[-3000:]+proc.stdout[-500:])
        lines=iter(proc.stdout.splitlines())
        nvalues=0
        maxabs=0.0
        errs=[]
        for index,(n,m,q,flavor,ext,x,y,a,b,left,right) in enumerate(trials,1):
            assert next(lines).split()==['CASE',str(index),'0']
            actual=np.array([[[float(v) for v in next(lines).split()]
                               for _ in q] for _ in range(5)])
            integ=np.array([float(v) for v in next(lines).split()])
            coeff=np.empty((4,n-1,m))
            for k in range(m):
                for j in range(n-1):
                    coeff[:,j,k]=[float(v) for v in next(lines).split()]
            ref=CubicSpline(x,y,axis=0,bc_type=bc_for(flavor,left,right),
                            extrapolate=None if ext==-1 else bool(ext))
            expected=np.array([ref(q,nu=nu) for nu in range(5)])
            for label,got,want in [('values',actual,expected),
                                   ('integral',integ,ref.integrate(a,b)),
                                   ('coefficients',coeff,ref.c[::-1])]:
                good=np.isclose(got,want,atol=1e-8,rtol=2e-9,equal_nan=True)
                if not good.all():
                    ix=tuple(np.argwhere(~good)[0])
                    errs.append((index,flavor,n,m,label,ix,got[ix],want[ix]))
                finite=np.isfinite(got)&np.isfinite(want)
                if finite.any(): maxabs=max(maxabs,float(abs(got[finite]-want[finite]).max()))
                nvalues+=got.size
        print(f'{len(trials)} SciPy matrix spline cases, {nvalues} values, '
              f'maximum absolute discrepancy {maxabs:.5g}')
        if errs:
            print(*errs[:12],sep='\n')
            raise AssertionError(f'{len(errs)} numerical discrepancies')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=20261031)
    parser.add_argument('--checked-only',action='store_true')
    args=parser.parse_args()
    run(['-O0','-g','-fcheck=all','-Wall','-Wextra'],args.seed)
    if not args.checked_only:
        run(['-O3','-march=native'],args.seed)
