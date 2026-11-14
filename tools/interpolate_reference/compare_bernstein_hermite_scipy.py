#!/usr/bin/env python3
"""Compare BPoly.from_derivatives, including varying degrees, with SciPy."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import BPoly

ROOT=Path(__file__).resolve().parents[2]


def examples(seed,number):
    rng=np.random.default_rng(seed)
    cases=[]
    for k in range(number):
        n=2+k%5
        nc=1+k%4
        counts=rng.integers(1,5,size=n)
        if k%5==0:
            counts[:]=1
        if k%5==1:
            counts[:]=2
        nd=int(counts.max())
        x=np.r_[rng.uniform(-1,1),rng.uniform(.3,1.4,size=n-1).cumsum()]
        x[1:]+=x[0]
        y=rng.normal(0,0.35,size=(n,nd,nc))
        q=np.r_[x, rng.uniform(x[0]-.25,x[-1]+.25,size=5), x[0]-.2, x[-1]+.2]
        ext=k%7!=0
        cases.append((x,counts,y,q,ext))
    return cases


def serialize(cases):
    s=[str(len(cases))]
    def dump(a):
        return ' '.join(f'{float(i):.17g}' for i in np.ravel(a,order='F'))
    for x,counts,y,q,ext in cases:
        s.append(f'{len(x)} {y.shape[1]} {y.shape[2]} {len(q)} {int(ext)}')
        s.extend((dump(x),dump(counts),dump(y),dump(q)))
    return '\n'.join(s)+'\n'


def check(flags,seed,number):
    cases=examples(seed,number)
    with tempfile.TemporaryDirectory(prefix='scifort_hermite_') as td:
        wd=Path(td)
        sources=['src/kinds.f90','src/interpolate_piecewise.f90','src/interpolate_bspline.f90',
                 'src/interpolate_regular_grid.f90','src/interpolate_barycentric.f90', 'src/interpolate_global.f90', 'src/interpolate_kdtree.f90', 'src/interpolate_rbf.f90', 'src/interpolate_legacy.f90', 'src/interpolate_scattered.f90', 'src/interpolate_clough_tocher.f90', 'src/interpolate_fitpack_bivariate.f90', 'src/interpolate_curfit_core.f90', 'src/interpolate_fitpack_compat.f90', 'src/interpolate.f90',
                 'tools/interpolate_reference/bernstein_hermite_driver.f90']
        cmd=['gfortran','-std=f2018','-ffree-line-length-none',*flags,'-J'+str(wd),'-I'+str(wd),
             *(str(ROOT/s) for s in sources),'-o',str(wd/'test')]
        cp=subprocess.run(cmd,capture_output=True,text=True,cwd=wd)
        if cp.returncode:
            raise RuntimeError('Compilation failed:\n'+cp.stderr[-4000:])
        cp=subprocess.run([str(wd/'test')],input=serialize(cases),text=True,capture_output=True,cwd=wd)
        if cp.returncode:
            raise RuntimeError('Fortran run failed:\n'+cp.stderr[-4000:])
        lines=iter(cp.stdout.splitlines())
        nchecks=0
        max_error=0.0
        for caseid,(x,counts,y,q,ext) in enumerate(cases,1):
            data=[[list(y[i,j]) for j in range(counts[i])] for i in range(len(x))]
            b=BPoly.from_derivatives(x,data,extrapolate=ext)
            header=next(lines).split()
            assert header==['CASE',str(caseid),str(b.c.shape[0])], (header,b.c.shape)
            coef=np.fromstring(next(lines),sep=' ').reshape(b.c.shape,order='F')
            checks=[(coef,b.c,'coefficients')]
            for nu in range(4):
                native=np.array([np.fromstring(next(lines),sep=' ') for _ in q])
                if nu>b.c.shape[0]:
                    from scipy.interpolate import PPoly
                    expected=PPoly(np.zeros((1,len(x)-1,y.shape[-1])),x,extrapolate=ext)(q)
                else:
                    expected=b(q,nu=nu)
                checks.append((native,expected,f'derivative {nu}'))
            for a,z in [(x[0],x[-1]),(x[-1],x[0])]:
                native=np.fromstring(next(lines),sep=' ')
                checks.append((native,b.integrate(a,z),'integration'))
            for actual,expected,kind in checks:
                actual=np.asarray(actual)
                expected=np.asarray(expected)
                if not np.allclose(actual,expected,rtol=2e-9,atol=2e-10,equal_nan=True):
                    index=tuple(np.argwhere(~np.isclose(actual,expected,rtol=2e-9,atol=2e-10,equal_nan=True))[0])
                    raise AssertionError(f'case {caseid}: {kind} at {index}, '
                                         f'got {actual[index]} expected {expected[index]}')
                nchecks+=actual.size
                if np.any(np.isfinite(expected)):
                    max_error=max(max_error,float(np.nanmax(np.abs(actual-expected))))
        assert next(lines,None) is None
        print(f'{flags}: {len(cases)} Hermite cases, {nchecks} entries, maximum abs error {max_error:.3g}')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=20261031)
    parser.add_argument('--cases',type=int,default=140)
    args=parser.parse_args()
    for flags in [['-O0','-g','-fcheck=all','-fbacktrace'],['-O3']]:
        check(flags,args.seed,args.cases)
