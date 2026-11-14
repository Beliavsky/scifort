#!/usr/bin/env python3
"""Test PPoly.solve, PPoly.roots and PPoly.from_spline against SciPy 1.17.

Run checked and optimized tests, two independent seeded random generators.
Usage: python tools/interpolate_reference/compare_piecewise_roots_scipy.py
"""
from __future__ import annotations
import argparse
import subprocess
from pathlib import Path
import tempfile
import numpy as np
from scipy.interpolate import PPoly, BSpline

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ['src/kinds.f90','src/interpolate_piecewise.f90',
           'src/interpolate_bspline.f90','src/interpolate_regular_grid.f90','src/interpolate_barycentric.f90', 'src/interpolate_global.f90', 'src/interpolate_kdtree.f90', 'src/interpolate_rbf.f90', 'src/interpolate_legacy.f90', 'src/interpolate_scattered.f90', 'src/interpolate_clough_tocher.f90', 'src/interpolate_fitpack_bivariate.f90', 'src/interpolate_curfit_core.f90', 'src/interpolate_fitpack_compat.f90', 'src/interpolate.f90']


def compile_driver(wd, name, flags):
    args = ['gfortran','-std=f2018','-ffree-line-length-none', *flags,
            '-J'+str(wd),'-I'+str(wd),*(str(ROOT/f) for f in SOURCES),
            str(ROOT/'tools/interpolate_reference'/f'{name}.f90'),'-o',str(wd/name)]
    cp=subprocess.run(args,cwd=wd,text=True,capture_output=True)
    if cp.returncode: raise RuntimeError(cp.stderr[-4500:])
    return wd/name


def fmt(a):
    return ' '.join(f'{float(v):.17g}' for v in np.ravel(a,order='F'))


def roots_cases(seed, n):
    rng=np.random.default_rng(seed)
    out=[]
    for i in range(n):
        degree=i%9
        nseg=1+(i%6)
        nch=1+(i%3)
        x=np.r_[0.,np.cumsum(rng.uniform(.3,1.8,nseg))]
        x+=rng.normal()
        if i%7 == 2: x=x[::-1].copy()
        c=rng.uniform(-2,2,(degree+1,nseg,nch))
        if i%4 == 0:
            c[:-1]*=0.15
            c[-1] = rng.normal(size=(nseg,nch))
        if i%11 == 0:
            c[:,(i//11)%nseg,0]=0.0
        if i%13 == 0 and degree >= 2:
            # Tangential double root, translated to a random local knot.
            loc=.2*(x[1]-x[0]); c[:,0,0]=0
            c[-3,0,0]=1.; c[-2,0,0]=-2*loc; c[-1,0,0]=loc*loc
        if i%19 == 0 and degree >= 1:
            # Root at an exact breakpoint.
            c[-1,0,0]=0.
        # Zero-length intervals are tested by a separate deterministic suite;
        # no zero spans are created here.
        y=0. if i%3 else float(rng.uniform(-1,1))
        ext=(i%4 != 1)
        periodic=(i%17 == 1)
        discontinuity=(i%4 != 2)
        out.append((x,c,y,ext,periodic,discontinuity))
    # Degenerate interior/end intervals are legal in PPoly and appear in
    # BSpline-to-PPoly conversions with repeated knot multiplicities.
    for x_values in ([0.,0.,1.],[0.,1.,1.,2.],[2.,1.,1.,0.],
                     [0.,0.,.5,1.,1.,2.]):
        x=np.asarray(x_values)
        for degree in range(5):
            for ext in (True,False):
                c=rng.uniform(-1,1,(degree+1,len(x)-1,1))
                if degree==0: c[0,0,0]=0.
                if degree==1: c[-1,:,0]=0.
                out.append((x.copy(),c,0.,ext,False,True))
    return out


def roots_serial(cases):
    lines=[str(len(cases))]
    for x,c,y,ext,periodic,discontinuity in cases:
        lines += [f'{len(c)-1} {len(x)-1} {c.shape[-1]} {int(ext)} {int(periodic)} {int(discontinuity)} {y:.17g}',fmt(x),fmt(c)]
    return '\n'.join(lines)+'\n'


def compare_roots(binary,cases):
    cp=subprocess.run([str(binary)],input=roots_serial(cases),capture_output=True,text=True)
    if cp.returncode: raise RuntimeError(cp.stderr[-2500:])
    lines=iter(cp.stdout.splitlines())
    count=0
    err=0.
    for idx,(x,c,y,ext,periodic,disc) in enumerate(cases,1):
        assert next(lines) == f'CASE {idx}'
        for j in range(c.shape[-1]):
            head=next(lines).split()
            assert head[:2] == ['ROOTS',str(j+1)],head
            count_values=int(head[2]); line=next(lines)
            actual=np.fromstring(line,sep=' ')
            assert len(actual)==count_values,(idx,j,count_values,actual)
            mode='periodic' if periodic else ext
            expected=PPoly(c[:,:,j],x,extrapolate=mode).solve(y=y,discontinuity=disc)
            if actual.shape != expected.shape:
                raise AssertionError(f'case={idx}, channel={j}, shape {actual} != {expected}; data c={c[:,:,j]}, x={x}, rhs={y}')
            # SciPy's polynomial-root ordering within a section is not
            # guaranteed, notably when roots lie outside the base interval.
            a_sort=np.sort(actual[np.isfinite(actual)])
            e_sort=np.sort(expected[np.isfinite(expected)])
            if (len(a_sort)!=len(e_sort) or np.isnan(actual).sum()!=np.isnan(expected).sum() or
                not np.allclose(a_sort,e_sort,rtol=1e-7,atol=1e-8)):
                raise AssertionError(f'case={idx}, channel={j}, actual {actual}, expected {expected}; x={x}, c={c[:,:,j]}, rhs={y}')
            count+=len(actual)
            if len(a_sort):
                err=max(err,float(np.max(np.abs(a_sort-e_sort))))
    assert next(lines,None) is None
    return count,err


def spline_cases(seed,n):
    rng=np.random.default_rng(seed)
    out=[]
    for i in range(n):
        k=i%7
        nc=2*k+1+(i%7)
        # Generate valid knot vectors with k+1 repeated endpoints and optional
        # internal multiplicity 2 or 3, as long as enough coefficients remain.
        internal_count=nc-k-1
        base=np.cumsum(rng.uniform(.2,1.2,internal_count+1))
        internal=base[:-1]
        if k>1 and internal_count>3 and i%3==1:
            internal[1]=internal[0]
        if k>2 and internal_count>4 and i%4==1:
            internal[3]=internal[2]
        t=np.r_[np.repeat(base[0]-0.4,k+1),internal,np.repeat(base[-1]+0.4,k+1)]
        # There are nc coefficients when len(t)=nc+k+1, and interior count nc-k-1.
        assert len(t)==nc+k+1,(len(t),nc,k)
        c=rng.normal(0,1,(nc,1+i%3))
        q=rng.uniform(t[k]-.5,t[nc]+.5,12)
        q=np.r_[q,t[k],t[nc],t[0],t[-1],t[1:-1]]
        out.append((t,c,k,q))
    return out


def spline_serial(cases):
    lines=[str(len(cases))]
    for t,c,k,q in cases:
        lines.extend([f'{k} {len(c)} {c.shape[-1]} {len(q)}',fmt(t),fmt(c),fmt(q)])
    return '\n'.join(lines)+'\n'


def compare_splines(binary,cases):
    cp=subprocess.run([str(binary)],input=spline_serial(cases),capture_output=True,text=True)
    if cp.returncode: raise RuntimeError(cp.stderr[-2000:])
    lines=iter(cp.stdout.splitlines())
    count=0; err=0.
    for idx,(t,c,k,q) in enumerate(cases,1):
        assert next(lines)==f'CASE {idx}'
        h=next(lines).split(); assert h[0]=='COEFF'; shape=tuple(map(int,h[1:]))
        a=np.fromstring(next(lines),sep=' ').reshape(shape,order='F')
        h=next(lines).split(); assert h[0]=='VALUES'; shape_v=tuple(map(int,h[1:]))
        av=np.fromstring(next(lines),sep=' ').reshape(shape_v,order='F')
        assert a.shape==(k+1,len(t)-1,len(c[0]))
        for j in range(c.shape[-1]):
            src=BSpline(t,c[:,j],k)
            ref=PPoly.from_spline(src)
            # Numerically unstable high-order coefficient representations
            # get a scale-sensitive tolerance; see README.
            if not np.allclose(a[:,:,j],ref.c,rtol=1e-7,atol=5e-7):
                d=np.abs(a[:,:,j]-ref.c); at=np.unravel_index(np.argmax(d),d.shape)
                raise AssertionError(f'spline={idx}, channel={j}, coefficient at {at}, actual={a[at+(j,)]} expected={ref.c[at]}')
            expected=ref(q)
            if not np.allclose(av[:,j],expected,rtol=1e-8,atol=2e-7):
                d=np.abs(av[:,j]-expected); at=np.argmax(d)
                raise AssertionError(f'spline={idx}, channel={j}, value at {at}, actual={av[at,j]} expected={expected[at]}')
            count+=a[:,:,j].size+av[:,j].size
            err=max(err,np.max(np.abs(a[:,:,j]-ref.c)),np.max(np.abs(av[:,j]-expected)))
    assert next(lines,None) is None
    return count,err


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=20261031)
    parser.add_argument('--roots',type=int,default=140)
    parser.add_argument('--splines',type=int,default=80)
    parser.add_argument('--optimized-only',action='store_true')
    args=parser.parse_args()
    root_tests=roots_cases(args.seed,args.roots)
    spline_tests=spline_cases(args.seed,args.splines)
    with tempfile.TemporaryDirectory(prefix='scifort_piecewise_roots_') as d:
        wd=Path(d)
        builds=[['-O3']] if args.optimized_only else [['-O0','-g','-fcheck=all','-fbacktrace'],['-O3']]
        for flags in builds:
            r=compile_driver(wd,'ppoly_roots_driver',flags)
            s=compile_driver(wd,'ppoly_from_spline_driver',flags)
            nr,er=compare_roots(r,root_tests)
            ns,es=compare_splines(s,spline_tests)
            print(f'{flags}: {len(root_tests)} root cases / {nr} roots; '
                  f'{len(spline_tests)} spline cases / {ns} entries; '
                  f'max abs root error {er:.3g}; spline error {es:.3g}')

if __name__=='__main__': main()
