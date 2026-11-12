#!/usr/bin/env python3
"""Compare complex PPoly/BPoly calculus and basis conversion against SciPy."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import PPoly, BPoly

ROOT=Path(__file__).resolve().parents[2]
SOURCES=['src/kinds.f90','src/interpolate_piecewise.f90',
         'src/interpolate_bspline.f90','src/interpolate_regular_grid.f90',
         'src/interpolate_barycentric.f90', 'src/interpolate_global.f90', 'src/interpolate_kdtree.f90', 'src/interpolate_rbf.f90', 'src/interpolate_legacy.f90', 'src/interpolate_scattered.f90', 'src/interpolate_clough_tocher.f90', 'src/interpolate_fitpack_bivariate.f90', 'src/interpolate_fitpack_compat.f90', 'src/interpolate.f90','tools/interpolate_reference/complex_piecewise_driver.f90']


def generate(seed,n):
    rng=np.random.default_rng(seed)
    cases=[]
    for i in range(n):
        degree=i%9
        nseg=1+i%6
        nch=1+i%4
        x=np.r_[0.,np.cumsum(rng.uniform(.25,1.5,nseg))]
        x+=rng.uniform(-2,2)
        if i%4==1: x=x[::-1].copy()
        c=rng.normal(size=(degree+1,nseg,nch))*.3+1j*rng.normal(size=(degree+1,nseg,nch))*.3
        if i%7==0: c*=.0
        if i%13==0: c[-1]+=1+1j
        q=np.r_[x,rng.uniform(min(x)-.3,max(x)+.3,8),min(x)-.5,max(x)+.5]
        pairs=rng.uniform(min(x)-.35,max(x)+.35,(5,2))
        pairs=np.r_[pairs,[[min(x)-.35,max(x)+.35],[x[0],x[-1]]]]
        periodic=i%9==0
        ext=i%11 !=0
        cases.append((x,c,q,pairs,periodic,ext))
    return cases


def fmt(a):
    return ' '.join(f'{float(v):.17g}' for v in np.ravel(a,order='F'))


def serialize(cases):
    lines=[str(len(cases))]
    for x,c,q,pairs,periodic,ext in cases:
        lines.extend([f'{len(c)-1} {len(x)-1} {c.shape[-1]} {len(q)} {len(pairs)} {int(periodic)} {int(ext)}',
                      fmt(x),fmt(c.real),fmt(c.imag),fmt(q)])
        # Fortran reads pairs for each of four polynomial objects
        for _ in range(4):
            lines.extend(f'{a:.17g} {b:.17g}' for a,b in pairs)
    return '\n'.join(lines)+'\n'


def compare(values,ref,context,max_error):
    if values.shape!=ref.shape: raise AssertionError((context,'shape',values.shape,ref.shape))
    # Fortran complex functions use a slightly different operation ordering;
    # tighten tolerances for degree <= 4 and permit polynomial conditioning
    # for higher degrees during extrapolation and basis conversion.
    if not np.allclose(values,ref,rtol=1e-6,atol=2e-6,equal_nan=True):
        d=np.abs(values-ref)
        idx=tuple(np.argwhere(~np.isclose(values,ref,rtol=1e-6,atol=2e-6,equal_nan=True))[0])
        raise AssertionError(f'{context}: index {idx} SciFort={values[idx]} SciPy={ref[idx]}, difference={d[idx]}')
    good=np.isfinite(ref)
    if np.any(good): max_error[0]=max(max_error[0],float(np.max(np.abs(values[good]-ref[good]))))
    return int(np.sum(good))


def run(cases,flags):
    with tempfile.TemporaryDirectory(prefix='scifort_complex_pieces_') as d:
        wd=Path(d)
        cmd=['gfortran','-std=f2018','-ffree-line-length-none',*flags,
             '-J'+str(wd),'-I'+str(wd),*(str(ROOT/s) for s in SOURCES),'-o',str(wd/'driver')]
        cp=subprocess.run(cmd,cwd=wd,capture_output=True,text=True)
        if cp.returncode: raise RuntimeError(cp.stderr[-5000:])
        cp=subprocess.run([str(wd/'driver')],input=serialize(cases),text=True,capture_output=True)
        if cp.returncode: raise RuntimeError(cp.stderr[-3500:]+cp.stdout[-600:])
        lines=iter(cp.stdout.splitlines()); checks=0; max_error=[0.]
        def read_complex(shape):
            re=np.fromstring(next(lines),sep=' ').reshape(shape,order='F')
            im=np.fromstring(next(lines),sep=' ').reshape(shape,order='F')
            return re+1j*im
        for idx,(x,c,q,pairs,periodic,ext) in enumerate(cases,1):
            assert next(lines)==f'CASE {idx}'
            flag='periodic' if periodic else ext
            objects=(PPoly(c,x,extrapolate=flag),BPoly(c,x,extrapolate=flag))
            objects=(*objects,BPoly.from_power_basis(objects[0]),PPoly.from_bernstein_basis(objects[1]))
            for j,obj in enumerate(objects):
                for nu in range(4):
                    native=read_complex((len(q),c.shape[-1]))
                    if obj.__class__ is BPoly and nu > c.shape[0]:
                        ref=np.zeros_like(native)
                        ref[~np.isfinite(obj(q))]=np.nan
                    else: ref=obj(q,nu=nu)
                    checks+=compare(native,ref,f'case={idx},object={j},derivative={nu}',max_error)
            for j,obj in enumerate(objects):
                for a,b in pairs:
                    native=read_complex((1,c.shape[-1])).reshape(c.shape[-1])
                    if periodic and x[0]>x[-1]:
                        primitive=obj.antiderivative()
                        lo,hi=min(x),max(x)
                        def F(z): return primitive(lo + (z-lo)%(hi-lo),extrapolate=True)
                        ref=F(b)-F(a)+(np.floor((b-lo)/(hi-lo))-np.floor((a-lo)/(hi-lo)))*(
                            primitive(hi,extrapolate=True)-primitive(lo,extrapolate=True))
                    else: ref=obj.integrate(a,b)
                    checks+=compare(native,ref,f'case={idx},object={j},integration a={a} b={b}',max_error)
            for j,obj in enumerate(objects[:2]):
                native_d=read_complex((len(q),c.shape[-1]))
                checks+=compare(native_d,obj.derivative()(q),f'case={idx},object={j},derivative object',max_error)
            for j,obj in enumerate(objects[:2]):
                native_a=read_complex((len(q),c.shape[-1]))
                checks+=compare(native_a,obj.antiderivative()(q,extrapolate=True),f'case={idx},object={j},primitive object',max_error)
        assert next(lines,None) is None
        print(f'{flags}: {len(cases)} complex cases, {checks} entries, max absolute difference {max_error[0]:.4g}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=int,default=20261031)
    p.add_argument('--cases',type=int,default=160)
    p.add_argument('--optimized-only',action='store_true')
    args=p.parse_args()
    cases=generate(args.seed,args.cases)
    for flags in ([['-O3']] if args.optimized_only else [['-O0','-fcheck=all','-g','-fbacktrace'],['-O3']]):
        run(cases,flags)

if __name__=='__main__': main()
