#!/usr/bin/env python3
"""Independent SciPy comparisons of Krogh and Floater-Hormann interpolation.

Compiles checked and optimized Fortran builds and compares randomized real and
complex, scalar/multichannel interpolators. The inputs are also checked against
known Hermite polynomial constraints and exact Floater nodal interpolation.
"""
from __future__ import annotations
import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import FloaterHormannInterpolator, KroghInterpolator

ROOT=Path(__file__).resolve().parents[2]
SOURCES=['src/kinds.f90','src/interpolate_piecewise.f90','src/interpolate_bspline.f90',
         'src/interpolate_regular_grid.f90','src/interpolate_barycentric.f90',
         'src/interpolate_global.f90','src/interpolate_kdtree.f90', 'src/interpolate_rbf.f90', 'src/interpolate_legacy.f90', 'src/interpolate_scattered.f90', 'src/interpolate_clough_tocher.f90', 'src/interpolate_fitpack_bivariate.f90', 'src/interpolate_curfit_core.f90', 'src/interpolate_fitpack_compat.f90', 'src/interpolate.f90']


def compile_fortran(work,flags,driver):
    cmd=['gfortran','-std=f2018','-ffree-line-length-none',*flags,
         '-J'+str(work),'-I'+str(work),*(str(ROOT/s) for s in SOURCES),
         str(ROOT/driver),'-o',str(work/'run')]
    p=subprocess.run(cmd,cwd=work,text=True,capture_output=True)
    if p.returncode:raise RuntimeError('Fortran compilation failed:\n'+p.stderr[-5000:])
    return work/'run'


def make_cases(seed,count):
    rng=np.random.default_rng(seed)
    cases=[]
    for i in range(count):
        algo=1 if i%2==0 else 2
        is_complex=(i//2)%2==1
        nchan=1+(i//5)%4
        if algo==1:
            ng=2+(i//3)%5
            ng=min(ng,7)
            nodes=np.linspace(-.85,.85,ng)+rng.uniform(-.018,.018,ng)
            repeat=np.array([1+(i+j*3//2)%3 for j in range(ng)])
            n=int(repeat.sum())
            deg=min(n-1,4+(i//29)%3)
            coeff=rng.uniform(-.5,.5,(deg+1,nchan))
            if is_complex:coeff=coeff+1j*rng.uniform(-.5,.5,coeff.shape)
            x=np.repeat(nodes,repeat)
            y=np.empty((n,nchan),dtype=np.complex128 if is_complex else np.float64)
            j=0
            for node,r in zip(nodes,repeat):
                for order in range(r):
                    dc=np.polynomial.polynomial.polyder(coeff,m=order,axis=0)
                    y[j,:]=[np.polynomial.polynomial.polyval(node,dc[:,c]) for c in range(nchan)]
                    j+=1
            d=0
            order=3
            q=np.r_[np.linspace(-1.05,1.05,19),nodes]
            cases.append((algo,is_complex,x,y,q,order,d,coeff))
        else:
            n=4+(i//7)%18
            x=np.linspace(-1.,1.,n)+rng.uniform(-.025,.025,n)
            # Both increasing and decreasing grids, when supported by SciPy.
            if i%8==1: x=x[::-1].copy()
            d=[0,1,2,3,4,6,n-1][(i//9)%7]
            d=min(d,n-1)
            if d==n-1 and n>12: d=6  # avoid ill-conditioned high-degree polynomial extrapolation
            y=rng.standard_normal((n,nchan))*.4 + np.sin(3*x[:,None])*rng.uniform(-1,1,(1,nchan))
            if is_complex:y=y+1j*(rng.standard_normal((n,nchan))*.4 + .6*np.cos(2*x[:,None]))
            order=0
            q=np.r_[np.linspace(-1.15,1.15,23),x[:min(n,4)]]
            cases.append((algo,is_complex,x,y,q,order,d,None))
    return cases


def serialize(cases):
    out=[str(len(cases))]
    def line(x):return ' '.join(f'{v:.17g}' for v in np.asarray(x).ravel(order='F'))
    for algo,cmplx,x,y,q,order,d,coeff in cases:
        out += [f'{algo} {2 if cmplx else 1} {len(x)} {y.shape[1]} {len(q)} {order} {d}',
                line(x),line(y.real),line(y.imag),line(q)]
    return '\n'.join(out)+'\n'


def compare(exe,cases):
    p=subprocess.run([str(exe)],input=serialize(cases),text=True,capture_output=True)
    if p.returncode:raise RuntimeError('Fortran execution failed: '+p.stderr[-4000:])
    lines=iter(p.stdout.splitlines())
    max_abs=0.;max_scaled=0.;nvals=0
    for index,(algo,cmplx,x,y,q,order,d,coeff) in enumerate(cases,1):
        assert next(lines)==f'CASE {index}',(index,'missing case marker')
        if algo==1:
            ref=KroghInterpolator(x,y,axis=0)
        else:
            ref=FloaterHormannInterpolator(x,y,d=d,axis=0)
        for nu in range(order+1):
            assert next(lines)==f'ORDER {nu}',(index,nu)
            real=np.fromstring(next(lines),sep=' ')
            imag=np.fromstring(next(lines),sep=' ')
            shape=(len(q),y.shape[1])
            if real.size!=int(np.prod(shape)) or imag.size!=real.size:
                raise AssertionError((index,nu,shape,real.size,imag.size))
            got=(real+1j*imag).reshape(shape,order='F')
            wanted=(ref(q) if algo==2 or nu==0 else ref.derivative(q,der=nu))
            wanted=np.asarray(wanted).reshape(shape)
            # Krogh's higher-order Newton divided differences are sensitive to
            # node spacing; judge relative to the magnitude of the result.
            atol=2e-8 if algo==1 and nu>=2 else 2e-10
            rtol=2e-8 if algo==1 and nu>=2 else 4e-10
            np.testing.assert_allclose(got,wanted,rtol=rtol,atol=atol,
                                       err_msg=f'case={index}, algo={algo}, derivative={nu}')
            if algo==1 and coeff is not None:
                dc=np.polynomial.polynomial.polyder(coeff,m=nu,axis=0)
                exact=np.stack([np.polynomial.polynomial.polyval(q,dc[:,k])
                                for k in range(y.shape[1])],axis=1)
                np.testing.assert_allclose(got,exact,rtol=2e-8,atol=2e-8,
                                           err_msg=f'polynomial constraint {index}/{nu}')
            scale=1+np.abs(wanted)
            max_abs=max(max_abs,float(np.max(abs(got-wanted))))
            max_scaled=max(max_scaled,float(np.max(abs(got-wanted)/scale)))
            nvals+=int(np.prod(shape))
        if algo==1:
            nord=min(len(x),4)
            assert next(lines)==f'ALL {nord}',(index,'derivative array marker')
            rr=np.fromstring(next(lines),sep=' ')
            ii=np.fromstring(next(lines),sep=' ')
            shape=(nord,len(q),y.shape[1])
            assert rr.size==int(np.prod(shape)) and ii.size==rr.size
            derivatives=(rr+1j*ii).reshape(shape,order='F')
            ref_all=np.asarray(ref.derivatives(q,der=nord)).reshape(shape)
            np.testing.assert_allclose(derivatives,ref_all,rtol=2e-8,atol=2e-8,
                                       err_msg=f'all derivative orders, case {index}')
            max_abs=max(max_abs,float(np.max(abs(derivatives-ref_all))))
            max_scaled=max(max_scaled,float(np.max(abs(derivatives-ref_all)/(1+abs(ref_all)))))
            nvals+=int(np.prod(shape))
    if next(lines,None) is not None:raise AssertionError('unexpected output lines')
    return nvals,max_abs,max_scaled


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--seed',type=int,default=20261031)
    ap.add_argument('--cases',type=int,default=240)
    args=ap.parse_args()
    cases=make_cases(args.seed,args.cases)
    for build,flags in [('checked',['-O0','-g','-fcheck=all','-fbacktrace']),
                        ('optimized',['-O3'])]:
        with tempfile.TemporaryDirectory() as tmp:
            work=Path(tmp)
            exe=compile_fortran(work,flags,'tools/interpolate_reference/global_interpolators_driver.f90')
            nvals,max_abs,max_scaled=compare(exe,cases)
            native=compile_fortran(work,flags,'test/test_interpolate_global.f90')
            result=subprocess.run([str(native)],text=True,capture_output=True)
            if result.returncode:raise RuntimeError(result.stderr)
        print(f'{build}: PASS {args.cases} SciPy cases, {nvals:,} outputs, '
              f'max absolute error {max_abs:.3g}, scaled error {max_scaled:.3g}; native PASS')

if __name__=='__main__':main()
