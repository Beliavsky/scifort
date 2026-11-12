#!/usr/bin/env python3
"""Compare real/complex multichannel barycentric interpolation with SciPy.

Tests value and first three derivatives before and after incremental node
insertion. Also compiles the new native regression in checked and optimized
GNU Fortran builds. Independent randomized seeds can be supplied by --seed.
"""
from __future__ import annotations
import argparse
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from scipy.interpolate import BarycentricInterpolator

ROOT=Path(__file__).resolve().parents[2]
SOURCES=['src/kinds.f90','src/interpolate_piecewise.f90',
         'src/interpolate_bspline.f90','src/interpolate_regular_grid.f90',
         'src/interpolate_barycentric.f90', 'src/interpolate_global.f90','src/interpolate_kdtree.f90', 'src/interpolate_rbf.f90', 'src/interpolate_legacy.f90', 'src/interpolate_scattered.f90', 'src/interpolate_clough_tocher.f90', 'src/interpolate_fitpack_bivariate.f90', 'src/interpolate_fitpack_compat.f90', 'src/interpolate.f90']


def compile_fortran(directory,flags,driver):
    cmd=['gfortran','-std=f2018','-ffree-line-length-none',*flags,
         '-I'+str(directory),'-J'+str(directory),*(str(ROOT/s) for s in SOURCES),
         str(ROOT/driver),'-o',str(directory/'run')]
    p=subprocess.run(cmd,cwd=directory,text=True,capture_output=True)
    if p.returncode:
        raise RuntimeError(p.stderr[-5000:])
    return directory/'run'


def write_values(v):
    v=np.asarray(v,dtype=np.float64)
    return ' '.join(f'{a:.17g}' for a in v.ravel(order='F'))


def make_cases(seed,count):
    rng=np.random.default_rng(seed)
    cases=[]
    for k in range(count):
        n=2+k%8
        nch=1+(k//7)%4
        cmplx=k%4==2 or k%4==3
        nadd=k%4
        degree=min(n-1,1+(k//19)%6)
        # Uniform separation reduces avoidable ill-conditioning of high-order
        # differentiation matrices; still includes irregular and unsorted nodes.
        x=np.linspace(-0.8,0.8,n)+rng.uniform(-0.015,0.015,n)
        x=x[rng.permutation(n)]
        coeff=rng.uniform(-0.9,0.9,(degree+1,nch))
        if cmplx: coeff=coeff+1j*rng.uniform(-0.9,0.9,coeff.shape)
        calc=lambda xx:np.stack([np.polynomial.polynomial.polyval(xx,coeff[:,j]) for j in range(nch)],axis=-1)
        y=calc(x)
        q=np.r_[np.linspace(-0.95,0.95,17),x[:min(4,n)]]
        nx=rng.uniform(-0.85,0.85,nadd)
        for a in range(len(nx)):
            while np.any(abs(x-nx[a])<.04) or np.any(abs(nx[:a]-nx[a])<.04):
                nx[a]=rng.uniform(-0.85,0.85)
        ny=calc(nx) if nadd else np.zeros((0,nch),dtype=y.dtype)
        cases.append((cmplx,x,y,q,nx,ny))
    return cases


def serialize(cases):
    out=[str(len(cases))]
    for complex_data,x,y,q,nx,ny in cases:
        out+=[f'{2 if complex_data else 1} {len(x)} {y.shape[1]} {len(q)} {len(nx)}',
              write_values(x),write_values(y.real),write_values(y.imag),write_values(q)]
        if len(nx):
            out +=[write_values(nx),write_values(ny.real),write_values(ny.imag)]
    return '\n'.join(out)+'\n'


def compare(exe,cases):
    p=subprocess.run([str(exe)],input=serialize(cases),text=True,capture_output=True)
    if p.returncode: raise RuntimeError(p.stderr[-3500:])
    lines=iter(p.stdout.splitlines()); nvalues=0; maxdiff=0.; sci_add_bugs=0
    for index,(complex_data,x,y,q,nx,ny) in enumerate(cases,1):
        assert next(lines)==f'CASE {index}'
        # Use a deterministic RNG for SciPy's internally randomized weights.
        spline=BarycentricInterpolator(x,y,axis=0,rng=np.random.default_rng(119))
        for phase in (0,1):
            if phase==1 and len(nx):
                spline.add_xi(nx,ny)
                fresh=BarycentricInterpolator(np.r_[x,nx],np.concatenate([y,ny],axis=0),
                                              axis=0,rng=np.random.default_rng(119))
                if not np.allclose(spline(q),fresh(q),rtol=1e-9,atol=1e-9):
                    # SciPy 1.17 add_xi incorrectly signs the weights of a
                    # batch of even length. A freshly constructed polynomial
                    # is the independent reference for those cases.
                    sci_add_bugs+=1
                spline=fresh
            for nu in range(4):
                assert next(lines)==f'D {phase} {nu}'
                r=np.fromstring(next(lines),sep=' ')
                im=np.fromstring(next(lines),sep=' ')
                shape=(len(q),y.shape[1])
                assert r.size==np.prod(shape) and im.size==r.size,(index,nu,r.size)
                got=(r+1j*im if complex_data else r).reshape(shape,order='F')
                should=spline(q) if nu==0 else spline.derivative(q,der=nu)
                should=np.asarray(should).reshape(shape)
                # A numerical differentiation matrix can amplify floating-point
                # roundoff at larger orders. Polynomial operands are <=degree 6.
                atol=3e-6 if nu>=3 else (2e-7 if nu==2 else 1e-8)
                rtol=3e-6 if nu>=3 else (2e-7 if nu==2 else 1e-8)
                np.testing.assert_allclose(got,should,rtol=rtol,atol=atol,
                                           err_msg=f'case={index},phase={phase},nu={nu}')
                maxdiff=max(maxdiff,float(np.max(np.abs(got-should))))
                nvalues+=got.size
    assert next(lines,None) is None
    return nvalues,maxdiff,sci_add_bugs


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--seed',type=int,default=20261031)
    ap.add_argument('--cases',type=int,default=360)
    args=ap.parse_args()
    cases=make_cases(args.seed,args.cases)
    for kind,flags in [('checked',['-O0','-g','-fcheck=all','-fbacktrace']),
                       ('optimized',['-O3'])]:
        with tempfile.TemporaryDirectory() as tmp:
            work=Path(tmp)
            exe=compile_fortran(work,flags,'tools/interpolate_reference/barycentric_driver.f90')
            n,maxdiff,sci_add_bugs=compare(exe,cases)
            native=compile_fortran(work,flags,'test/test_interpolate_barycentric.f90')
            test=subprocess.run([str(native)],cwd=work,text=True,capture_output=True)
            if test.returncode:raise RuntimeError(test.stderr)
        print(f'{kind}: PASS {len(cases)} randomized cases, {n:,} complex/real values, '
              f'max absolute discrepancy {maxdiff:.3g}; native PASS; '+
              f'SciPy 1.17 even-batch add_xi discrepancies={sci_add_bugs}')


if __name__=='__main__':
    main()
