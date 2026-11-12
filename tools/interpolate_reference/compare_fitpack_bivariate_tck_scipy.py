#!/usr/bin/env python3
"""Independently compare SciFort's pure-Fortran bisplev against SciPy.

Does not require FITPACK, FPM, or an installed optional backend. Uses random
open tensor-product knot vectors and coefficients with all degrees 1 through 5.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import bisplev

ROOT=Path(__file__).resolve().parents[2]
SOURCES=('kinds','interpolate_bspline','interpolate_fitpack_bivariate')


def run(args,cwd,source=None):
    p=subprocess.run(list(map(str,args)),cwd=cwd,input=source,text=True,capture_output=True)
    if p.returncode:
        raise RuntimeError(f'{args}:\n{p.stdout[-1500:]}\n{p.stderr[-3000:]}')
    return p.stdout


def build(work,flags):
    common=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    for name in SOURCES:
        run([*common,'-I'+str(work),'-J'+str(work),'-c',ROOT/'src'/f'{name}.f90',
             '-o',work/f'{name}.o'],work)
    exe=work/'fitpack_bivariate_tck_driver'
    run([*common,'-I'+str(work),ROOT/'tools/interpolate_reference/fitpack_bivariate_tck_driver.f90',
         *(work/f'{name}.o' for name in SOURCES),'-o',exe],work)
    return exe


def check(exe,rng,i):
    kx,ky=1+i%5,1+(i//5)%5
    lo,hi=-1.2,1.7
    nintx=int(rng.integers(0,6)); ninty=int(rng.integers(0,6))
    tx=np.r_[np.full(kx+1,lo),np.sort(rng.uniform(lo+.05,hi-.05,nintx)),
             np.full(kx+1,hi)]
    ty=np.r_[np.full(ky+1,-.8),np.sort(rng.uniform(-.7,1.9,ninty)),
             np.full(ky+1,2.)]
    nx=len(tx)-kx-1;ny=len(ty)-ky-1
    c=rng.standard_normal(nx*ny)
    qx=rng.uniform(lo-.6,hi+.6,25)
    qy=rng.uniform(-1.2,2.4,25)
    dx=rng.integers(0,kx,len(qx));dy=rng.integers(0,ky,len(qy))
    lines=[f'{len(tx)} {len(ty)} {kx} {ky}']
    lines += [' '.join(f'{v:.17g}' for v in arr) for arr in (tx,ty,c)]
    lines += [str(len(qx))]
    lines += [f'{a:.17g} {b:.17g} {j} {k}' for a,b,j,k in zip(qx,qy,dx,dy)]
    result=run([exe],exe.parent,'\n'.join(lines)+'\n').strip().splitlines()
    if len(result)!=len(qx):raise AssertionError((i,'output length'))
    got=[]
    for line in result:
        a,b=line.split();got.append(float(a))
        if int(b):raise AssertionError((i,'Fortran error code',b))
    tck=(tx,ty,c,kx,ky)
    scipy_values=np.array([bisplev(float(x),float(y),tck,dx=int(j),dy=int(k))
                           for x,y,j,k in zip(qx,qy,dx,dy)])
    error=float(np.max(abs(scipy_values-np.asarray(got))))
    np.testing.assert_allclose(got,scipy_values,atol=1.e-9,rtol=5.e-12,
                               err_msg=f'case={i}, degree={(kx,ky)}')
    return error


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cases',type=int,default=200)
    ap.add_argument('--seed',type=int,default=149173)
    ap.add_argument('--only',choices=('both','checked','optimized'),default='both')
    args=ap.parse_args()
    modes=[]
    if args.only!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all']))
    if args.only!='checked':modes.append(('optimized',['-O3']))
    for mode,flags in modes:
        with tempfile.TemporaryDirectory() as directory:
            exe=build(Path(directory),flags)
            rng=np.random.default_rng(args.seed)
            error=max(check(exe,rng,i) for i in range(args.cases))
            print(f'{mode}: {args.cases}/{args.cases} FITPACK TCK cases PASS;'
                  f' max absolute error {error:.4g}',flush=True)

if __name__=='__main__':main()
