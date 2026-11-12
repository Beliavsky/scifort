#!/usr/bin/env python3
"""Independent SciPy comparisons for SciFort's legacy interp1d implementation.

Includes exact knot and nearest-midpoint tie behavior, unsorted sites,
complex and multichannel ordinates, all supported kinds, asymmetric fill,
and extrapolation; also checks real-valued and scalar public APIs natively.
"""
from __future__ import annotations
import argparse
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d

ROOT=Path(__file__).resolve().parents[2]
KINDS=['linear','nearest','nearest-up','previous','next','zero','slinear','quadratic','cubic']


def run(cmd,where):
    p=subprocess.run(cmd,cwd=where,text=True,capture_output=True)
    if p.returncode: raise RuntimeError(' '.join(map(str,cmd))+'\n'+p.stderr[-1500:])
    return p.stdout


def build(w,flags):
    common=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    for name in ('kinds','interpolate_bspline','interpolate_legacy'):
        run([*common,'-I'+str(w),'-J'+str(w),'-c',str(ROOT/'src'/f'{name}.f90'),
             '-o',str(w/f'{name}.o')],w)
    exe=w/'interp1d_driver'
    run([*common,'-I'+str(w),str(ROOT/'tools/interpolate_reference/interp1d_driver.f90'),
         *(str(w/f'{n}.o') for n in ('interpolate_legacy','interpolate_bspline','kinds')),
         '-o',str(exe)],w)
    return exe


def check(exe,rng,index):
    kind=KINDS[index%len(KINDS)]
    n=4+index%19
    channels=1+(index//5)%4
    x=np.cumsum(rng.uniform(0.14,1.1,n))-2.0
    ys=rng.normal(size=(n,channels))
    z=ys+1j*rng.normal(size=(n,channels))
    if index%4==0:z=ys.astype(complex)
    q=np.concatenate((np.linspace(x[0]-0.8,x[-1]+0.8,14),x.copy(),
                     (x[1:]+x[:-1])/2))
    # Unsorted input is deliberately tested.
    if index%5==0:
        order=rng.permutation(n)
        x=x[order]
        z=z[order]
    extrapolate=bool(index%3==0)
    low=complex(-99.,17.)
    high=complex(123.,-19.)
    rows=[f'{n} {channels} {len(q)} {int(extrapolate)}',
          f'{kind} {low.real:.17g} {low.imag:.17g} {high.real:.17g} {high.imag:.17g}',
          ' '.join(f'{v:.17g}' for v in x)]
    rows.extend(f'{v.real:.17g} {v.imag:.17g}' for a in z for v in a)
    rows.append(' '.join(f'{v:.17g}' for v in q))
    p=subprocess.run([str(exe)],input='\n'.join(rows)+'\n',text=True,capture_output=True)
    if p.returncode: raise AssertionError(f'case {index}: Fortran failure: {p.stderr}')
    lines=p.stdout.splitlines()
    if len(lines)!=2+len(q)*channels or lines[0]!='0' or lines[1]!='0':
        raise AssertionError(f'case {index}: status or line count: {p.stdout[:150]}')
    values=np.array([list(map(float,s.split())) for s in lines[2:]])
    got=(values[:,0]+1j*values[:,1]).reshape(len(q),channels)
    f=interp1d(x,z,kind=kind,axis=0,assume_sorted=False,bounds_error=False,
               fill_value='extrapolate' if extrapolate else (low,high))
    expected=f(q)
    err=np.nanmax(np.abs(got-expected))
    scaled=np.nanmax(np.abs(got-expected)/(1+np.abs(expected)))
    if not np.allclose(got,expected,rtol=3e-10,atol=2e-10,equal_nan=True):
        at=np.unravel_index(np.nanargmax(np.abs(got-expected)),got.shape)
        raise AssertionError(f'case {index} kind {kind} extrapolate={extrapolate}, '
            f'diff={err:.4g}, at={at}, q={q[at[0]]}, sorted_x={np.sort(x)}, got={got[at]}, wanted={expected[at]}')
    return float(err),float(scaled),got.size


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cases',type=int,default=360)
    ap.add_argument('--seed',type=int,default=20261031)
    ap.add_argument('--mode',choices=['both','checked','optimized'],default='both')
    args=ap.parse_args()
    modes=[]
    if args.mode!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.mode!='checked':modes.append(('optimized',['-O3']))
    for name,flags in modes:
        rng=np.random.default_rng(args.seed)
        with tempfile.TemporaryDirectory() as tmp:
            exe=build(Path(tmp),flags)
            maxabs=maxscaled=0.0
            entries=0
            for i in range(args.cases):
                try:a,s,num=check(exe,rng,i)
                except Exception as exc:raise RuntimeError(f'{name}, seed={args.seed}, index={i}') from exc
                entries+=num
                maxabs=max(maxabs,a)
                maxscaled=max(maxscaled,s)
            print(f'{name}: PASS {args.cases} scenarios, {entries} complex entries; '
                  f'max absolute error {maxabs:.5g}; max scaled error {maxscaled:.5g}',flush=True)

if __name__=='__main__':main()
