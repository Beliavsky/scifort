#!/usr/bin/env python3
"""Compare actual optional FITPACK 2.0.0 bivariate fitting against SciPy.

The Fortran executable is built by check_optional_fitpack_local.py or locally:
  python tools/interpolate_reference/compare_optional_fitpack_bivariate_scipy.py \
      --executable /path/to/fitpack_bivariate_driver --cases 100

Tests independent SciPy bisplev evaluations and derivatives of the returned
TCK, the FITPACK reported weighted residual, and SciPy LSQBivariateSpline
on identical preassigned interior knots. Automatic knot selection can differ
between FITPACK 2.0.0 and SciPy, and is not required to match exactly.
"""
from __future__ import annotations
import argparse
from collections import Counter
from pathlib import Path
import subprocess
import warnings
import numpy as np
from scipy.interpolate import bisplev, LSQBivariateSpline, bisplrep


def call_driver(exe: Path, mode: int, x: np.ndarray, y: np.ndarray,
                z: np.ndarray, w: np.ndarray, kx: int, ky: int,
                s: float, tx: np.ndarray, ty: np.ndarray,
                qx: np.ndarray, qy: np.ndarray):
    lines = [f'{mode} {len(x)} {kx} {ky} {s:.17g} {len(tx)} {len(ty)}']
    lines += [f'{a:.17g} {b:.17g} {c:.17g} {d:.17g}'
              for a,b,c,d in zip(x,y,z,w)]
    if len(tx):lines.append(' '.join(f'{v:.17g}' for v in tx))
    if len(ty):lines.append(' '.join(f'{v:.17g}' for v in ty))
    lines += [str(len(qx))]
    lines += [f'{a:.17g} {b:.17g}' for a,b in zip(qx,qy)]
    process = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                             capture_output=True, text=True, timeout=30)
    if process.returncode:
        raise RuntimeError(f'FITPACK driver exited {process.returncode}: {process.stderr[-1500:]}')
    response = process.stdout.strip().splitlines()
    code,resid = response[0].split()
    if int(code):
        if len(response)!=1: raise AssertionError('Invalid adapter response')
        return int(code), float(resid), None
    a,b,kxx,kyy = map(int,response[1].split())
    if (kxx,kyy)!=(kx,ky): raise AssertionError('Degrees changed unexpectedly')
    data = np.array([float(line.strip()) for line in response[2:]])
    ncoef=(a-kx-1)*(b-ky-1)
    if len(data)!=(a+b+ncoef+len(qx)*(2 if min(kx,ky)>1 else 1)):
        raise AssertionError(f'Invalid TCK output length: {len(data)}')
    tck = (data[:a],data[a:a+b],data[a+b:a+b+ncoef],kx,ky)
    evaluated=data[a+b+ncoef:a+b+ncoef+len(qx)]
    mixed=data[a+b+ncoef+len(qx):]
    return int(code),float(resid),(tck,evaluated,mixed)


def test_case(exe:Path, rng:np.random.Generator, i:int, counters:Counter):
    mode=i%2
    kx=1+(i//2)%5
    ky=kx if mode==0 else 1+(i//10)%5
    n=65+(i*17)%24
    x=rng.uniform(-1.5,1.0,n)
    y=rng.uniform(-0.5,1.7,n)
    w=rng.uniform(0.7,1.6,n)
    z=np.sin(3.3*x)+np.cos(2.6*y)+0.17*x*y+0.02*rng.standard_normal(n)
    s=float(0.12+0.07*(i%5))
    if mode:
        tx=np.array([np.min(x)+0.47*np.ptp(x)])
        ty=np.array([np.min(y)+0.54*np.ptp(y)])
    else:
        tx=ty=np.empty(0)
    qx=np.r_[rng.uniform(x.min(),x.max(),20),x.min()-.2,x.max()+.2]
    qy=np.r_[rng.uniform(y.min(),y.max(),20),y.min()-.25,y.max()+.25]
    code,resid,response=call_driver(exe,mode,x,y,z,w,kx,ky,s,tx,ty,qx,qy)
    if code:
        counters['upstream_failures']+=1
        return 0.,0.
    counters['success']+=1
    if mode:counters['least_squares']+=1
    else:counters['automatic']+=1
    tck,got,mixed=response
    eval_ref=np.array([bisplev(float(a),float(b),tck) for a,b in zip(qx,qy)])
    max_eval=float(np.max(np.abs(eval_ref-got)))
    np.testing.assert_allclose(got,eval_ref,atol=3.e-10,rtol=1.e-10)
    if min(kx,ky)>1:
        dref=np.array([bisplev(float(a),float(b),tck,dx=1,dy=1) for a,b in zip(qx,qy)])
        np.testing.assert_allclose(mixed,dref,atol=2.e-8,rtol=2.e-10)
        max_eval=max(max_eval,float(np.max(np.abs(mixed-dref))))
    # Independent residual reconstruction using SciPy's bisplev on every site
    fit=np.array([bisplev(float(a),float(b),tck) for a,b in zip(x,y)])
    recomputed=float(np.sum((w*(z-fit))**2))
    err_resid=abs(resid-recomputed)
    if err_resid > 2.e-7*(1+recomputed):
        raise AssertionError(f'Incorrect weighted residual: {resid}, {recomputed}')
    if mode:
        scipy_fit=LSQBivariateSpline(x,y,z,tx,ty,w=w,kx=kx,ky=ky)
        scipy_eval=scipy_fit.ev(qx,qy)
        # The two core FITPACK implementations should agree when the knots are fixed.
        np.testing.assert_allclose(got,scipy_eval,atol=3.e-8,rtol=2.e-8)
        counters['scipy_fixed_knot_matches']+=1
    else:
        if resid > s*1.003+1.e-9:
            raise AssertionError(f'Incorrect FITPACK residual target: {resid} > {s}')
        # Diagnostic only: upstream and SciPy may select different knots.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                sci=bisplrep(x,y,z,w=w,kx=kx,ky=ky,s=s,quiet=1)
            if len(sci[0])==len(tck[0]) and len(sci[1])==len(tck[1]) \
               and np.allclose(sci[0],tck[0],atol=1.e-10) \
               and np.allclose(sci[1],tck[1],atol=1.e-10):
                counters['scipy_smoothing_knots_match']+=1
        except (ValueError,RuntimeError,OverflowError):
            counters['scipy_smoothing_reference_unavailable']+=1
    return max_eval,err_resid


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',required=True,type=Path)
    parser.add_argument('--cases',type=int,default=100)
    parser.add_argument('--seed',type=int,default=15231)
    args=parser.parse_args()
    args.executable = args.executable.resolve()
    rng=np.random.default_rng(args.seed)
    counts:Counter=Counter()
    max_diff=max_resid=0.
    for i in range(args.cases):
        ed,rd=test_case(args.executable,rng,i,counts)
        max_diff=max(max_diff,ed)
        max_resid=max(max_resid,rd)
    if counts['success'] < int(args.cases*.85):
        raise AssertionError(f'Too many upstream failures: {dict(counts)}')
    if counts['least_squares'] < int(args.cases*.4):
        raise AssertionError(f'Insufficient fixed-knot cases: {dict(counts)}')
    print(f'{counts["success"]}/{args.cases} successful bivariate FITPACK cases;'
          f' fixed knots {counts["least_squares"]}, adaptive {counts["automatic"]};'
          f' failures={counts["upstream_failures"]};'
          f' maximum TCK evaluation difference={max_diff:.4g};'
          f' maximum residual difference={max_resid:.4g};'
          f' smoothing knot matches={counts["scipy_smoothing_knots_match"]}',flush=True)

if __name__=='__main__': main()
