#!/usr/bin/env python3
"""Independently compare SciFort's RBFInterpolator with scipy.interpolate.RBFInterpolator.

Native GNU Fortran compilation is done once per mode. This test includes global and
local fits, all named SciPy-1.17 kernels, complex data, multichannel values,
per-site smoothing, and 1D-4D scattered sites. Does not require FPM or BLAS.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import RBFInterpolator

ROOT=Path(__file__).resolve().parents[2]
KERNELS=["linear","thin_plate_spline","cubic","quintic", "multiquadric",
         "inverse_multiquadric","inverse_quadratic","gaussian"]
MIND={'linear':0,'thin_plate_spline':1,'cubic':1,'quintic':2,'multiquadric':0}


def run(cmd, cwd):
    proc=subprocess.run(cmd,cwd=cwd,text=True,capture_output=True)
    if proc.returncode:
        raise RuntimeError(f"{' '.join(map(str,cmd))}\n{proc.stderr[-2500:]}\n{proc.stdout[-1000:]}")
    return proc.stdout


def compile_driver(work, flags):
    common=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    for name in ['kinds','interpolate_kdtree','interpolate_rbf']:
        run([*common,'-I'+str(work),'-J'+str(work),'-c',str(ROOT/'src'/f'{name}.f90'),
             '-o',str(work/f'{name}.o')],work)
    exe=work/'rbf_driver'
    run([*common,'-I'+str(work),str(ROOT/'tools/interpolate_reference/rbf_driver.f90'),
         str(work/'interpolate_rbf.o'),str(work/'interpolate_kdtree.o'),str(work/'kinds.o'),'-o',str(exe)],work)
    return exe


def reference_case(rng, idx):
    kernel=KERNELS[idx%len(KERNELS)]
    dim=1+(idx//len(KERNELS))%4
    n=18+idx%21
    channels=1+idx%4
    degree=MIND.get(kernel,0)
    if idx%9==0:
        degree=1 if kernel not in ('quintic',) else 2
    if idx%11==0:
        degree=-1
    if degree==2 and dim==4:
        n=max(n,24)
    x=rng.normal(size=(n,dim))*rng.uniform(0.5,2.0,size=(1,dim))
    z=rng.normal(size=(n,channels)) + 1j*rng.normal(size=(n,channels))
    if idx%7==0:
        z=np.real(z).astype(complex)
    q=rng.uniform(-1.4,1.4,size=(7,dim))
    eps=float(rng.uniform(0.8,1.7))
    lam=float(rng.uniform(1e-4,0.05))
    if idx%4==0 and kernel in ('linear','thin_plate_spline','cubic','quintic'):
        lam=0.0
    smoothing=np.full(n,lam)
    if idx%5==0: smoothing*=rng.uniform(0.1,1.9,size=n)
    neighbors=0
    if idx%3==0:
        neighbors=max((dim+1)*(dim+2)//2+2, n//2)
        neighbors=min(n,neighbors)
    return kernel,dim,n,channels,degree,eps,x,z,q,smoothing,neighbors


def test_one(exe, case, idx):
    kernel,dim,n,channels,degree,eps,x,z,q,smoothing,neighbors=case
    rows=[f'{n} {dim} {channels} {len(q)} {degree} {neighbors}',f'{kernel} {eps:.17g}']
    rows.extend(' '.join([f'{smoothing[i]:.17g}',*[f'{t:.17g}' for t in x[i]]]) for i in range(n))
    rows.extend(f'{v.real:.17g} {v.imag:.17g}' for row in z for v in row)
    rows.extend(' '.join(f'{t:.17g}' for t in row) for row in q)
    proc=subprocess.run([str(exe)],input='\n'.join(rows)+'\n',text=True,capture_output=True)
    if proc.returncode:
        raise AssertionError(f'case {idx}: native crash: {proc.stderr[-1000:]}')
    lines=proc.stdout.splitlines()
    if len(lines)<2 or lines[0].strip()!='0' or lines[1].strip()!='0':
        raise AssertionError(f'case {idx}: invalid status: {proc.stdout[:200]}')
    native=np.array([list(map(float,row.split())) for row in lines[2:]])
    expected_shape=(q.shape[0]*channels,2)
    if native.shape!=expected_shape:
        raise AssertionError(f'case {idx}: got {native.shape} instead of {expected_shape}')
    native=(native[:,0]+1j*native[:,1]).reshape((len(q),channels))
    reference=RBFInterpolator(x,z,neighbors=neighbors if neighbors else None,
        smoothing=smoothing,kernel=kernel,epsilon=eps,degree=degree)(q)
    err=np.max(np.abs(native-reference))
    scaled=np.max(np.abs(native-reference)/(1+np.abs(reference)))
    # Very small regularization and almost coinciding points can magnify the
    # condition number; only well-conditioned random sites are used here.
    tol=1e-7
    if scaled>tol:
        raise AssertionError(f'case {idx} ({kernel} dim={dim} degree={degree} '
            f'neighbors={neighbors}): error={err:.3g}, scaled={scaled:.3g}\n'
            f'native={native[:3]}, expected={reference[:3]}')
    return err,scaled,native.size


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--seed',type=int,default=20261031)
    ap.add_argument('--cases',type=int,default=220)
    ap.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    args=ap.parse_args()
    modes=[]
    if args.mode!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.mode!='checked':modes.append(('optimized',['-O3']))
    for title,flags in modes:
        rng=np.random.default_rng(args.seed)
        with tempfile.TemporaryDirectory() as tmp:
            exe=compile_driver(Path(tmp),flags)
            max_abs=max_scaled=0.0
            entries=0
            for idx in range(args.cases):
                case=reference_case(rng,idx)
                try:
                    a,s,n=test_one(exe,case,idx)
                except Exception as exc:
                    raise RuntimeError(f'{title} seed={args.seed}, index={idx}') from exc
                max_abs=max(max_abs,a)
                max_scaled=max(max_scaled,s)
                entries+=n
            print(f'{title}: PASS {args.cases} cases, {entries} complex entries; '
                  f'max abs {max_abs:.5g}; max scaled {max_scaled:.5g}',flush=True)

if __name__=='__main__':main()
