#!/usr/bin/env python3
"""Independent randomized SciPy comparison of the native 2-D C1 cubic interpolator.

Triangulation is native Fortran Bowyer-Watson; SciPy uses Qhull. Avoid repeated
sites and nearly co-circular/collinear configurations in randomized comparisons.
"""
import argparse
from pathlib import Path
import tempfile
import subprocess
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator
from scipy.spatial import Delaunay

ROOT=Path(__file__).resolve().parents[2]
SOURCES=('kinds','interpolate_kdtree','interpolate_scattered','interpolate_clough_tocher')


def create_cases(rng, nc):
    for k in range(nc):
        n=int(rng.integers(7, 65))
        channels=int(rng.integers(1, 5))
        rescale=bool(k%2)
        sites=rng.random((n,2))
        if k%5==0 and rescale:
            sites[:,0]*=2.0e-2
            sites[:,1]*=40.0
            sites+=np.array([20.,-16.])
        # Each construction is also tested at all input sites and at locations
        # strictly inside Qhull's triangles; outside-hull fill is tested too.
        triang=Delaunay(sites)
        ids=rng.integers(0,len(triang.simplices),size=18)
        weights=rng.dirichlet([2.,2.,2.],size=18)
        query=np.einsum('qt,qtd->qd', weights, sites[triang.simplices[ids]])
        query=np.concatenate([sites[:min(n,5)],query,(sites.max(axis=0)+np.ptp(sites,axis=0))[None,:]],axis=0)
        v=rng.normal(size=(n,channels))
        if k%4==0:
            v=np.column_stack([2+3*sites[:,0]-4*sites[:,1] for _ in range(channels)])
        elif k%4==1:
            v=np.sin(3.0*sites[:,0,None]+1.5*sites[:,1,None]+np.arange(channels)[None,:])
        if k%3==0:
            v=v+1j*rng.normal(size=(n,channels))
        else:
            v=v.astype(complex)
        yield sites, v, query, rescale


def run(cases,exe):
    lines=[]
    for sites,y,queries,rescale in cases:
        n,m=y.shape
        lines.append(f'{n} {m} {len(queries)} {int(rescale)}\n')
        lines += [f'{p[0]:.17e} {p[1]:.17e}\n' for p in sites]
        lines += [f'{z.real:.17e} {z.imag:.17e}\n' for row in y for z in row]
        lines += [f'{p[0]:.17e} {p[1]:.17e}\n' for p in queries]
    lines.append('0 0 0 0\n')
    completed=subprocess.run([str(exe)],input=''.join(lines),capture_output=True,text=True,timeout=170)
    if completed.returncode:
        raise RuntimeError(f'Fortran execution failed ({completed.returncode})\n'+completed.stderr[-3000:])
    output=iter(completed.stdout.splitlines())
    compared=0
    max_abs=max_rel=0.0
    differences=[]
    for i,(sites,y,queries,rescale) in enumerate(cases):
        code,converged=next(output).split()
        if code!='0' or converged!='T':
            raise RuntimeError(f'Fortran case {i} nonconvergent/error: code={code} converged={converged}')
        got=np.array([list(map(float,next(output).split())) for _ in range(len(queries)*y.shape[1])]).reshape(len(queries),y.shape[1],2)
        got=got[...,0]+1j*got[...,1]
        reference=CloughTocher2DInterpolator(sites,y,rescale=rescale,tol=1e-11,maxiter=1600)(queries)
        if np.any(np.isnan(got)!=np.isnan(reference)):
            raise AssertionError(f'case {i}: hull disagreement index={np.argwhere(np.isnan(got)!=np.isnan(reference)).tolist()}, got={got[np.isnan(got)!=np.isnan(reference)]}, reference={reference[np.isnan(got)!=np.isnan(reference)]}')
        mask=~np.isnan(reference)
        abs_err=np.abs(got[mask]-reference[mask])
        scale=1+np.abs(reference[mask])
        if len(abs_err):
            max_abs=max(max_abs,float(max(abs_err)))
            max_rel=max(max_rel,float(max(abs_err/scale)))
            if max(abs_err/scale)>1e-8:
                differences.append((i,float(max(abs_err/scale)),float(max(abs_err))))
        compared+=int(mask.sum())
    return compared,max_abs,max_rel,differences


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cases',type=int,default=80)
    ap.add_argument('--seed',type=int,default=20261031)
    ap.add_argument('--only',choices=['checked','optimized','both'],default='both')
    args=ap.parse_args()
    cases=list(create_cases(np.random.default_rng(args.seed),args.cases))
    modes=[]
    if args.only!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.only!='checked':modes.append(('optimized',['-O3']))
    for name,flags in modes:
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp)
            common=['gfortran','-std=f2018','-ffree-line-length-none',*flags,'-I'+str(d),'-J'+str(d)]
            for stem in SOURCES:
                p=subprocess.run([*common,'-c',str(ROOT/'src'/f'{stem}.f90'),'-o',str(d/f'{stem}.o')],cwd=d,capture_output=True,text=True)
                if p.returncode:raise RuntimeError(p.stderr)
            exe=d/'ct_test'
            p=subprocess.run([*common,str(ROOT/'tools/interpolate_reference/clough_tocher_driver.f90'),
                              *[str(d/f'{stem}.o') for stem in SOURCES],'-o',str(exe)],cwd=d,capture_output=True,text=True)
            if p.returncode:raise RuntimeError(p.stderr)
            compared,max_abs,max_rel,diff=run(cases,exe)
            print(f'{name}: {args.cases} scenarios; {compared} finite complex values; '
                  f'max abs={max_abs:.6g}; max normalized={max_rel:.6g}; above tolerance={len(diff)}',flush=True)
            if diff:raise AssertionError(f'{name}: largest mismatches: {sorted(diff,key=lambda a:a[1],reverse=True)[:10]}')

if __name__=='__main__':main()
