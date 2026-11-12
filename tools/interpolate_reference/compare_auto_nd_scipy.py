#!/usr/bin/env python3
"""Compare native 4D/5D Delaunay and indexed LinearND against SciPy 1.17.

The native geometry and SciPy/Qhull simplex sets must match in nondegenerate
cases. On extremely ill-conditioned inputs Qhull can return NaN even at a known
convex combination of simplex vertices; these cases are checked against their
independent barycentric construction and counted separately.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

ROOT=Path(__file__).resolve().parents[2]
SRC=('kinds','interpolate_kdtree','interpolate_scattered')


def run(cmd, *,cwd,stdin=None):
    p=subprocess.run(cmd,cwd=cwd,input=stdin,text=True,capture_output=True)
    if p.returncode:
        raise RuntimeError(f"{cmd}: {p.stderr[-2000:]}\n{p.stdout[-1000:]}")
    return p.stdout


def case_data(rng,dim,index):
    n=dim+3+(index*7)%(24 if dim==5 else 30)
    points=rng.uniform(-2,2,(n,dim))
    mode=('normal','thin','offset','cluster','wide','shell')[index%6]
    if mode=='thin': points[:,-1]*=1e-3
    if mode=='offset': points+=rng.uniform(-2e5,2e5,size=dim)
    if mode=='cluster': points[:n//2]*=0.01
    if mode=='wide': points*=10**rng.uniform(-2,2,size=dim)
    if mode=='shell':
        points/=np.linalg.norm(points,axis=1)[:,None]
        points*=1+0.12*rng.normal(size=(n,1))
    return points,mode


def check_case(exe,rng,dim,index,work):
    points,mode=case_data(rng,dim,index)
    n=len(points)
    # Public facade must also work when the axes are individually rescaled.
    rescale=index%2==0
    channels=1+index%3
    values=rng.standard_normal((n,channels))+1j*rng.standard_normal((n,channels))
    center=np.mean(points,axis=0)
    tri_geom=Delaunay(points-center)
    axes=(points-np.min(points,axis=0))/np.ptp(points,axis=0) if rescale else points-center
    tri_interp=Delaunay(axes)
    n_inside=min(12,len(tri_interp.simplices))
    chosen=tri_interp.simplices[rng.integers(len(tri_interp.simplices),size=n_inside)]
    weights=rng.dirichlet(0.9*np.ones(dim+1),size=n_inside)
    q=np.einsum('nj,njk->nk',weights,points[chosen])
    test_nodes=[0,n//2]
    q=np.vstack([q,points[test_nodes],points.mean(axis=0),points.max(axis=0)+2])
    lines=[f'{n} {dim} {channels} {len(q)} {int(rescale)}']
    lines+=[' '.join(f'{v:.17g}' for v in row) for row in points]
    lines+=[' '.join(f'{v:.17g}' for v in [*row.real,*row.imag]) for row in values]
    lines+=[' '.join(f'{v:.17g}' for v in row) for row in q]
    out=run([str(exe)],cwd=work,stdin='\n'.join(lines)+'\n').splitlines()
    status,ncell=map(int,out[0].split())
    if status:raise AssertionError(f'dim={dim} case={index} native triangulation status={status}')
    cells=np.array([list(map(int,s.split())) for s in out[1:ncell+1]])-1
    native_set={tuple(sorted(c)) for c in cells}
    ref_set={tuple(sorted(c)) for c in tri_geom.simplices}
    if native_set!=ref_set:
        raise AssertionError(f'dim={dim} case={index} {mode} simplex set mismatch; '
                             f'native={len(native_set)} Qhull={len(ref_set)} '
                             f'extra={len(native_set-ref_set)} missing={len(ref_set-native_set)}')
    if [int(s) for s in out[1+ncell:3+ncell]] != [0,0]:
        raise AssertionError(f'public constructor/evaluate error: {out[1+ncell:3+ncell]}')
    observed=np.array([[float(v) for v in s.split()] for s in out[3+ncell:]])
    if observed.shape!=(len(q),2*channels):
        raise AssertionError(f'unexpected output shape {observed.shape}')
    got=observed[:,:channels]+1j*observed[:,channels:]
    expected=LinearNDInterpolator(points,values,fill_value=-99+7j,rescale=rescale)(q)
    if expected.ndim==1:expected=expected[:,None]
    # Qhull occasionally marks a known interior query as outside when its
    # barycentric system is ill conditioned. In that case use its defining
    # convex combination as the reference, rather than accepting a NaN.
    fallback=0
    for j in range(n_inside):
        if not np.all(np.isfinite(expected[j,:])):
            expected[j,:]=np.sum(weights[j,:,None]*values[chosen[j],:],axis=0)
            fallback+=1
    for j,node in enumerate(test_nodes,start=n_inside):
        if not np.all(np.isfinite(expected[j,:])):
            expected[j,:]=values[node,:]
            fallback+=1
    # A mean of all nodes is always in their convex hull.  The independent
    # reference for this one query is provided only when SciPy returns finite.
    mask=np.isfinite(expected)
    if np.any(mask & ~np.isfinite(got)):
        raise AssertionError(f'dim={dim} {mode} unexpected fill at finite SciPy interior query')
    bad=np.abs(got-expected)[mask]
    err=float(np.max(bad,initial=0))
    if err>2e-8:
        raise AssertionError(f'dim={dim} case={index} {mode}, max error={err:.3g}')
    if not np.all(got[-1,:]==(-99+7j)):
        raise AssertionError(f'dim={dim} case={index} out-of-hull fill wrong')
    return len(ref_set),int(np.sum(mask))*2,err,fallback


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--seed',type=int,default=194105)
    ap.add_argument('--cases',type=int,default=48,help='number of cases for EACH dimension')
    ap.add_argument('--only',choices=('both','checked','optimized'),default='both')
    args=ap.parse_args()
    for name,flags in (('checked',['-O0','-g','-fcheck=all','-fbacktrace']),('optimized',['-O3'])):
        if args.only not in ('both',name):continue
        rng=np.random.default_rng(args.seed)
        with tempfile.TemporaryDirectory() as td:
            work=Path(td)
            compiler=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
            for module in SRC:
                run([*compiler,'-I'+str(work),'-J'+str(work),'-c',str(ROOT/'src'/f'{module}.f90'),
                     '-o',str(work/f'{module}.o')],cwd=work)
            exe=work/'auto_nd_driver'
            run([*compiler,'-I'+str(work),str(ROOT/'tools/interpolate_reference/auto_nd_driver.f90'),
                 *(str(work/f'{m}.o') for m in reversed(SRC)),'-o',str(exe)],cwd=work)
            count=0;values=0;worst=0.;fallback=0;ncells=0
            for d in (4,5):
                for i in range(args.cases):
                    n,nv,err,nfallback=check_case(exe,rng,d,i,work)
                    count+=1;values+=nv;worst=max(worst,err);fallback+=nfallback;ncells+=n
            print(f'{name}: {count} native 4D/5D triangulations, {ncells} simplices, '
                  f'{values} real/complex component comparisons, max error {worst:.3g}, '
                  f'{fallback} independent convex-reference fallbacks: PASS',flush=True)

if __name__=='__main__':main()
