#!/usr/bin/env python3
"""Compare SciFort's balanced KD-tree against scipy.spatial.cKDTree.

This checks exact k-nearest indices and distances and closed-radius queries
in 1..7 dimensions, with 1-based Fortran indices and squared distances.
Tie-breaking on duplicate sites is checked by the native Fortran test instead;
SciPy may return a different equally near site for exact ties.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[2]


def run(cmd,cwd,input=None):
    p=subprocess.run(cmd,cwd=cwd,input=input,text=True,capture_output=True)
    if p.returncode:
        raise RuntimeError(f"Command failed: {' '.join(map(str,cmd))}\n{p.stderr[-3000:]}\n{p.stdout[-1000:]}")
    return p.stdout


def build(tmp,flags):
    fc=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    for name in ('kinds','interpolate_kdtree'):
        run([*fc,'-J'+str(tmp),'-I'+str(tmp),'-c',str(ROOT/'src'/f'{name}.f90'),
             '-o',str(tmp/f'{name}.o')],tmp)
    exe=tmp/'kdtree_driver'
    run([*fc,'-I'+str(tmp),str(ROOT/'tools/interpolate_reference/kdtree_driver.f90'),
         str(tmp/'interpolate_kdtree.o'),str(tmp/'kinds.o'),'-o',str(exe)],tmp)
    return exe


def one_case(exe,rng,k):
    dim=1+k%7
    n=20+(k*23)%230
    neighbors=min(n,1+(k*7)%40)
    count=3+(k%7)
    points=rng.normal(size=(n,dim))*rng.uniform(0.1,3.0,size=(1,dim))
    query=rng.normal(size=(count,dim))
    query[0]=points[k%n]  # zero-distance exact match
    radius=float(rng.uniform(0.15,1.5))
    lines=[f'{n} {dim} {count} {neighbors} {radius:.17g}']
    lines.extend(' '.join(f'{v:.17g}' for v in row) for row in points)
    lines.extend(' '.join(f'{v:.17g}' for v in row) for row in query)
    output=iter(run([str(exe)],exe.parent,'\n'.join(lines)+'\n').splitlines())
    ref=cKDTree(points)
    distances,indices=ref.query(query,k=neighbors)
    distances=np.asarray(distances).reshape(count,neighbors)
    indices=np.asarray(indices).reshape(count,neighbors)
    comparisons=0
    peak=0.0
    for j in range(count):
        found=[]
        squared=[]
        for _ in range(neighbors):
            si,di=next(output).split()
            found.append(int(si)-1)
            squared.append(float(di))
        found=np.asarray(found)
        squared=np.asarray(squared)
        if not np.array_equal(indices[j],found):
            raise AssertionError(f'case {k} query {j}: different nearest ids: {found} vs {indices[j]}')
        err=float(np.max(np.abs(squared-distances[j]**2)))
        peak=max(err,peak)
        if err>2e-12*(1+max(squared)):
            raise AssertionError(f'case {k} query {j}: squared distances differ: {err}')
        comparisons+=neighbors
        ball=np.array([int(next(output)) for _ in range(int(next(output)))],dtype=int)-1
        reference_ball=np.asarray(sorted(ref.query_ball_point(query[j],radius)),dtype=int)
        if not np.array_equal(ball,reference_ball):
            raise AssertionError(f'case {k} query {j}: radius indices mismatch')
        comparisons+=len(ball)
    return comparisons,peak


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--seed',type=int,default=20261031)
    ap.add_argument('--cases',type=int,default=180)
    ap.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    args=ap.parse_args()
    modes=[]
    if args.mode!='optimized': modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.mode!='checked': modes.append(('optimized',['-O3']))
    for label,flags in modes:
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td)
            exe=build(tmp,flags)
            rng=np.random.default_rng(args.seed)
            comparisons=0
            maxerr=0
            for k in range(args.cases):
                n,err=one_case(exe,rng,k)
                comparisons+=n
                maxerr=max(maxerr,err)
            print(f'{label}: PASS {args.cases} KD-tree scenarios; '
                  f'{comparisons} individual neighbor/radius indices; max squared-distance error {maxerr:.3g}',flush=True)

if __name__=='__main__':
    main()
