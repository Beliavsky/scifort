#!/usr/bin/env python3
"""Measure n-dimensional simplex-BVH query speed, optionally versus a baseline ZIP.

Comparisons use IDENTICAL externally supplied SciPy simplices and queries,
not different automatically constructed triangulations. Reports one-run CPU
measurements, excluding triangulation, construction and program start-up.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import zipfile
import numpy as np
from scipy.spatial import Delaunay

ROOT=Path(__file__).resolve().parents[2]
SRC=('kinds','interpolate_kdtree','interpolate_scattered')


def run(cmd, cwd, stdin=None):
    p=subprocess.run(cmd,cwd=cwd,input=stdin,capture_output=True,text=True)
    if p.returncode:
        raise RuntimeError(f'Failed {cmd}: {p.stderr[-1600:]} {p.stdout[-400:]}')
    return p.stdout


def build(root,dest):
    cc=['gfortran','-std=f2018','-ffree-line-length-none','-O3']
    for module in SRC:
        run([*cc,'-I'+str(dest),'-J'+str(dest),'-c',str(root/'src'/f'{module}.f90'),
             '-o',str(dest/f'{module}.o')],dest)
    exe=dest/'bench'
    run([*cc,'-I'+str(dest),str(ROOT/'benchmarks/bench_interpolate_simplex_bvh.f90'),
         *(str(dest/f'{module}.o') for module in reversed(SRC)),'-o',str(exe)],dest)
    return exe


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--baseline-zip',type=Path)
    ap.add_argument('--queries',type=int,default=2500)
    ap.add_argument('--repetitions',type=int,default=3)
    args=ap.parse_args()
    rng=np.random.default_rng(20261031)
    with tempfile.TemporaryDirectory() as temp:
        td=Path(temp)
        newdir=td/'new';newdir.mkdir()
        cur=build(ROOT,newdir)
        older=None
        if args.baseline_zip:
            with zipfile.ZipFile(args.baseline_zip) as archive:
                archive.extractall(td/'oldroot')
            entries=list((td/'oldroot').glob('*/src/interpolate_scattered.f90'))
            if len(entries)!=1:raise RuntimeError('Could not find baseline source in ZIP')
            olddir=td/'oldbuild';olddir.mkdir()
            older=build(entries[0].parents[1],olddir)
        for dim,n in ((4,65),(5,50)):
            pts=rng.random((n,dim))
            tri=Delaunay(pts)
            vals=rng.normal(size=n)
            ids=rng.integers(len(tri.simplices),size=args.queries)
            weights=rng.dirichlet(np.ones(dim+1),size=args.queries)
            q=np.einsum('ij,ijk->ik',weights,pts[tri.simplices[ids]])
            rows=[f'{n} {dim} {len(tri.simplices)} {len(q)} {args.repetitions}']
            rows.extend(' '.join(f'{v:.17g}' for v in row) for row in pts)
            rows.append(' '.join(f'{v:.17g}' for v in vals))
            rows.extend(' '.join(str(int(v)+1) for v in row) for row in tri.simplices)
            rows.extend(' '.join(f'{v:.17g}' for v in row) for row in q)
            inp='\n'.join(rows)+'\n'
            ti,ch=map(float,run([str(cur)],newdir,inp).split())
            print(f'{dim}D / {n} sites / {len(tri.simplices)} simplices / '
                  f'{len(q)} queries x {args.repetitions}: BVH {ti:.6g}s; checksum {ch:.12g}',end='')
            if older:
                tb,cb=map(float,run([str(older)],older.parent,inp).split())
                if abs(ch-cb)>1e-8*max(1,abs(ch)):
                    raise RuntimeError(f'mismatched checksums: {ch} vs {cb}')
                print(f'; baseline {tb:.6g}s; speedup {tb/ti:.1f}x',flush=True)
            else:
                print(flush=True)

if __name__=='__main__':main()
