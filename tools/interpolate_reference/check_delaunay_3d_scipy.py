#!/usr/bin/env python3
"""Check native 3D tetrahedra against Qhull hull volume and topology.

Randomized general-position comparisons check convex-hull completeness and
number/index validity. Exact simplex sets are also checked on general-position clouds. Qhull inputs
are centered before triangulation to avoid large-offset roundoff artifacts.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.spatial import ConvexHull, Delaunay

ROOT=Path(__file__).resolve().parents[2]


def run(command: list[str], cwd: Path, data: str | None = None) -> str:
    p=subprocess.run(command,input=data,cwd=cwd,capture_output=True,text=True)
    if p.returncode:
        raise RuntimeError('Command failed: '+ ' '.join(command) + '\n' + p.stderr[-2000:])
    return p.stdout


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=83455)
    parser.add_argument('--cases',type=int,default=150)
    parser.add_argument('--only',choices=('checked','optimized','both'),default='both')
    args=parser.parse_args()
    for mode,flags in (('checked',['-O0','-g','-fcheck=all','-fbacktrace']),
                       ('optimized',['-O3'])):
        if args.only!='both' and args.only!=mode:
            continue
        rng=np.random.default_rng(args.seed)
        relative_worst=0.0
        exact_meshes=0
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td)
            common=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
            for name in ('kinds','interpolate_kdtree','interpolate_scattered'):
                run([*common,'-I'+str(tmp),'-J'+str(tmp),'-c',
                     str(ROOT/'src'/f'{name}.f90'),'-o',str(tmp/f'{name}.o')],tmp)
            exe=tmp/'delaunay_3d_driver'
            run([*common,'-I'+str(tmp),str(ROOT/'tools/interpolate_reference/delaunay_3d_driver.f90'),
                 *[str(tmp/f'{name}.o') for name in ('interpolate_scattered','interpolate_kdtree','kinds')],
                 '-o',str(exe)],tmp)
            for k in range(args.cases):
                n=(4,5,7,12,20,35,60,100,150)[k%9]
                points=rng.normal(size=(n,3))
                if k%8==0:
                    points*=np.array([10.0,1.0,0.1])
                if k%11==0:
                    points+=1.0e6
                if k%17==0:
                    points[:,2]*=0.003
                data='\n'.join([str(n), *[' '.join(f'{v:.17g}' for v in p) for p in points]])+'\n'
                result=run([str(exe)],tmp,data).splitlines()
                if int(result[0])!=0:
                    raise AssertionError(f'{mode} case {k}: construction status {result[0]}')
                nt=int(result[1])
                simplices=np.array([[int(v)-1 for v in line.split()] for line in result[2:]],int)
                if simplices.shape!=(nt,4) or np.any(simplices<0) or np.any(simplices>=n):
                    raise AssertionError(f'{mode} case {k}: invalid simplex shape or vertex index')
                verts=points[simplices]
                edges=np.stack((verts[:,1]-verts[:,0], verts[:,2]-verts[:,0],
                                verts[:,3]-verts[:,0]),axis=1)
                vol=float(np.sum(np.abs(np.linalg.det(edges))/6.0))
                reference=ConvexHull(points).volume
                relerr=abs(vol-reference)/reference
                relative_worst=max(relative_worst,relerr)
                # Large coordinate offsets intentionally cause some roundoff
                # in *independent* double-precision hull-volume estimates.
                if relerr>1e-8:
                    raise AssertionError(f'{mode} case {k}: convex hull volume error {relerr:.3g}')
                expected={tuple(sorted(s)) for s in Delaunay(points-points.mean(axis=0)).simplices}
                actual={tuple(sorted(s)) for s in simplices}
                if actual!=expected:
                    raise AssertionError(f'{mode} case {k}: different Delaunay mesh from centered Qhull')
                exact_meshes+=1
        print(f'{mode}: {args.cases} complete 3D hulls, {exact_meshes} exact Qhull meshes; '
              f'max relative volume error {relative_worst:.3g} PASS',flush=True)

if __name__=='__main__':
    main()
