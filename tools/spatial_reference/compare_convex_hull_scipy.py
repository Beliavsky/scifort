#!/usr/bin/env python3
"""Independent SciPy geometry reference for native 2D/3D ConvexHull.

Example:
    python tools/spatial_reference/compare_convex_hull_scipy.py --mode both --count 200
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial import ConvexHull, QhullError

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("kinds", "interpolate_kdtree", "spatial_distance", "spatial_hausdorff",
           "spatial_slerp", "spatial_procrustes", "spatial_voronoi", "spatial_convex_hull", "spatial_halfspace", "spatial")


def compile_driver(folder: Path, flags: list[str]) -> Path:
    fc = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    objects = []
    for source in SOURCES:
        output = folder / f"{source}.o"
        subprocess.run([*fc, "-I", str(folder), "-J", str(folder), "-c",
                        str(ROOT / "src" / f"{source}.f90"), "-o", str(output)],
                       cwd=folder, capture_output=True, check=True)
        objects.append(str(output))
    executable = folder / "convex_hull_driver"
    subprocess.run([*fc, "-I", str(folder),
                    str(ROOT / "tools/spatial_reference/convex_hull_driver.f90"),
                    *objects, "-o", str(executable)], cwd=folder, capture_output=True, check=True)
    return executable


def make_cases(seed: int, count: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    cases: list[np.ndarray] = []
    for k in range(count):
        d = 2 if k % 2 == 0 else 3
        n = int(rng.integers(d + 2, 110 if d == 3 else 160))
        mode = (k // 2) % 12
        if mode == 0:
            points = rng.uniform(-1, 1, (n, d))
        elif mode == 1:
            points = rng.normal(size=(n, d))
        elif mode == 2:
            points = rng.normal(size=(n, d))
            points /= np.linalg.norm(points, axis=1, keepdims=True)
        elif mode == 3:
            points = rng.uniform(-1, 1, (n, d))
            points *= rng.choice([1.0, .25, 2.0], size=d)
        elif mode == 4:
            points = rng.uniform(0, 1, (n, d))
            points = points * rng.uniform(.1, 5, size=d) + rng.normal(size=d)
        elif mode == 5:
            # Many interior points and a simplex containing them.
            points = rng.uniform(-.15, .15, (n, d))
            corners = np.eye(d + 1, d)
            corners[0] = -1
            points[:d+1] = corners * 2
        elif mode == 6:
            points = rng.normal(size=(n, d))
            if n > 12:
                points[:5] = points[0]
        elif mode == 7:
            # Clustered 3D clouds and arc-like 2D clouds (nondegenerate).
            points = rng.normal(size=(n,d)) * (np.array([1., .2, .04])[:d])
        elif mode == 8:
            points = rng.uniform(-1,1,(n,d))
            points *= rng.uniform(.2,2,size=d)
            points += 10000 * rng.normal(size=d)
        elif mode == 9:
            points = rng.normal(size=(n,d))
            points *= 1e-7
        elif mode == 10:
            points = rng.normal(size=(n,d))
            points *= 1e5
        else:
            points = rng.normal(size=(n,d))
            if n > 15:
                points[:n//3] = points[0]
        # Qhull can fail to select an initial full-rank simplex when the
        # input contains many exact duplicates, even for a valid hull.
        # Use a fresh nondegenerate sample for oracle comparisons instead.
        try:
            ConvexHull(points)
        except QhullError:
            points = rng.normal(size=(max(12,n), d))
        cases.append(points)
    cases.append(np.array([[x,y] for x in (-1.,1.) for y in (-1.,1.)]))
    cases.append(np.array([[x,y,z] for x in (-1.,1.) for y in (-1.,1.)
                           for z in (-1.,1.)]))
    cases.append(np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.],[.5,0.],[.5,.5]]))
    # Regression: near-coplanar 3D dual hull previously rejected due to
    # unmatched facet edges in the incremental construction.
    dual = np.loadtxt(Path(__file__).with_name("convex_hull_nearcoplanar_26.txt"))
    cases.append(dual)
    for k in range(8):
        cases.append(dual[rng.permutation(len(dual))])
    return cases


def run_cases(executable: Path, cases: list[np.ndarray]) -> tuple[int,float,float]:
    data = [str(len(cases))]
    for pts in cases:
        n,d = pts.shape
        data.append(f"{n} {d}")
        data.extend(" ".join(f"{x:.17e}" for x in row) for row in pts)
    process = subprocess.run([str(executable)], input="\n".join(data)+"\n", text=True,
                             capture_output=True, check=True)
    stdout = iter(process.stdout.splitlines())
    max_area = max_volume = 0.0
    faces_checked = 0
    for case_id, pts in enumerate(cases):
        d = pts.shape[1]
        ref = ConvexHull(pts)
        status,nvertices,nfaces = map(int,next(stdout).split())
        if status != 0:
            raise AssertionError((case_id, 'Fortran rejected SciPy hull',status,pts.shape))
        area,volume = map(float,next(stdout).split())
        vertices = np.array([int(i)-1 for i in next(stdout).split()],dtype=int)
        if nvertices != len(ref.vertices):
            raise AssertionError((case_id,'vertex count',nvertices,len(ref.vertices),pts.shape))
        if sorted(vertices) != sorted(ref.vertices):
            # Qhull and SciFort may legitimately choose different input
            # indices for *identical* duplicate points.
            ours = pts[vertices]
            theirs = pts[ref.vertices]
            ds = np.linalg.norm(ours[:,None,:]-theirs[None,:,:],axis=2)
            tolerance = 1e-12*max(1.,np.ptp(pts,axis=0).max())
            if max(ds.min(axis=0).max(),ds.min(axis=1).max()) > tolerance:
                raise AssertionError((case_id,'vertex coordinates',vertices,ref.vertices,pts.shape))
        # Area means perimeter in 2D, and surface area in 3D.
        area_error = abs(area-ref.area)/max(abs(ref.area), 1e-300)
        volume_error = abs(volume-ref.volume)/max(abs(ref.volume), 1e-300)
        max_area = max(max_area,area_error)
        max_volume = max(max_volume,volume_error)
        if area_error > 2e-9 or volume_error > 2e-9:
            raise AssertionError((case_id,'area volume', area,volume,
                                  ref.area,ref.volume,area_error,volume_error))
        eq = []
        faces = []
        adj = []
        for _ in range(nfaces):
            faces.append(np.array([int(v)-1 for v in next(stdout).split()]))
            eq.append(np.array([float(v) for v in next(stdout).split()]))
            adj.append(np.array([int(v)-1 for v in next(stdout).split()]))
        faces = np.stack(faces)
        eq = np.stack(eq)
        adj = np.stack(adj)
        faces_checked += len(faces)
        if faces.shape != (nfaces,d) or eq.shape != (nfaces,d+1):
            raise AssertionError((case_id,'shape'))
        norm = np.linalg.norm(eq[:,:d],axis=1)
        if np.max(abs(norm-1)) > 2e-12:
            raise AssertionError((case_id,'nonunit outward normals'))
        scale = max(1.,np.ptp(pts,axis=0).max())
        if np.max(eq[:,:d]@pts.T+eq[:,d,None]) > 1e-10*scale:
            raise AssertionError((case_id,'uncontained site'))
        if np.any(adj<0) or np.any(adj>=nfaces):
            raise AssertionError((case_id,'missing adjacency'))
        for facet in range(nfaces):
            for vertex in range(d):
                other = adj[facet,vertex]
                # Each edge in 3D (point in 2D) shares d-1 vertices.
                common = set(faces[facet]) & set(faces[other])
                if common != (set(faces[facet])-{faces[facet,vertex]}):
                    raise AssertionError((case_id,'neighbor geometry',facet,vertex))
        if d == 2:
            p = pts[vertices]
            signed_area = np.sum(p[:,0]*np.roll(p[:,1],-1) -
                                 p[:,1]*np.roll(p[:,0],-1))
            if signed_area <= 0:
                raise AssertionError((case_id,'non-CCW hull'))
        else:
            if nfaces != 2*nvertices-4:
                raise AssertionError((case_id,'Euler triangulation',nfaces,nvertices))
    return faces_checked,max_area,max_volume


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--seed',type=int,default=6919)
    parser.add_argument('--count',type=int,default=200)
    args = parser.parse_args()
    cases = make_cases(args.seed,args.count)
    modes = []
    if args.mode != 'optimized':
        modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.mode != 'checked':
        modes.append(('optimized',['-O3']))
    for label,flags in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-hull-') as directory:
            executable = compile_driver(Path(directory),flags)
            total,da,dv = run_cases(executable,cases)
            print(f'{label}: {len(cases)} SciPy ConvexHull cases, {total} facet records '
                  f'PASS; max relative area {da:.3g}, volume {dv:.3g}',flush=True)


if __name__ == '__main__':
    main()
