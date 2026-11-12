#!/usr/bin/env python3
"""Independent SciPy oracle for SciFort spherical Voronoi construction.

Compares the geometric vertex sets and the area of each generator's region.
Facet order and the cyclic starting position may legitimately differ.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial import SphericalVoronoi, cKDTree

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("kinds", "interpolate_kdtree", "spatial_distance", "spatial_hausdorff",
           "spatial_slerp", "spatial_procrustes", "spatial_voronoi", "spatial_convex_hull", "spatial_halfspace", "spatial")


def compile_driver(folder: Path, flags: list[str]) -> Path:
    compiler = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    objects = []
    for source in SOURCES:
        output = folder / f"{source}.o"
        subprocess.run([*compiler, "-I", str(folder), "-J", str(folder), "-c",
                        str(ROOT / "src" / f"{source}.f90"), "-o", str(output)],
                       cwd=folder, capture_output=True, check=True)
        objects.append(str(output))
    exe = folder / "voronoi_driver"
    subprocess.run([*compiler, "-I", str(folder),
                    str(ROOT / "tools/spatial_reference/voronoi_driver.f90"),
                    *objects, "-o", str(exe)], cwd=folder, capture_output=True, check=True)
    return exe


def cases(seed: int, count: int) -> list[tuple[np.ndarray, float, np.ndarray]]:
    rng = np.random.default_rng(seed)
    result = []
    for k in range(count):
        d = 2 if k % 5 == 0 else 3
        n = int(rng.integers(d+1, 81 if d == 3 else 130))
        mode = (k // 5) % 8
        if mode == 0 and d == 3:
            # Clustered hemisphere input: center of the sphere need not be in hull.
            x = rng.normal(size=(n,3)) * .35 + np.array([3.0,0.0,0.0])
        elif mode == 1 and d == 3:
            # Slightly clustered points with small spherical cells.
            x = rng.normal(size=(n,3)) * .01 + np.array([1.0,0.0,0.0])
        elif mode == 2 and d == 3:
            # Highly anisotropic sample distribution.
            x = rng.normal(size=(n,3)) * np.array([1.0, .2, .05])
        elif mode == 3 and d == 2:
            # 2-D points concentrated on a short arc of the circle.
            while True:
                t = rng.uniform(0, .35, size=n)
                if np.min(np.diff(np.sort(t))) > 2e-5:
                    break
            x = np.column_stack([np.cos(t),np.sin(t)])
        else:
            x = rng.normal(size=(n,d))
        x /= np.linalg.norm(x,axis=1,keepdims=True)
        radius = float((.1,1.,3.,12.,100.)[(k//3)%5])
        center = rng.normal(size=d) * radius * (0 if k % 3 else 2.0)
        result.append((center + radius*x,radius,center))
    return result


def run_cases(exe: Path, data: list[tuple[np.ndarray,float,np.ndarray]]) -> tuple[int,float,float]:
    lines = [str(len(data))]
    for pts,r,c in data:
        n,d = pts.shape
        lines.append(f"{n} {d} {r:.17e} 1e-6")
        lines.append(" ".join(f"{v:.17e}" for v in c))
        lines.extend(" ".join(f"{v:.17e}" for v in row) for row in pts)
    proc = subprocess.run([str(exe)], input="\n".join(lines)+"\n", text=True,
                          capture_output=True, check=True)
    output = iter(proc.stdout.splitlines())
    max_dist = max_area = 0.0
    for k,(pts,r,c) in enumerate(data):
        status,nverts = map(int,next(output).split())
        if status != 0:
            raise AssertionError((k,'invalid Fortran hull',status,pts.shape))
        nv = nverts
        fverts = np.array([[float(z) for z in next(output).split()]
                          for _ in range(nv)]).reshape(nv,pts.shape[1])
        fareas = np.array([float(next(output)) for _ in range(len(pts))])
        fregions = []
        for _ in range(len(pts)):
            row = [int(v) for v in next(output).split()]
            assert len(row)-1 == row[0], (k,'CSR length')
            fregions.append(np.array(row[1:])-1)
        ref = SphericalVoronoi(pts,radius=r,center=c)
        ref.sort_vertices_of_regions() if pts.shape[1] == 3 else None
        rareas = ref.calculate_areas()
        assert nv == len(ref.vertices), (k,'vertex count',nv,len(ref.vertices))
        tree_f = cKDTree(fverts)
        tree_ref = cKDTree(ref.vertices)
        dist = max(float(max(tree_ref.query(fverts)[0])),
                   float(max(tree_f.query(ref.vertices)[0]))) / r
        max_dist = max(max_dist,dist)
        assert dist < 2e-8, (k,'vertices',dist)
        if pts.shape[1] == 2:
            # SciPy 2-D calculate_areas reports hull-edge arc lengths, not
            # per-generator Voronoi region lengths; hull edge order differs.
            rel_area = float(np.max(abs(np.sort(fareas)-np.sort(rareas))))/r
        else:
            rel_area = float(np.max(abs(fareas-rareas)))/(r*r)
        max_area = max(max_area,rel_area)
        assert rel_area < 5e-8, (k,'areas',rel_area, fareas, rareas)
        measure = 2*np.pi*r if pts.shape[1]==2 else 4*np.pi*r*r
        assert abs(sum(fareas)/measure - 1) < 2e-9, (k,'total area',sum(fareas),measure)
        for i,reg in enumerate(fregions):
            assert len(reg) == len(ref.regions[i]), (k,'degree',i,len(reg),len(ref.regions[i]))
            assert np.all((reg>=0)&(reg<nv)), (k,'range')
            coords = fverts[reg]
            reference_coords = ref.vertices[ref.regions[i]]
            a = float(max(cKDTree(reference_coords).query(coords)[0]))/r
            b = float(max(cKDTree(coords).query(reference_coords)[0]))/r
            assert max(a,b) < 2e-8, (k,'region',i,max(a,b))
    assert not list(output), 'Unexpected extra Fortran output'
    return len(data),max_dist,max_area


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=19237)
    parser.add_argument('--count',type=int,default=180)
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    options = parser.parse_args()
    data = cases(options.seed,options.count)
    modes = []
    if options.mode != 'optimized':
        modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if options.mode != 'checked':
        modes.append(('optimized',['-O3']))
    for mode,flags in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-voronoi-') as folder:
            exe = compile_driver(Path(folder),flags)
            count,vertex_error,area_error=run_cases(exe,data)
            print(f'{mode}: {count} SciPy Voronoi cases PASS; '
                  f'max vertex/radius error={vertex_error:.3g}, '
                  f'max area/radius^(d-1) error={area_error:.3g}',flush=True)


if __name__ == '__main__':
    main()
