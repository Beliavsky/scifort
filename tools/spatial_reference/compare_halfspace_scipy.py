#!/usr/bin/env python3
"""Compare SciFort's bounded 2D/3D halfspace intersections with SciPy/Qhull.

Example:
    python tools/spatial_reference/compare_halfspace_scipy.py --mode both --count 150
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial import HalfspaceIntersection, ConvexHull, QhullError
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("kinds", "interpolate_kdtree", "spatial_distance", "spatial_hausdorff",
           "spatial_slerp", "spatial_procrustes", "spatial_voronoi", "spatial_convex_hull",
           "spatial_halfspace", "spatial")


def compile_driver(folder: Path, flags: list[str]) -> Path:
    fc = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    objects = []
    for source in SOURCES:
        output = folder / f"{source}.o"
        run = subprocess.run([*fc, "-I", str(folder), "-J", str(folder), "-c",
                              str(ROOT / "src" / f"{source}.f90"), "-o", str(output)],
                             cwd=folder, capture_output=True, text=True)
        if run.returncode:
            raise RuntimeError(run.stderr)
        objects.append(str(output))
    executable = folder / "halfspace_driver"
    run = subprocess.run([*fc, "-I", str(folder),
                          str(ROOT / "tools/spatial_reference/halfspace_driver.f90"),
                          *objects, "-o", str(executable)],
                         cwd=folder, capture_output=True, text=True)
    if run.returncode:
        raise RuntimeError(run.stderr)
    return executable


def make_cases(seed: int, count: int) -> list[tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(seed)
    cases = []
    for k in range(count):
        d = 2 + k % 2
        scale = 10.0 ** rng.uniform(-2.0, 2.0)
        center = rng.normal(size=d) * rng.choice([0.0, 1.0, 100.0])
        if (k // 2) % 3 == 0:
            # Random supporting planes, with explicit bounding-box constraints.
            n = int(rng.integers(5, 34))
            a = rng.normal(size=(n,d))
            length = np.linalg.norm(a,axis=1)
            a /= length[:,None]
            a = np.vstack((a, np.eye(d), -np.eye(d)))
            radius = rng.uniform(.6, 2.5, len(a)) * scale
            a *= rng.uniform(.25, 4.0, size=(len(a),1))
            b = -(a @ center) - radius*np.linalg.norm(a,axis=1)
            hs = np.column_stack([a, b])
            if (k // 2) % 6 == 0:
                hs = np.vstack([hs, hs[0], np.r_[hs[1,:d],hs[1,d]-scale*3]])
        elif (k // 2) % 3 == 1:
            # Supporting halfspaces from a random convex polytope.
            points = center + scale*rng.normal(size=(int(rng.integers(6,40)),d))
            ch = ConvexHull(points)
            hs = ch.equations.copy()
            center = points.mean(axis=0)
        else:
            # Rotated and scaled box, including non-simplicial duals.
            q,_ = np.linalg.qr(rng.normal(size=(d,d)))
            a = np.vstack((q, -q))
            radius = rng.uniform(.4, 2.0, size=d)*scale
            b = -(a @ center) - np.r_[radius, radius]
            hs = np.column_stack([a,b])
        cases.append((hs,center))
    # Exact configurations with degeneracy, redundant faces and different scales.
    for d in (2,3):
        a = np.vstack((np.eye(d),-np.eye(d)))
        center = np.arange(1,d+1)*31.0
        hs = np.column_stack((a, -(a @ center)-np.ones(2*d)))
        cases.append((hs,center))
    return cases


def run_cases(executable: Path, cases: list[tuple[np.ndarray,np.ndarray]]) -> tuple[float,float,int,int]:
    text_input = [str(len(cases))]
    for h,p in cases:
        n,d = h.shape[0],h.shape[1]-1
        text_input.append(f"{n} {d}")
        text_input.append(" ".join(f"{v:.17e}" for v in p))
        text_input += [" ".join(f"{v:.17e}" for v in row) for row in h]
    run = subprocess.run([str(executable)], input="\n".join(text_input)+"\n", text=True,
                         capture_output=True, timeout=160)
    if run.returncode:
        raise RuntimeError(f"Fortran driver failed: {run.stderr[-2500:]}")
    lines = iter(run.stdout.splitlines())
    max_abs_error, max_scaled_error = 0.0, 0.0
    checked = facet_count = 0
    for k,(h,center) in enumerate(cases):
        try:
            reference = HalfspaceIntersection(h,center)
        except QhullError as exc:
            raise AssertionError((k,"unexpected SciPy failure",str(exc))) from exc
        status,nv,nfacets = map(int,next(lines).split())
        if status != 0:
            raise AssertionError((k,'Fortran rejected bounded polytope',h.shape))
        dual_area,dual_volume = map(float,next(lines).split())
        coords = np.array([[float(v) for v in next(lines).split()] for _ in range(nv)])
        duals = np.array([[float(v) for v in next(lines).split()] for _ in range(len(h))])
        expected = np.asarray(reference.intersections)
        # SciPy/Qhull can return repeated vertices when several simplicial
        # polar facets lie on the same polygonal face. Compare unique sites,
        # not Qhull's triangulation-dependent multiplicity.
        ref_span = max(1e-300,np.ptp(expected,axis=0).max())
        unique = []
        for p in expected:
            if not any(np.max(np.abs(p-v)) < 2e-8*ref_span for v in unique):
                unique.append(p)
        if len(coords) != len(unique):
            raise AssertionError((k,'unique vertex count',len(coords),len(unique),
                                  len(expected),h.shape))
        ds = cdist(coords,expected)
        error = max(ds.min(axis=0).max(),ds.min(axis=1).max())
        ref_span = max(1e-300,np.ptp(expected,axis=0).max())
        scaled_error = error/ref_span
        max_abs_error = max(max_abs_error,error)
        max_scaled_error = max(max_scaled_error,scaled_error)
        # 3D hull triangulation need not agree in facet indexing or order.
        if scaled_error > 2e-8:
            raise AssertionError((k,'intersection coordinates',error,scaled_error,
                                  h.shape,coords,expected))
        qerr = np.max(np.abs(duals-reference.dual_points))
        if qerr > 2e-8*max(1.0,np.max(np.abs(reference.dual_points))):
            raise AssertionError((k,'dual coordinates',qerr))
        for key,ours,theirs in (('dual area',dual_area,reference.dual_area),
                                ('dual volume',dual_volume,reference.dual_volume)):
            if abs(ours-theirs)> 4e-8*max(abs(theirs),1e-300):
                raise AssertionError((k,key,ours,theirs))
        if len(coords) and np.max((h[:,:-1] @ coords.T + h[:,-1,None])) > \
                5e-8*max(1.0,np.max(abs(h[:,-1]))):
            raise AssertionError((k,'violated inequality'))
        checked += 1
        facet_count += nfacets
    return max_abs_error,max_scaled_error,checked,facet_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'),default='both')
    parser.add_argument('--count', type=int, default=200)
    parser.add_argument('--seed',type=int,default=20261102)
    options = parser.parse_args()
    cases = make_cases(options.seed,options.count)
    mode_flags = {
        'checked': ['-O0','-g','-fcheck=bounds','-fbacktrace'],
        'optimized': ['-O3','-march=native']
    }
    modes = list(mode_flags) if options.mode=='both' else [options.mode]
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix=f'scifort_halfspace_{mode}_') as temp:
            executable = compile_driver(Path(temp),mode_flags[mode])
            ae,se,checked,facets = run_cases(executable,cases)
            print(f'{mode}: PASS cases={checked} facets={facets} '
                  f'max_abs_vertex_error={ae:.4g} max_relative_vertex_error={se:.4g}')


if __name__=='__main__':
    main()
