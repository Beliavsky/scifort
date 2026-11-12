#!/usr/bin/env python3
"""Validate SciFort's new KD-tree pair APIs independently against SciPy.

Compile checked and optimized Fortran with no external libraries; check
query_pairs, query_ball_tree (CSR), count_neighbors, sparse_distance_matrix.

    python tools/spatial_reference/compare_kdtree_pairs_scipy.py
    python tools/spatial_reference/compare_kdtree_pairs_scipy.py --seed 8675309
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]


def build(folder: Path, flags: list[str]) -> Path:
    fc = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    sources = ("kinds.f90", "interpolate_kdtree.f90", "spatial_distance.f90", "spatial_hausdorff.f90", "spatial_slerp.f90", "spatial_procrustes.f90", "spatial_voronoi.f90", "spatial_convex_hull.f90", "spatial_halfspace.f90", "spatial.f90")
    objects = []
    for source in sources:
        obj = folder / source.replace(".f90", ".o")
        subprocess.run([*fc, "-J", str(folder), "-I", str(folder), "-c",
                        str(ROOT / "src" / source), "-o", str(obj)], cwd=folder, check=True)
        objects.append(str(obj))
    exe = folder / "kdtree_pairs_driver"
    subprocess.run([*fc, "-I", str(folder),
                    str(ROOT / "tools/spatial_reference/kdtree_pairs_driver.f90"),
                    *objects, "-o", str(exe)], cwd=folder, check=True)
    return exe


def run_case(exe: Path, points: np.ndarray, other: np.ndarray,
             radius: float, case_num: int) -> tuple[int, float]:
    n, d = points.shape
    m = len(other)
    lines = [f"{n} {m} {d} {radius:.17g}"]
    lines += [" ".join(f"{v:.17g}" for v in row) for row in points]
    lines += [" ".join(f"{v:.17g}" for v in row) for row in other]
    proc = subprocess.run([str(exe)], input="\n".join(lines) + "\n", text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    output = iter(proc.stdout.split())

    def header(key: str) -> int:
        got = next(output)
        if got != key:
            raise AssertionError(f"case {case_num}: expected {key}, got {got}")
        return int(next(output))

    first = cKDTree(points)
    second = cKDTree(other)
    npairs = header("P")
    pairs = {(int(next(output))-1, int(next(output))-1) for _ in range(npairs)}
    expected_pairs = first.query_pairs(radius)
    if pairs != expected_pairs or len(pairs) != npairs:
        raise AssertionError(f"case {case_num}: query_pairs mismatch: {pairs ^ expected_pairs}")

    nrows = header("R")
    row_ptr = [int(next(output)) for _ in range(nrows)]
    ncols = header("J")
    indices = [int(next(output))-1 for _ in range(ncols)]
    if nrows != n + 1 or row_ptr[0] != 1 or row_ptr[-1] != ncols+1:
        raise AssertionError(f"case {case_num}: invalid CSR row pointers")
    cross_expected = first.query_ball_tree(second, radius)
    for i, expected in enumerate(cross_expected):
        got = indices[row_ptr[i]-1:row_ptr[i+1]-1]
        if got != sorted(expected):
            raise AssertionError(f"case {case_num}, row {i}: query_ball_tree mismatch")

    if header("C") != int(first.count_neighbors(second, radius)):
        raise AssertionError(f"case {case_num}: count_neighbors mismatch")
    nentries = header("S")
    matrix = {(int(next(output))-1, int(next(output))-1): float(next(output))
              for _ in range(nentries)}
    reference = first.sparse_distance_matrix(second, radius, output_type="dict")
    if set(matrix) != set(reference) or nentries != len(reference):
        raise AssertionError(f"case {case_num}: sparse matrix index mismatch")
    if list(matrix) != sorted(matrix):
        raise AssertionError(f"case {case_num}: sparse matrix not in row-major order")
    worst = 0.0
    for key, value in matrix.items():
        ref = reference[key]
        rel = abs(value-ref)/max(1.0, abs(ref))
        if rel > 4e-13:
            raise AssertionError(f"case {case_num}: sparse matrix {key} {value} vs {ref}")
        worst = max(worst,rel)
    if next(output, None) is not None:
        raise AssertionError(f"case {case_num}: trailing reference driver output")
    return npairs + ncols + nentries, worst


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--cases", type=int, default=200)
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    args = parser.parse_args()
    modes = []
    if args.mode != "optimized":
        modes.append(("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]))
    if args.mode != "checked":
        modes.append(("optimized", ["-O3"]))
    for label, flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-kdtree-pairs-") as tmp:
            exe = build(Path(tmp), flags)
            rng = np.random.default_rng(args.seed)
            total = 0
            worst = 0.0
            for scenario in range(args.cases):
                ndim = 1 + scenario % 7
                n = 3 + (scenario*17) % 54
                m = 2 + (scenario*11) % 46
                scale = rng.uniform(0.1, 3.0, size=ndim)
                xa = rng.normal(size=(n,ndim))*scale
                xb = rng.normal(size=(m,ndim))*scale
                if scenario % 4 == 0:  # duplicates within and across trees
                    xa[1] = xa[0]
                    xb[0] = xa[0]
                elif scenario % 4 == 1:  # exact lattice and ties
                    xa = rng.integers(-2, 3, size=(n,ndim)).astype(float)
                    xb = rng.integers(-2, 3, size=(m,ndim)).astype(float)
                radius = (0.0, 0.25, 1.0, 3.0,
                          float(np.sqrt(np.sum((np.ptp(xa,axis=0)+np.ptp(xb,axis=0))**2))+1.0)
                         )[scenario % 5]
                comparisons, err = run_case(exe, xa, xb, radius, scenario)
                total += comparisons
                worst = max(worst,err)
            print(f"{label}: PASS {args.cases} scenarios, {total} matched pair/matrix records, "
                  f"max scale-relative distance error {worst:.3g}", flush=True)


if __name__ == "__main__":
    main()
