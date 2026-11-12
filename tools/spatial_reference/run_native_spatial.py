#!/usr/bin/env python3
"""Compile and run all self-contained native scipy.spatial regressions.

Runs with bounds checks and optimizations; no dependency on FPM/BLAS/LAPACK.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("kinds", "interpolate_kdtree", "spatial_distance", "spatial_hausdorff", "spatial_slerp", "spatial_procrustes", "spatial_voronoi", "spatial_convex_hull", "spatial_halfspace", "spatial_transform", "spatial_rotation_spline", "spatial_rotation_stack", "spatial_rotation_grid", "spatial_rotation_group", "spatial_rigid_transform", "spatial_rigid_transform_grid", "spatial_transform_nd", "spatial")
TESTS = ("test_spatial_distance", "test_spatial_kdtree_pairs", "test_spatial_kdtree_periodic", "test_spatial_kdtree_weighted_count", "test_spatial_distance_covariance", "test_spatial_distance_boolean", "test_spatial_hausdorff", "test_spatial_slerp", "test_spatial_procrustes", "test_spatial_voronoi", "test_spatial_convex_hull", "test_spatial_convex_hull_nearcoplanar", "test_spatial_halfspace", "test_spatial_rotation", "test_spatial_rotation_slerp", "test_spatial_rotation_euler_mrp", "test_spatial_rotation_davenport", "test_spatial_rotation_estimation", "test_spatial_rotation_spline", "test_spatial_rotation_stack", "test_spatial_rotation_random_grid", "test_spatial_rotation_group", "test_spatial_rigid_transform", "test_spatial_rigid_transform_grid", "test_spatial_transform_nd")


def run(command: list[str], cwd: Path) -> str:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError("Failed: " + " ".join(command) +
                           "\nstdout:\n" + result.stdout[-2500:] +
                           "\nstderr:\n" + result.stderr[-3000:])
    return result.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    args = parser.parse_args()
    modes = []
    if args.mode != "optimized":
        modes.append(("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]))
    if args.mode != "checked":
        modes.append(("optimized", ["-O3"]))
    for label, flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-spatial-native-") as td:
            tmp = Path(td)
            fc = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
            objects = []
            for name in SOURCES:
                obj = tmp / f"{name}.o"
                run([*fc, "-I", str(tmp), "-J", str(tmp), "-c",
                     str(ROOT / "src" / f"{name}.f90"), "-o", str(obj)], tmp)
                objects.append(str(obj))
            for name in TESTS:
                exe = tmp / name
                run([*fc, "-I", str(tmp), str(ROOT / "test" / f"{name}.f90"),
                     *objects, "-o", str(exe)], tmp)
                print(f"{label}: {name}: {run([str(exe)], tmp)}", flush=True)
        print(f"{label}: {len(TESTS)}/{len(TESTS)} native spatial programs PASS", flush=True)


if __name__ == "__main__":
    main()
