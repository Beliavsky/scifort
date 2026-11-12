#!/usr/bin/env python3
"""Independent SciPy 1.17 oracle for real-valued Procrustes fits.

Checks the two point-cloud matrices, disparity, optimal orthogonal residual,
unitarity, and nuclear norm. Rank-deficient rotations are compared by
objective values rather than entries because the optimizer is not unique.

Run: python tools/spatial_reference/compare_procrustes_scipy.py --mode both
     python tools/spatial_reference/compare_procrustes_scipy.py --seed 23791
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.linalg import orthogonal_procrustes
from scipy.spatial import procrustes

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("kinds", "interpolate_kdtree", "spatial_distance", "spatial_hausdorff",
           "spatial_slerp", "spatial_procrustes", "spatial_voronoi", "spatial_convex_hull", "spatial_halfspace", "spatial")


def compile_driver(tmp: Path, flags: list[str]) -> Path:
    command = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    objects = []
    for name in SOURCES:
        output = tmp / f"{name}.o"
        subprocess.run([*command, "-I", str(tmp), "-J", str(tmp), "-c",
                        str(ROOT / "src" / f"{name}.f90"), "-o", str(output)],
                       cwd=tmp, check=True, capture_output=True)
        objects.append(str(output))
    exe = tmp / "procrustes_driver"
    subprocess.run([*command, "-I", str(tmp),
                    str(ROOT / "tools/spatial_reference/procrustes_driver.f90"),
                    *objects, "-o", str(exe)], cwd=tmp, check=True, capture_output=True)
    return exe


def make_cases(seed: int, count: int) -> list[tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(seed)
    result = []
    for i in range(count):
        d = (1, 2, 3, 4, 6, 10, 15)[i % 7]
        n = int(rng.integers(2, 30))
        a = rng.normal(size=(n, d))
        b = rng.normal(size=(n, d))
        mode = (i // 7) % 13
        if mode == 0:
            # Rotated and scaled, with independently translated coordinates.
            q, _ = np.linalg.qr(rng.normal(size=(d, d)))
            b = (a @ q) * float(rng.uniform(.04, 15.0)) + rng.normal(size=(1, d))
        elif mode == 1:
            # Reflection and translation.
            q, _ = np.linalg.qr(rng.normal(size=(d, d)))
            q[:, 0] *= -1
            b = a @ q + rng.normal(size=(1, d)) * 25
        elif mode == 2:
            # Identical standardized point configurations.
            b = a.copy()
        elif mode == 3:
            # Nonconstant but rank-one geometry.
            a = rng.normal(size=(n, 1)) @ rng.normal(size=(1, d))
            b = a + 5.0 * rng.normal(size=(1, d))
        elif mode == 4:
            # Both sets rank-deficient in different subspaces.
            rank = max(1, min(d, n-1) // 2)
            a = rng.normal(size=(n, rank)) @ rng.normal(size=(rank, d))
            b = rng.normal(size=(n, rank)) @ rng.normal(size=(rank, d))
        elif mode == 5:
            # Repeated points, which are permitted by SciPy.
            a[1::3] = a[0]
            b[1::4] = b[0]
            if n == 2:
                a[1,0] += 1.0
                b[1,0] += 1.0
        elif mode == 6:
            # Mildly ill-conditioned coordinate configurations.
            a *= 10.0 ** (-np.arange(d) / 1.5)
            b *= 10.0 ** (-np.arange(d) / 1.5)
        elif mode == 7:
            a *= 1e180
            b *= 1e-180
        elif mode == 8:
            a *= 1e-180
            b *= 1e180
        elif mode == 9:
            a *= 1e-150
            b *= 1e-150
        elif mode == 10:
            a *= 1e150
            b *= 1e150
        elif mode == 11:
            # Tiny positive signal riding a translation.
            a = a * 1e-8 + 2.0
            b = b * 1e-8 - 3.0
        else:
            # Near-rank-one data; nonunique rotations can differ.
            a = a[:, :1] @ rng.normal(size=(1, d)) + a * 1e-9
            b = b[:, :1] @ rng.normal(size=(1, d)) + b * 1e-9
        result.append((a, b))
    return result


def run_cases(exe: Path, cases: list[tuple[np.ndarray, np.ndarray]]) -> tuple[int, float, float]:
    lines = [str(len(cases))]
    for a, b in cases:
        lines.append(f"{a.shape[0]} {a.shape[1]}")
        for row in (*a, *b):
            lines.append(" ".join(f"{v:.17e}" for v in row))
    proc = subprocess.run([str(exe)], input="\n".join(lines) + "\n", text=True,
                          capture_output=True, check=True)
    output = iter(proc.stdout.splitlines())
    max_disparity = max_matrix = 0.0
    for k, (a, b) in enumerate(cases):
        n, d = a.shape
        status, disparity = next(output).split()
        assert int(status) == 0, (k, status)
        x = np.array([[float(z) for z in next(output).split()] for _ in range(n)])
        y = np.array([[float(z) for z in next(output).split()] for _ in range(n)])
        status, scale = next(output).split()
        assert int(status) == 0, (k, status)
        r = np.array([[float(z) for z in next(output).split()] for _ in range(d)])
        disparity = float(disparity)
        scale = float(scale)
        # SciPy's raw norm2 over/underflows beyond 1e154/1e-154;
        # independent positive rescalings leave these fits unchanged.
        sa = np.max(abs(a))
        sb = np.max(abs(b))
        a_unit = a / sa
        b_unit = b / sb
        ref_x, ref_y, ref_disp = procrustes(a_unit-a_unit[:1], b_unit-b_unit[:1])
        norm_a = np.linalg.norm(a_unit)
        norm_b = np.linalg.norm(b_unit)
        ref_r, ref_scale_unit = orthogonal_procrustes(a_unit/norm_a, b_unit/norm_b)
        ref_scale = ((ref_scale_unit*norm_a)*sa)*(norm_b*sb)
        disp_err = abs(disparity - ref_disp)
        max_disparity = max(max_disparity, disp_err)
        assert disp_err < 4e-11, (k, 'disparity', disparity, ref_disp, disp_err)
        assert abs(disparity - np.sum((x - y)**2)) < 3e-13, (k, 'disparity identity')
        assert np.max(np.abs(x-ref_x)) < 3e-12, (k, 'reference point matrix')
        assert abs(np.linalg.norm(x) - 1) < 3e-12, (k, 'unit trace')
        assert np.max(np.abs(np.sum(x, axis=0))) < 3e-12, (k, 'centered')
        if abs(float(ref_scale)) > 0:
            scale_err = abs(scale/ref_scale - 1)
        else:
            scale_err = abs(scale - ref_scale)
        assert scale_err < 3e-11, (k, 'nuclear norm', scale, ref_scale)
        assert np.max(abs(r.T @ r - np.eye(d))) < 3e-11, (k, 'orthogonal factor')
        # Avoid overflow in the objective by normalizing input matrices.
        na = np.linalg.norm(a/np.max(abs(a)))
        nb = np.linalg.norm(b/np.max(abs(b)))
        ua = a / np.max(abs(a)) / na
        ub = b / np.max(abs(b)) / nb
        residual = np.linalg.norm(ua @ r - ub)**2
        ref_residual = np.linalg.norm(ua @ ref_r - ub)**2
        assert abs(residual - ref_residual) < 3e-11, (k, 'optimality', residual, ref_residual)
        # Full-rank solutions are unique: check the actual transformed points.
        normalized = b/np.max(abs(b))
        normalized -= np.mean(normalized, axis=0)
        normalized /= np.linalg.norm(normalized)
        x_centered = a/np.max(abs(a))
        x_centered -= np.mean(x_centered, axis=0)
        x_centered /= np.linalg.norm(x_centered)
        singular = np.linalg.svd(normalized.T @ x_centered, compute_uv=False)
        if singular[-1] > singular[0]*1e-8:
            matrix_err = float(np.max(abs(y-ref_y)))
            max_matrix = max(max_matrix, matrix_err)
            assert matrix_err < 5e-10, (k, 'transformed coordinates', matrix_err)
        raw_singular = np.linalg.svd((a_unit/norm_a).T @ (b_unit/norm_b),
                                     compute_uv=False)
        if raw_singular[-1] > raw_singular[0]*1e-8:
            assert np.max(abs(r - ref_r)) < 5e-10, (k, 'rotation matrix')
    assert not list(output), 'Unexpected extra driver output'
    return len(cases), max_disparity, max_matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=81723)
    parser.add_argument("--count", type=int, default=350)
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    options = parser.parse_args()
    cases = make_cases(options.seed, options.count)
    modes = []
    if options.mode != "optimized":
        modes.append(("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]))
    if options.mode != "checked":
        modes.append(("optimized", ["-O3"]))
    for label, flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-procrustes-") as directory:
            exe = compile_driver(Path(directory), flags)
            count, disp_err, matrix_err = run_cases(exe, cases)
            print(f"{label}: PASS {count} SciPy Procrustes comparisons, "
                  f"maximum disparity error {disp_err:.3g}, "
                  f"max full-rank coordinate error {matrix_err:.3g}", flush=True)


if __name__ == "__main__":
    main()
