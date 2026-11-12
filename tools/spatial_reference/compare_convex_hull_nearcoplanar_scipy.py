#!/usr/bin/env python3
"""Compare difficult nearly-coplanar 3D hulls with SciPy's Qhull.

For points within a few rounding units of a supporting plane, exact vertex
indices are not numerically well determined. Compare *geometry* instead:
complete facet adjacency, outward containment, surface area and volume.

Example:
    python tools/spatial_reference/compare_convex_hull_nearcoplanar_scipy.py --mode both
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial import ConvexHull, QhullError

from compare_convex_hull_scipy import compile_driver


def make_cases(seed: int, count: int) -> list[tuple[np.ndarray, ConvexHull]]:
    rng = np.random.default_rng(seed)
    corners = np.array([[x, y, z] for x in (-1.0, 1.0)
                        for y in (-1.0, 1.0) for z in (-1.0, 1.0)])
    cases: list[tuple[np.ndarray, ConvexHull]] = []
    for k in range(count):
        n = int(rng.integers(18, 70))
        points = rng.uniform(-1, 1, size=(n, 3))
        faces = rng.integers(0, 6, n)
        axes = faces // 2
        signs = np.where(faces % 2 == 0, -1.0, 1.0)
        jitter = (0.0, 1e-14, 1e-12, 1e-10, 1e-8)[k % 5]
        points[np.arange(n), axes] = signs + rng.normal(size=n) * jitter
        points = np.vstack([points, corners])
        try:
            reference = ConvexHull(points)
        except QhullError:
            continue
        cases.append((points, reference))
    return cases


def check_cases(executable: Path, cases: list[tuple[np.ndarray, ConvexHull]]) -> tuple[float, float]:
    lines = [str(len(cases))]
    for points, _ in cases:
        lines.append(f"{len(points)} 3")
        lines.extend(" ".join(f"{v:.17e}" for v in point) for point in points)
    result = subprocess.run([str(executable)], input="\n".join(lines) + "\n",
                            capture_output=True, text=True, check=True)
    output = iter(result.stdout.splitlines())
    max_area = max_volume = 0.0
    for case_number, (points, ref) in enumerate(cases):
        status, nvertices, nfaces = map(int, next(output).split())
        if status != 0:
            raise AssertionError((case_number, "hull rejected", len(points)))
        area, volume = map(float, next(output).split())
        vertices = np.array([int(v)-1 for v in next(output).split()])
        if nfaces != 2*nvertices-4 or len(np.unique(vertices)) != nvertices:
            raise AssertionError((case_number, "Euler topology", nfaces, nvertices))
        area_error = abs(area-ref.area)/ref.area
        volume_error = abs(volume-ref.volume)/ref.volume
        max_area = max(max_area, area_error)
        max_volume = max(max_volume, volume_error)
        if area_error > 2e-9 or volume_error > 2e-9:
            raise AssertionError((case_number, "area/volume", area_error, volume_error))
        for f in range(nfaces):
            simplex = np.array([int(v)-1 for v in next(output).split()])
            equation = np.array([float(v) for v in next(output).split()])
            neighbors = np.array([int(v)-1 for v in next(output).split()])
            if np.any(simplex < 0) or np.any(simplex >= len(points)):
                raise AssertionError((case_number, "simplex indices", f))
            if np.any(neighbors < 0) or np.any(neighbors >= nfaces):
                raise AssertionError((case_number, "neighbors", f))
            # Allow roundoff-scale shifts in the nearly coplanar facets.
            if np.max(points @ equation[:3] + equation[3]) > 2e-10:
                raise AssertionError((case_number, "facet not supporting", f))
    return max_area, max_volume


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=512)
    parser.add_argument("--count", type=int, default=100)
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    args = parser.parse_args()
    cases = make_cases(args.seed, args.count)
    modes = []
    if args.mode != "optimized":
        modes.append(("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]))
    if args.mode != "checked":
        modes.append(("optimized", ["-O3"]))
    for label, flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-near-hull-") as td:
            executable = compile_driver(Path(td), flags)
            area_error, volume_error = check_cases(executable, cases)
            print(f"{label}: {len(cases)} near-coplanar 3D SciPy hull comparisons PASS; "
                  f"max relative area error {area_error:.3g}, "
                  f"volume error {volume_error:.3g}", flush=True)


if __name__ == "__main__":
    main()
