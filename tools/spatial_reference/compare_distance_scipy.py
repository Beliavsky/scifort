#!/usr/bin/env python3
"""Independent SciPy comparison for scifort_spatial pdist and cdist.

Run from any working directory::

    python tools/spatial_reference/compare_distance_scipy.py
    python tools/spatial_reference/compare_distance_scipy.py --seed 8675309

Requires NumPy, SciPy, and gfortran. Builds are isolated in temporary folders.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy.spatial.distance import cdist, pdist

ROOT = Path(__file__).resolve().parents[2]
METRICS = (
    "euclidean", "sqeuclidean", "cityblock", "chebyshev", "minkowski",
    "cosine", "correlation", "hamming", "jaccard", "canberra", "braycurtis",
)


def cases(seed: int) -> list[tuple[str, np.ndarray, np.ndarray, np.ndarray | None, float]]:
    rng = np.random.default_rng(seed)
    scenarios = []
    for metric in METRICS:
        for ndim in (1, 2, 3, 5, 8):
            for scenario in range(9):
                n = int(rng.integers(2, 9))
                m = int(rng.integers(1, 7))
                if scenario == 0:
                    xa, xb = rng.normal(size=(n, ndim)), rng.normal(size=(m, ndim))
                elif scenario == 1:
                    xa, xb = rng.integers(0, 2, size=(n, ndim)).astype(float), rng.integers(
                        0, 2, size=(m, ndim)).astype(float)
                elif scenario == 2:
                    xa, xb = rng.integers(-2, 3, size=(n, ndim)).astype(float), rng.integers(
                        -2, 3, size=(m, ndim)).astype(float)
                elif scenario == 3:
                    xa, xb = rng.uniform(0, 3, (n, ndim)), rng.uniform(0, 3, (m, ndim))
                elif scenario == 4:
                    xa, xb = rng.normal(size=(n, ndim)), rng.normal(size=(m, ndim))
                    xa[0, :] = 0.0
                    xb[0, :] = 0.0
                elif scenario == 5:
                    xa, xb = rng.normal(size=(n, ndim)), rng.normal(size=(m, ndim))
                    xa[0, :] = 1.0
                    xb[0, :] = 3.0
                    if metric == "correlation" and ndim > 1:
                        xa[0, 1] = 1.5
                        xb[0, 1] = 3.5
                elif scenario == 6:
                    xa, xb = rng.uniform(-10, 10, (n, ndim)), rng.uniform(-10, 10, (m, ndim))
                    xa[0, :] = xb[0, :]
                elif scenario == 7:
                    scale = 1e-3 if metric in ("cosine", "correlation") else 1e-90
                    xa, xb = rng.normal(size=(n, ndim))*scale, rng.normal(size=(m, ndim))*scale
                else:
                    scale = 1e3 if metric in ("cosine", "correlation") else 1e90
                    xa, xb = rng.normal(size=(n, ndim))*scale, rng.normal(size=(m, ndim))*scale
                w = None
                if scenario % 3:
                    w = rng.uniform(0.01, 3.0, size=ndim)
                    if scenario % 2 and ndim > 1:
                        w[0] = 0.0
                p = (0.75, 1.2, 2.0, 3.3, 6.0)[scenario % 5]
                scenarios.append((metric, xa, xb, w, p))
    return scenarios


def run_fortran(executable: Path, inputs: list[tuple]) -> list[str]:
    lines = [str(len(inputs))]
    for metric, xa, xb, w, p in inputs:
        lines.append(f"{len(xa)} {len(xb)} {xa.shape[1]} {metric} {int(w is not None)} {p:.17g}")
        if w is not None:
            lines.append(" ".join(f"{a:.17e}" for a in w))
        for row in np.concatenate((xa, xb)):
            lines.append(" ".join(f"{a:.17e}" for a in row))
    proc = subprocess.run([str(executable)], input="\n".join(lines)+"\n", text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return proc.stdout.splitlines()


def compare(output: list[str], inputs: list[tuple]) -> tuple[int, float]:
    offset = 0
    count = 0
    max_error = 0.0
    for case_id, (metric, xa, xb, w, p) in enumerate(inputs):
        status = output[offset].split()
        offset += 1
        if status != ["0", "0"]:
            raise AssertionError(f"case {case_id} metric={metric}: unexpected status {status}")
        n, m = len(xa), len(xb)
        size_a, size_b = n*m, n*(n-1)//2
        got_c = np.array([float(s) for s in output[offset:offset+size_a]]).reshape(n, m)
        offset += size_a
        got_p = np.array([float(s) for s in output[offset:offset+size_b]])
        offset += size_b
        kw = {}
        if w is not None:
            kw["w"] = w
        if metric == "minkowski":
            kw["p"] = p
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore", RuntimeWarning)
            ref_c = cdist(xa, xb, metric=metric, **kw)
            ref_p = pdist(xa, metric=metric, **kw)
        for label, got, ref in (("cdist", got_c, ref_c), ("pdist", got_p, ref_p)):
            if not np.array_equal(np.isnan(got), np.isnan(ref)):
                raise AssertionError(f"{label} NaN mismatch in case {case_id} ({metric})\n{got}\n{ref}")
            if not np.allclose(got, ref, rtol=5e-12, atol=2e-12, equal_nan=True):
                where = np.unravel_index(np.argmax(np.nan_to_num(abs(got-ref), nan=0)), got.shape)
                raise AssertionError(f"{label} error case {case_id} ({metric}) at {where}: "
                                     f"got={got[where]}, SciPy={ref[where]}, weights={w}, p={p}")
            finite = np.isfinite(got) & np.isfinite(ref)
            error = abs(got[finite]-ref[finite])/np.maximum(1.0,abs(ref[finite]))
            max_error = max(max_error, float(error.max(initial=0.0)))
            count += int(got.size)
    if offset != len(output):
        raise AssertionError(f"output not fully consumed: {len(output)-offset} lines")
    return count, max_error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--compiler", default=os.environ.get("FC", "gfortran"))
    args = parser.parse_args()
    inputs = cases(args.seed)
    for label, flags in (("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]),
                         ("optimized", ["-O3"])):
        with tempfile.TemporaryDirectory(prefix="scifort-spatial-") as temp:
            path = Path(temp)
            exe = path / "distance_driver"
            sources = [ROOT / "src" / x for x in ("kinds.f90", "interpolate_kdtree.f90", "spatial_distance.f90", "spatial_hausdorff.f90", "spatial_slerp.f90", "spatial_procrustes.f90", "spatial_voronoi.f90", "spatial_convex_hull.f90", "spatial_halfspace.f90", "spatial.f90")]
            sources.append(ROOT / "tools/spatial_reference/distance_driver.f90")
            subprocess.run([args.compiler, "-std=f2018", *flags, "-J", str(path),
                            "-o", str(exe), *map(str, sources)], cwd=path, check=True)
            result = run_fortran(exe, inputs)
            count, err = compare(result, inputs)
            print(f"{label}: {len(inputs)} cases, {count} distance values, "
                  f"max scale-relative error={err:.3g}: PASS")


if __name__ == "__main__":
    main()
