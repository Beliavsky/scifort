#!/usr/bin/env python3
"""Independent weighted Boolean pdist/cdist comparison against SciPy 1.17+.

Run::
    python tools/spatial_reference/compare_boolean_distance_scipy.py
    python tools/spatial_reference/compare_boolean_distance_scipy.py --seed 86123

Native executable is compiled in a temporary directory for each build mode.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import os
import subprocess
import tempfile
import warnings

import numpy as np
from scipy.spatial import distance

ROOT = Path(__file__).resolve().parents[2]
METRICS = ("dice", "rogerstanimoto", "russellrao", "sokalsneath", "yule")


def cases(seed: int) -> list[tuple]:
    rng = np.random.default_rng(seed)
    inputs = []
    for metric in METRICS:
        for ndim in (1, 2, 3, 6, 12, 31):
            for scenario in range(15):
                n, m = int(rng.integers(2, 10)), int(rng.integers(1, 8))
                if scenario in (0, 7, 8):
                    prob = {0: 0.5, 7: 0.05, 8: 0.95}[scenario]
                    xa = (rng.random((n, ndim)) < prob).astype(float)
                    xb = (rng.random((m, ndim)) < prob).astype(float)
                elif scenario == 1:
                    xa, xb = np.zeros((n, ndim)), np.zeros((m, ndim))
                elif scenario == 2:
                    xa, xb = np.ones((n, ndim)), np.ones((m, ndim))
                elif scenario == 3:
                    xa = (rng.random((n, ndim)) < 0.5).astype(float)
                    xb = 1.0-xa[:1].repeat(m, axis=0)
                elif scenario in (4, 12):
                    xa = rng.integers(0, 2, size=(n, ndim)).astype(float)
                    xb = xa[:1].repeat(m, axis=0)
                elif scenario == 5:
                    xa = rng.integers(-3, 4, size=(n, ndim)).astype(float)
                    xb = rng.integers(-3, 4, size=(m, ndim)).astype(float)
                elif scenario == 6:
                    xa, xb = rng.normal(size=(n, ndim)), rng.normal(size=(m, ndim))
                elif scenario == 9:
                    xa = np.ones((n, ndim))
                    xb = np.zeros((m, ndim))
                elif scenario == 10:
                    xa = np.zeros((n, ndim))
                    xb = np.ones((m, ndim))
                elif scenario == 11:
                    xa = np.eye(n, ndim)
                    xb = np.eye(m, ndim)
                elif scenario == 13:
                    xa, xb = rng.integers(-1, 2, size=(n, ndim)).astype(float), rng.integers(
                        -1, 2, size=(m, ndim)).astype(float)
                else:
                    xa, xb = rng.integers(0, 2, size=(n, ndim)).astype(float), rng.integers(
                        0, 2, size=(m, ndim)).astype(float)
                    xa[:, ::2] = 0.0
                    xb[:, ::2] = 0.0
                mode = (scenario + ndim) % 8
                w = None
                if mode == 1:
                    w = rng.uniform(0.001, 10.0, ndim)
                elif mode == 2:
                    w = np.zeros(ndim)
                elif mode == 3:
                    w = rng.integers(0, 5, size=ndim).astype(float)
                elif mode == 4:
                    w = np.zeros(ndim)
                    w[int(rng.integers(0, ndim))] = 1.0
                elif mode == 5:
                    w = rng.uniform(0.1, 10.0, ndim)*1.0e70
                elif mode == 6:
                    w = rng.uniform(0.1, 10.0, ndim)*1.0e-70
                elif mode == 7:
                    w = rng.uniform(0, 3.0, ndim)
                    w[rng.random(ndim) < 0.5] = 0.0
                inputs.append((metric, xa, xb, w))
    return inputs


def run_fortran(executable: Path, inputs: list[tuple]) -> list[str]:
    text = [str(len(inputs))]
    for metric, xa, xb, w in inputs:
        n, ndim = xa.shape
        text.append(f"{n} {len(xb)} {ndim} {metric} {int(w is not None)} 2.0")
        if w is not None:
            text.append(" ".join(f"{val:.17e}" for val in w))
        for row in np.concatenate((xa, xb)):
            text.append(" ".join(f"{val:.17e}" for val in row))
    proc = subprocess.run([str(executable)], input="\n".join(text)+"\n",
                          text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return proc.stdout.splitlines()


def compare(lines: list[str], inputs: list[tuple]) -> tuple[int, float]:
    offset, count, max_error = 0, 0, 0.0
    for index, (metric, xa, xb, w) in enumerate(inputs):
        codes = lines[offset].split()
        offset += 1
        if codes != ["0", "0"]:
            raise AssertionError(f"case {index} ({metric}) bad status {codes}")
        n, m = len(xa), len(xb)
        nc, npairs = n*m, n*(n-1)//2
        got_c = np.array([float(s) for s in lines[offset:offset+nc]]).reshape(n, m)
        offset += nc
        got_p = np.array([float(s) for s in lines[offset:offset+npairs]])
        offset += npairs
        keywords = {} if w is None else {"w": w}
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore", RuntimeWarning)
            # Force the documented Boolean domain. For numeric values other
            # than 0/1, some SciPy fast paths (notably Dice) evaluate their
            # formulas on raw numeric values instead of Boolean indicators.
            expected_c = distance.cdist(xa != 0.0, xb != 0.0, metric=metric, **keywords)
            expected_p = distance.pdist(xa != 0.0, metric=metric, **keywords)
        for which, actual, expected in (("cdist", got_c, expected_c), ("pdist", got_p, expected_p)):
            if not np.array_equal(np.isnan(actual), np.isnan(expected)):
                raise AssertionError(f"{which} case {index} ({metric}) NaN mismatch:\n{actual}\n{expected}")
            if not np.allclose(actual, expected, rtol=2.0e-12, atol=2.0e-13, equal_nan=True):
                raise AssertionError(f"{which} case {index} ({metric}):\n{actual}\n{expected}")
            mask = np.isfinite(actual) & np.isfinite(expected)
            if np.any(mask):
                max_error = max(max_error, float(np.max(
                    np.abs(actual[mask]-expected[mask])/np.maximum(1.0, np.abs(expected[mask])))))
            count += actual.size
    if offset != len(lines):
        raise AssertionError(f"unparsed output: {len(lines)-offset} lines")
    return count, max_error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    args = parser.parse_args()
    input_cases = cases(args.seed)
    modes = []
    if args.mode != "optimized":
        modes.append(("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]))
    if args.mode != "checked":
        modes.append(("optimized", ["-O3"]))
    for label, flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-spatial-boolean-") as tempdir:
            work = Path(tempdir)
            obj = []
            for name in ("kinds", "interpolate_kdtree", "spatial_distance", "spatial_hausdorff", "spatial_slerp", "spatial_procrustes", "spatial_voronoi", "spatial_convex_hull", "spatial_halfspace", "spatial"):
                output = work/f"{name}.o"
                subprocess.run([os.environ.get("FC", "gfortran"), "-std=f2018", *flags,
                                "-ffree-line-length-none", "-I", str(work), "-J", str(work),
                                "-c", str(ROOT/"src"/f"{name}.f90"), "-o", str(output)],
                               cwd=work, check=True, capture_output=True)
                obj.append(str(output))
            exe = work/"boolean_distance_driver"
            subprocess.run([os.environ.get("FC", "gfortran"), "-std=f2018", *flags,
                            "-ffree-line-length-none", "-I", str(work),
                            str(ROOT/"tools"/"spatial_reference"/"distance_driver.f90"),
                            *obj, "-o", str(exe)], cwd=work, check=True, capture_output=True)
            total, error = compare(run_fortran(exe, input_cases), input_cases)
            print(f"{label}: {len(input_cases)} cases, {total} distance values, "
                  f"maximum scale-relative discrepancy {error:.3e} PASS", flush=True)


if __name__ == "__main__":
    main()
