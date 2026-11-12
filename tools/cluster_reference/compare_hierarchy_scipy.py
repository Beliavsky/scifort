#!/usr/bin/env python3
"""Independent SciPy 1.17 oracle for native hierarchical clustering.

Builds a standalone Fortran driver (no FPM or external numerical libraries).
Tests merge topology/height, cophenetic distances, leaf ordering, flat cluster
partitions, and several linkage methods/metrics.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, cophenet, leaves_list
from scipy.spatial.distance import pdist

ROOT = Path(__file__).resolve().parents[2]
METHODS = ("single", "complete", "average", "weighted", "centroid", "median", "ward")


def compile_driver(directory: Path, mode: str) -> Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    command = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags, "-J", str(directory), "-I", str(directory)]
    objects = []
    for name in ("kinds", "spatial_distance", "cluster_hierarchy"):
        obj = directory / (name + ".o")
        subprocess.run([*command, "-c", str(ROOT / "src" / (name + ".f90")), "-o", str(obj)],
                       cwd=directory, check=True, capture_output=True, text=True)
        objects.append(str(obj))
    exe = directory / "cluster_driver"
    subprocess.run([*command, str(ROOT / "tools/cluster_reference/cluster_hierarchy_driver.f90"),
                    *objects, "-o", str(exe)], cwd=directory, check=True,
                   capture_output=True, text=True)
    return exe


def partitions_equivalent(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    return np.array_equal(a[:, None] == a[None, :], b[:, None] == b[None, :])


def run_case(exe: Path, x: np.ndarray, method: str, metric: str, cutoff: float, k: int):
    n, dim = x.shape
    s = f"{n} {dim} {method} {metric}\n{cutoff:.17e} {k}\n"
    s += "\n".join(" ".join(f"{v:.17e}" for v in row) for row in x) + "\n"
    proc = subprocess.run([str(exe)], input=s, text=True, capture_output=True, check=True)
    lines = proc.stdout.strip().splitlines()
    if lines[0].strip() == 'INVALID':
        raise AssertionError(f'Fortran rejected valid {method} {metric} input')
    z_f = np.array([[float(v) for v in line.split()] for line in lines[:n - 1]])
    cp = np.fromstring(lines[n - 1], sep=' ')
    dist_flat = np.fromstring(lines[n], sep=' ', dtype=int)
    max_flat = np.fromstring(lines[n + 1], sep=' ', dtype=int)
    leaves = np.fromstring(lines[n + 2], sep=' ', dtype=int)
    z_p = linkage(x, method=method, metric=metric)
    # Offsetting the SciPy node numbers by one yields the SciFort convention.
    expected_z = z_p.copy()
    expected_z[:, :2] += 1
    err = float(np.max(np.abs(z_f[:, 2:] - expected_z[:, 2:])))
    if not np.allclose(z_f[:, 2:], expected_z[:, 2:], rtol=6e-12, atol=1e-11):
        raise AssertionError(f'Linkage distance/count mismatch {method} {metric}: {z_f} vs {expected_z}')
    if not np.array_equal(z_f[:, :2], expected_z[:, :2]):
        # With exact equal distances multiple valid trees may exist.
        # The randomized tests avoid ties; when the algorithms still choose a
        # different topology, catch it rather than hiding it via sorted heights.
        raise AssertionError(f'Merge topology mismatch {method} {metric}: {z_f} vs {expected_z}')
    if not np.allclose(cp, cophenet(z_p), rtol=6e-12, atol=1e-11):
        raise AssertionError(f'Cophenetic distances mismatch {method} {metric}')
    if not np.array_equal(leaves, leaves_list(z_p) + 1):
        raise AssertionError(f'Leaf order mismatch {method} {metric}')
    if not partitions_equivalent(dist_flat, fcluster(z_p, cutoff, criterion='distance')):
        raise AssertionError(f'Distance partition mismatch {method} {metric}')
    # scipy.fcluster(... maxclust) has specialized small-n behavior for n<4,
    # so only compare directly on reasonably sized trees.
    if n >= 6 and k <= n - 2 and not partitions_equivalent(max_flat, fcluster(z_p, k, criterion='maxclust')):
        raise AssertionError(f'Maxclust partition mismatch {method} {metric}, k={k}')
    return err, n * (n - 1) // 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['checked', 'optimized', 'both'], default='both')
    ap.add_argument('--seed', type=int, default=20261103)
    ap.add_argument('--cases', type=int, default=140)
    args = ap.parse_args()
    modes = ('checked', 'optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        rng = np.random.default_rng(args.seed)
        max_err = 0.0
        pairs = 0
        with tempfile.TemporaryDirectory(prefix='scifort-hierarchy-') as temp:
            exe = compile_driver(Path(temp), mode)
            for c in range(args.cases):
                n = int(rng.integers(4, 17))
                p = int(rng.integers(1, 7))
                method = METHODS[c % len(METHODS)]
                metric = ('cityblock' if method in METHODS[:4] and c % 3 == 0 else 'euclidean')
                x = rng.normal(size=(n,p)) * rng.uniform(0.1, 4.0, size=p)
                if c % 13 == 0:
                    x += 1000.0
                z = linkage(x, method=method, metric=metric)
                cutoff = float(z[int(rng.integers(n-1)), 2])
                # Avoid cutoff exactly on a floating-point borderline to
                # prevent minute cross-language rounding differences.
                cutoff += 1.0e-9
                k = int(rng.integers(2, min(n - 1, 7)))
                err, amount = run_case(exe, x, method, metric, cutoff, k)
                max_err = max(max_err, err)
                pairs += amount
            print(f'{mode}: {args.cases} SciPy linkage / flat-cluster comparisons PASS; '
                  f'{pairs} cophenetic values; max height error {max_err:.3g}', flush=True)


if __name__ == '__main__':
    main()
