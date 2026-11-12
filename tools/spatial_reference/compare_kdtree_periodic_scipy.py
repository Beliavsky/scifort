#!/usr/bin/env python3
"""Compare periodic and Minkowski SciFort KD-tree APIs against SciPy.

Includes periodic/mixed periodic axes, queries outside the fundamental box,
nonperiodic trees, P=1, 1.5, 2, 3, and infinity, duplicate sites, and
zero/dense/sparse radii.  Builds Fortran with checked and optimized flags.

    python tools/spatial_reference/compare_kdtree_periodic_scipy.py
    python tools/spatial_reference/compare_kdtree_periodic_scipy.py --seed 8675309
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
    objects = []
    for name in ("kinds", "interpolate_kdtree", "spatial_distance", "spatial_hausdorff", "spatial_slerp", "spatial_procrustes", "spatial_voronoi", "spatial_convex_hull", "spatial_halfspace", "spatial"):
        obj = folder / f"{name}.o"
        subprocess.run([*fc, "-J", str(folder), "-I", str(folder), "-c",
                        str(ROOT / "src" / f"{name}.f90"), "-o", str(obj)],
                       cwd=folder, capture_output=True, check=True)
        objects.append(str(obj))
    exe = folder / "kdtree_periodic_driver"
    subprocess.run([*fc, "-I", str(folder),
                    str(ROOT / "tools/spatial_reference/kdtree_periodic_driver.f90"),
                    *objects, "-o", str(exe)], cwd=folder, capture_output=True, check=True)
    return exe


def run_case(exe: Path, xa: np.ndarray, xb: np.ndarray, box: np.ndarray,
             q: np.ndarray, radius: float, p: float, case: int) -> tuple[int, float]:
    n, d = xa.shape
    m = len(xb)
    lines = [f"{n} {m} {d} {radius:.17g} {p if np.isfinite(p) else 'inf'}"]
    lines += [" ".join(f"{v:.17g}" for v in row) for row in (box, q)]
    lines += [" ".join(f"{v:.17g}" for v in row) for row in xa]
    lines += [" ".join(f"{v:.17g}" for v in row) for row in xb]
    proc = subprocess.run([str(exe)], input="\n".join(lines)+"\n", text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode:
        raise RuntimeError(f"case {case}: {proc.stderr}\n{proc.stdout}")
    tokens = iter(proc.stdout.split())

    def header(tag: str) -> int:
        got = next(tokens)
        assert got == tag, (case, tag, got)
        return int(next(tokens))

    first = cKDTree(xa, boxsize=box)
    second = cKDTree(xb, boxsize=box)
    worst = 0.0
    k = header("Q")
    kref = min(6, n)
    assert k == kref
    ids = np.empty(k, dtype=int)
    got_dist = np.empty(k, dtype=float)
    for j in range(k):
        ids[j] = int(next(tokens))-1
        got_dist[j] = float(next(tokens))
    ref_dist, ref_ids = first.query(q, k=k, p=p)
    ref_dist = np.atleast_1d(ref_dist)
    ref_ids = np.atleast_1d(ref_ids)
    if not np.allclose(got_dist, ref_dist, rtol=5e-13, atol=2e-14):
        raise AssertionError(f"case {case}: nearest distances {got_dist} vs {ref_dist}")
    all_dist, _ = first.query(q, k=n, p=p)
    for j in range(k):
        if (sum(abs(all_dist - ref_dist[j]) < 2e-13) == 1) and ids[j] != ref_ids[j]:
            raise AssertionError(f"case {case}: nearest index {ids[j]} vs {ref_ids[j]}")
    worst = max(worst, float(np.max(abs(got_dist-ref_dist)/np.maximum(1,abs(ref_dist)))))

    nb = header("B")
    ball = [int(next(tokens))-1 for _ in range(nb)]
    expected_ball = sorted(first.query_ball_point(q, radius, p=p))
    if ball != expected_ball:
        raise AssertionError(f"case {case}: query_ball_point {ball} vs {expected_ball}")

    npairs = header("P")
    pairs = [(int(next(tokens))-1, int(next(tokens))-1) for _ in range(npairs)]
    exp_pairs = sorted(first.query_pairs(radius, p=p))
    if pairs != exp_pairs:
        raise AssertionError(f"case {case}: query_pairs: {set(pairs)^set(exp_pairs)}")

    nrp = header("R")
    rp = [int(next(tokens)) for _ in range(nrp)]
    ncol = header("J")
    cols = [int(next(tokens))-1 for _ in range(ncol)]
    if nrp != n+1 or rp[0] != 1 or rp[-1] != ncol+1:
        raise AssertionError(f"case {case}: malformed CSR pointers")
    expect_cross = first.query_ball_tree(second, radius, p=p)
    for row, expected in enumerate(expect_cross):
        if cols[rp[row]-1:rp[row+1]-1] != sorted(expected):
            raise AssertionError(f"case {case}: CSR row {row} mismatch")

    got_count = header("C")
    expect_count = first.count_neighbors(second, radius, p=p)
    if got_count != int(expect_count):
        raise AssertionError(f"case {case}: count_neighbors {got_count} vs {expect_count}")

    nentries = header("S")
    sparse = [(int(next(tokens))-1, int(next(tokens))-1, float(next(tokens)))
              for _ in range(nentries)]
    expected = first.sparse_distance_matrix(second, radius, p=p, output_type="dict")
    if len(sparse) != len(expected) or [(i,j) for i,j,_ in sparse] != sorted(expected):
        raise AssertionError(f"case {case}: sparse pair positions mismatch")
    for i,j,value in sparse:
        ref = expected[i,j]
        err = abs(value-ref)/max(1.0,abs(ref))
        if err > 5e-13:
            raise AssertionError(f"case {case}: sparse {(i,j)} distance {value} vs {ref}")
        worst = max(worst,err)
    if next(tokens,None) is not None:
        raise AssertionError(f"case {case}: extra output")
    return k+nb+npairs+ncol+nentries, worst


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--cases", type=int, default=150)
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    args = parser.parse_args()
    modes = []
    if args.mode != "optimized":
        modes.append(("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]))
    if args.mode != "checked":
        modes.append(("optimized", ["-O3"]))
    for label, flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-kdtree-periodic-") as tmp:
            exe = build(Path(tmp), flags)
            rng = np.random.default_rng(args.seed)
            total, worst = 0, 0.0
            for scenario in range(args.cases):
                dim = 1+scenario%7
                n = 4+(scenario*13)%48
                m = 3+(scenario*7)%34
                p = (1.0, 1.5, 2.0, 3.0, np.inf)[scenario%5]
                # Vary scales, all/none/mixed periodic axes, and box widths.
                box = rng.uniform(0.3, 5.0, size=dim)
                periodic_mask = np.array([((i+scenario)%4 != 0) for i in range(dim)])
                if scenario%7 == 0:
                    periodic_mask[:] = False
                elif scenario%7 == 1:
                    periodic_mask[:] = True
                box *= periodic_mask
                xa = rng.random((n,dim))
                xb = rng.random((m,dim))
                widths = np.where(box>0,box,2.5)
                xa *= widths
                xb *= widths
                if scenario%3 == 0:
                    xa[1] = xa[0]
                    xb[0] = xa[0]
                if scenario%9 == 0:
                    xa = np.round(xa*3)/3
                    xb = np.round(xb*3)/3
                    for j in range(dim):
                        if box[j] > 0:
                            xa[:,j] = np.minimum(xa[:,j],np.nextafter(box[j],0))
                            xb[:,j] = np.minimum(xb[:,j],np.nextafter(box[j],0))
                q = rng.random(dim)*widths*3 - widths
                r = (0., 0.17, 0.6, 2., 20.)[scenario%5]
                matches, err = run_case(exe, xa, xb, box, q, r, p, scenario)
                total += matches
                worst = max(worst,err)
            print(f"{label}: PASS {args.cases} periodic/Minkowski cases; {total} "
                  f"matched index/distance records; max relative discrepancy {worst:.3g}", flush=True)


if __name__ == "__main__":
    main()
