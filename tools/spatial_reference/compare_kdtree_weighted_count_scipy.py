#!/usr/bin/env python3
"""Compare SciFort's weighted/multi-radius KD-tree counts against SciPy.

Checked and optimized runs independently exercise vector cumulative counts,
noncumulative bins, arbitrary-order and repeated radii, scalar compatibility,
unit/one-sided/signed weights, periodic coordinates, and Minkowski p norms.
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
        cmd = [*fc, "-J", str(folder), "-I", str(folder), "-c",
               str(ROOT / "src" / f"{name}.f90"), "-o", str(obj)]
        subprocess.run(cmd, cwd=folder, capture_output=True, text=True, check=True)
        objects.append(str(obj))
    exe = folder / "kdtree_weighted_count_driver"
    cmd = [*fc, "-I", str(folder), str(ROOT / "tools/spatial_reference/kdtree_weighted_count_driver.f90"),
           *objects, "-o", str(exe)]
    subprocess.run(cmd, cwd=folder, capture_output=True, text=True, check=True)
    return exe


def run_case(exe: Path, xa: np.ndarray, xb: np.ndarray, box: np.ndarray,
             radii: np.ndarray, p: float, cumulative: bool, wa: np.ndarray,
             wb: np.ndarray, case: int) -> tuple[int, float]:
    n, dim = xa.shape
    m = xb.shape[0]
    numbers = [f"{n} {m} {dim} {len(radii)} {p if np.isfinite(p) else 'inf'} {int(cumulative)}"]
    for row in (box, radii, *xa, *xb, wa, wb):
        if len(row):
            numbers.append(" ".join(f"{v:.17g}" for v in row))
    proc = subprocess.run([str(exe)], input="\n".join(numbers)+"\n",
                          text=True, capture_output=True)
    if proc.returncode:
        raise RuntimeError(f"scenario {case}, stderr={proc.stderr[-3000:]}, stdout={proc.stdout[-2000:]}")
    stream = iter(proc.stdout.split())

    def arr(tag: str, dtype) -> np.ndarray:
        observed = next(stream)
        if observed != tag:
            raise AssertionError(f"scenario {case} expected {tag}, got {observed}")
        length = int(next(stream))
        return np.array([dtype(next(stream)) for _ in range(length)], dtype=dtype)

    got_count = arr("CI", int)
    got_weight = arr("WW", float)
    got_a = arr("WA", float)
    got_b = arr("WB", float)
    got_unit = arr("WU", float)
    a = cKDTree(xa, boxsize=box)
    b = cKDTree(xb, boxsize=box)
    ref_count = a.count_neighbors(b, radii, p=p, cumulative=cumulative)
    ref_weight = a.count_neighbors(b, radii, p=p, cumulative=cumulative, weights=(wa, wb))
    ref_a = a.count_neighbors(b, radii, p=p, cumulative=cumulative, weights=(wa, None))
    ref_b = a.count_neighbors(b, radii, p=p, cumulative=cumulative, weights=(None, wb))
    ref_unit = np.asarray(ref_count, dtype=float)
    if not np.array_equal(got_count, ref_count):
        raise AssertionError(f"scenario {case} integer counts: {got_count} vs {ref_count}")
    worst = 0.0
    for label, got, expected in (("weighted", got_weight, ref_weight), ("source", got_a, ref_a),
                                 ("destination", got_b, ref_b), ("unit", got_unit, ref_unit)):
        # Absolute tolerance covers cancellation in mixed-sign examples.
        if not np.allclose(got, expected, atol=3e-11, rtol=5e-12):
            raise AssertionError(f"scenario {case} {label} counts: {got} vs {expected}")
        if len(got):
            worst = max(worst, np.max(abs(got-expected)/np.maximum(1,abs(expected))))
    if len(radii):
        if next(stream) != "SC":
            raise AssertionError("missing scalar integer")
        got_scalar = int(next(stream))
        if next(stream) != "SW":
            raise AssertionError("missing scalar weighted")
        got_scalar_weight = float(next(stream))
        ref_scalar = int(a.count_neighbors(b, float(radii[0]), p=p))
        ref_scalar_weight = float(a.count_neighbors(b, float(radii[0]), p=p, weights=(wa,wb)))
        if got_scalar != ref_scalar:
            raise AssertionError(f"scenario {case} scalar integer: {got_scalar} vs {ref_scalar}")
        if not np.isclose(got_scalar_weight, ref_scalar_weight, atol=3e-11, rtol=5e-12):
            raise AssertionError(f"scenario {case} scalar weighted: {got_scalar_weight} vs {ref_scalar_weight}")
        worst = max(worst, abs(got_scalar_weight-ref_scalar_weight)/max(1,abs(ref_scalar_weight)))
    if next(stream,None) is not None:
        raise AssertionError("unexpected trailing output")
    return len(radii)*5+(2 if len(radii) else 0), worst


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=150)
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    args = parser.parse_args()
    modes = []
    if args.mode != "optimized":
        modes.append(("checked", ["-O0", "-g", "-fcheck=all", "-fbacktrace"]))
    if args.mode != "checked":
        modes.append(("optimized", ["-O3"]))
    for label, flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-count-neighbors-") as temp:
            exe = build(Path(temp), flags)
            rng = np.random.default_rng(args.seed)
            total = 0
            worst = 0.0
            for case in range(args.cases):
                dim = 1+case%7
                n = 3+(7*case)%55
                m = 2+(11*case)%49
                p = (1.0,1.5,2.0,3.0,np.inf)[case%5]
                box = rng.uniform(0.3,2.5,dim)
                box *= np.asarray([((case+j)%4 != 0) for j in range(dim)],dtype=float)
                if case%7==0:
                    box[:]=0
                if case%7==1:
                    box[:]=1.0
                scale = np.where(box>0,box,2.0)
                xa = rng.random((n,dim))*scale
                xb = rng.random((m,dim))*scale
                if n>1 and case%3==0:
                    xa[1]=xa[0]
                if case%3==0:
                    xb[0]=xa[0]
                if case%9==0:
                    xa=np.round(xa*3)/3
                    xb=np.round(xb*3)/3
                    for j in range(dim):
                        if box[j]>0:
                            xa[:,j] = np.minimum(xa[:,j],np.nextafter(box[j],0))
                            xb[:,j] = np.minimum(xb[:,j],np.nextafter(box[j],0))
                wa=rng.normal(0.2,1.2,n)
                wb=rng.normal(-0.1,0.8,m)
                if case%4==0:
                    wa[0]=0
                    wb[0]=0
                if case%17==0:
                    wa[:]=1.0
                nr=(0,1,2,5,8,12)[case%6]
                r=np.array([0.,0.12,0.6,0.6,1.1,6.0,0.31,2.0,0.0,0.31,0.9,3.0][:nr],dtype=float)
                cumulative=case%2==0
                if not cumulative:
                    r.sort()
                matches,error=run_case(exe,xa,xb,box,r,p,cumulative,wa,wb,case)
                total+=matches
                worst=max(worst,error)
            print(f"{label}: PASS {args.cases} randomized cases; {total} SciPy comparisons; "
                  f"max scale-relative discrepancy {worst:.3g}",flush=True)


if __name__ == "__main__":
    main()
