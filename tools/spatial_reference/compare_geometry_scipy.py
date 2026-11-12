#!/usr/bin/env python3
"""Compare SciFort Hausdorff distances and geometric_slerp with SciPy 1.17+.

Usage::
    python tools/spatial_reference/compare_geometry_scipy.py
    python tools/spatial_reference/compare_geometry_scipy.py --seed 84321
    python tools/spatial_reference/compare_geometry_scipy.py --mode checked

SciPy is used as an independent numerical oracle, not an implementation source.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.spatial import geometric_slerp
from scipy.spatial.distance import cdist, directed_hausdorff

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("kinds", "interpolate_kdtree", "spatial_distance", "spatial_hausdorff",
           "spatial_slerp", "spatial_procrustes", "spatial_voronoi", "spatial_convex_hull", "spatial_halfspace", "spatial")


def compile_drivers(tmp: Path, flags: list[str]) -> tuple[Path, Path]:
    command = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    objects = []
    for name in SOURCES:
        obj = tmp / (name + ".o")
        subprocess.run([*command, "-I", str(tmp), "-J", str(tmp), "-c",
                        str(ROOT / "src" / (name + ".f90")), "-o", str(obj)],
                       cwd=tmp, check=True, capture_output=True)
        objects.append(str(obj))
    executables = []
    for name in ("hausdorff_driver", "slerp_driver"):
        exe = tmp / name
        subprocess.run([*command, "-I", str(tmp),
                        str(ROOT / "tools/spatial_reference" / (name + ".f90")),
                        *objects, "-o", str(exe)], cwd=tmp,
                       check=True, capture_output=True)
        executables.append(exe)
    return executables[0], executables[1]


def make_hausdorff_cases(seed: int) -> list[tuple[np.ndarray, np.ndarray, float]]:
    rng = np.random.default_rng(seed)
    cases = []
    scales = (1.0, 1.0e-200, 1.0e200, 1.0e-70, 1.0e70)
    for scenario in range(420):
        dim = (1, 2, 3, 4, 8, 13)[scenario % 6]
        n, m = int(rng.integers(1, 16)), int(rng.integers(1, 17))
        scale = scales[(scenario // 6) % len(scales)]
        a, b = rng.normal(size=(n, dim)), rng.normal(size=(m, dim))
        mod = (scenario // 5) % 9
        if mod == 0:
            a[:]=0.0
        elif mod == 1:
            b[:]=0.0
        elif mod == 2:
            a[:]=b[:1]
        elif mod == 3:
            b[:]=a[:1]
        elif mod == 4:
            a = np.round(a)
            b = np.round(b)
        elif mod == 5:
            a[::2] = a[:1]
        elif mod == 6:
            b[::2] = b[:1]
        elif mod == 7:
            a = np.abs(a)
            b = np.abs(b)
        else:
            a[0]=100.0
        cases.append((a*scale,b*scale,scale))
    return cases


def run_hausdorff(exe: Path, cases: list[tuple]) -> tuple[int, float]:
    lines = [str(len(cases))]
    for a, b, _scale in cases:
        lines.append(f"{len(a)} {len(b)} {a.shape[1]}")
        for row in np.vstack((a,b)):
            lines.append(" ".join(f"{v:.17e}" for v in row))
    out = subprocess.run([str(exe)], input="\n".join(lines)+"\n", text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    records = out.stdout.splitlines()
    assert len(records) == 4*len(cases), (len(records),len(cases))
    checked, largest = 0, 0.0
    for k, (a, b, scale) in enumerate(cases):
        ab = list(map(float,records[4*k].split()))
        ba = list(map(float,records[4*k+1].split()))
        sym = list(map(float,records[4*k+2].split()))
        full = list(map(float,records[4*k+3].split()))
        assert all(int(r[-1]) == 0 for r in (ab,ba,sym,full)), (k, ab, ba, sym, full)
        ar, br = a/scale, b/scale
        scipy_ab = directed_hausdorff(ar,br,rng=0)[0]
        scipy_ba = directed_hausdorff(br,ar,rng=0)[0]
        pairwise = cdist(ar,br)
        expected_ab = pairwise.min(axis=1).max()
        expected_ba = pairwise.min(axis=0).max()
        assert abs(scipy_ab-expected_ab)<1e-12*max(1,expected_ab), (k,scipy_ab,expected_ab)
        assert abs(scipy_ba-expected_ba)<1e-12*max(1,expected_ba), (k,scipy_ba,expected_ba)
        for got, expected in ((ab[0]/scale,scipy_ab),
                              (ba[0]/scale,scipy_ba),
                              (sym[0]/scale,max(scipy_ab,scipy_ba)),
                              (full[0]/scale,scipy_ab)):
            err = abs(got-expected)/max(1.0,abs(expected))
            if err > 8e-13:
                raise AssertionError(f"Hausdorff case {k}: got={got} expected={expected}")
            largest = max(largest,err)
            checked += 1
        for record, dist_rows in ((ab,pairwise),(ba,pairwise.T),(full,pairwise)):
            i, j = int(record[1])-1, int(record[2])-1
            assert 0<=i<len(dist_rows) and 0<=j<len(dist_rows[0]), (k,record)
            ref_nearest = dist_rows.min(axis=1)
            assert abs(dist_rows[i,j]-ref_nearest[i]) <= 8e-13*max(1.0,ref_nearest[i]), (k,record)
            assert abs(ref_nearest[i]-ref_nearest.max()) <= 8e-13*max(1.0,ref_nearest.max()), (k,record)
        # The two independent implementations of the directed search should
        # select exactly the same deterministic first witness.
        assert ab[1:]==full[1:], (k,ab,full)
        iu, jv = int(sym[1])-1, int(sym[2])-1
        assert 0<=iu<len(a) and 0<=jv<len(b), (k,sym)
        assert abs(pairwise[iu,jv]-sym[0]/scale) <= 8e-13*max(1.0,sym[0]/scale)
        if int(sym[3]) == 1:
            assert ab[1:3]==sym[1:3], (k,ab,sym)
        else:
            assert ba[1]==sym[2] and ba[2]==sym[1], (k,ba,sym)
    return checked, largest


def make_slerp_cases(seed: int) -> list[tuple]:
    rng=np.random.default_rng(seed+722)
    cases=[]
    for scenario in range(360):
        dim=(2,3,4,5,9,16)[scenario%6]
        n=int(rng.integers(1,13))
        a=rng.normal(size=dim)
        a/=np.linalg.norm(a)
        b=rng.normal(size=dim)
        # Some deliberately narrow angles, but avoid near-antipodal ambiguity.
        if scenario%17==0:
            b=a+rng.normal(size=dim)*1e-4
        elif scenario%17==1:
            b=a.copy()
        elif scenario%17==2:
            direction=rng.normal(size=dim)
            direction-=np.dot(a,direction)*a
            direction/=np.linalg.norm(direction)
            b=-a+0.03*direction
        b/=np.linalg.norm(b)
        times=rng.uniform(-1.5,2.5,n)
        if n>1: times[0],times[-1]=0.0,1.0
        cases.append((a,b,times))
    return cases


def run_slerp(exe: Path, cases: list[tuple]) -> tuple[int,float]:
    lines=[str(len(cases))]
    for a,b,times in cases:
        lines.append(f"{len(a)} {len(times)}")
        lines.append(" ".join(f"{v:.17e}" for v in a))
        lines.append(" ".join(f"{v:.17e}" for v in b))
        lines.append(" ".join(f"{v:.17e}" for v in times))
    proc=subprocess.run([str(exe)],input="\n".join(lines)+"\n",text=True,
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
    records=proc.stdout.splitlines()
    at=0
    total,max_err=0,0.0
    for k,(a,b,times) in enumerate(cases):
        codes=records[at].split()
        at+=1
        assert codes==["0","0"], (k,codes)
        n=len(times)
        actual=np.array([list(map(float,records[at+j].split())) for j in range(n)])
        at+=n
        scalar=np.array(list(map(float,records[at].split())))
        at+=1
        ref=np.asarray(geometric_slerp(a,b,times),dtype=float)
        assert ref.shape==actual.shape, (k,ref.shape,actual.shape)
        assert np.max(np.abs(scalar-actual[0])) < 2e-14, (k,scalar,actual[0])
        err=np.max(np.abs(ref-actual))
        if err>2e-11:
            raise AssertionError(f"slerp case {k}: max error {err}; a={a} b={b}")
        norm_err = np.max(np.abs(np.linalg.norm(actual,axis=1)-1.0))
        assert norm_err<1e-11, (k,norm_err)
        max_err=max(max_err,float(err))
        total+=n
    assert at==len(records), (at,len(records))
    return total,max_err


def main()->None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--seed",type=int,default=20261102)
    parser.add_argument("--mode",choices=("both","checked","optimized"),default="both")
    args=parser.parse_args()
    cases_h=make_hausdorff_cases(args.seed)
    cases_s=make_slerp_cases(args.seed)
    modes=[]
    if args.mode!="optimized": modes.append(("checked",["-O0","-g","-fcheck=all","-fbacktrace"]))
    if args.mode!="checked": modes.append(("optimized",["-O3"]))
    for name, flags in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-spatial-geometry-") as tmp:
            h_exe,s_exe=compile_drivers(Path(tmp),flags)
            count_h,error_h=run_hausdorff(h_exe,cases_h)
            count_s,error_s=run_slerp(s_exe,cases_s)
            print(f"{name}: {len(cases_h)} Hausdorff scenarios ({count_h} values), "
                  f"max relative error={error_h:.3g}",flush=True)
            print(f"{name}: {len(cases_s)} geometric_slerp scenarios ({count_s} points), "
                  f"max absolute error={error_s:.3g}",flush=True)


if __name__ == "__main__":
    main()
