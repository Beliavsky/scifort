#!/usr/bin/env python3
"""Numerical SciPy 1.17 reference comparisons for rank-two SE(3) grids.

Compiles the standalone Fortran test driver without external Fortran libraries;
checks two-dimensional broadcasting, SE(3) representations and weighted means.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial.transform import RigidTransform

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    "kinds", "spatial_transform", "spatial_rotation_stack", "spatial_rotation_grid",
    "spatial_rigid_transform", "spatial_rigid_transform_grid",
)


def compile_driver(mode: str, out_dir: Path) -> Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    fc = ["gfortran", "-std=f2018", "-ffree-line-length-none", *flags]
    objects = []
    for source in SOURCES:
        obj = out_dir / (source + ".o")
        subprocess.run([*fc, "-I", str(out_dir), "-J", str(out_dir), "-c",
                        str(ROOT / "src" / (source + ".f90")), "-o", str(obj)],
                       check=True, cwd=out_dir, capture_output=True, text=True)
        objects.append(str(obj))
    exe = out_dir / "rigid_grid_driver"
    subprocess.run([*fc, "-I", str(out_dir),
                    str(ROOT / "tools/spatial_reference/rigid_grid_driver.f90"),
                    *objects, "-o", str(exe)],
                   check=True, cwd=out_dir, capture_output=True, text=True)
    return exe


def emit_rows(arr: np.ndarray) -> list[str]:
    """Write Fortran's i-fastest batch indexing, then the trailing vector."""
    return [" ".join(f"{v:.17e}" for v in arr[i, j])
            for j in range(arr.shape[1]) for i in range(arr.shape[0])]


def read_matrices(values: np.ndarray, offset: int, dims: tuple[int, int]) -> tuple[np.ndarray, int]:
    n, m = dims
    count = 16 * n * m
    data = values[offset:offset + count].reshape((m, n, 4, 4)).transpose(1, 0, 2, 3)
    return data, offset + count


def read_vectors(values: np.ndarray, offset: int, dims: tuple[int, int]) -> tuple[np.ndarray, int]:
    n, m = dims
    count = 3 * n * m
    data = values[offset:offset + count].reshape((m, n, 3)).transpose(1, 0, 2)
    return data, offset + count


def compare(a: np.ndarray, b: np.ndarray, label: str, state: dict) -> None:
    if a.shape != b.shape:
        raise AssertionError(f"{label}: shapes {a.shape} != {b.shape}")
    diff = np.max(np.abs(a - b))
    scale = max(1.0, float(np.max(np.abs(b))))
    relative = diff / scale
    state["max_rel"] = max(state["max_rel"], relative)
    state["values"] += a.size
    if not np.isfinite(relative) or relative > 5e-10:
        raise AssertionError(f"{label}: max relative error {relative:.3e}")


def run_case(exe: Path, rng: np.random.Generator, state: dict, index: int) -> None:
    n, m = rng.choice([1, 2, 3, 4], size=2)
    p = int(rng.choice([1, n]))
    q = int(rng.choice([1, m]))
    c = int(rng.choice([1, n]))
    d = int(rng.choice([1, m]))
    scale = float(rng.choice([1e-9, 0.03, 0.5, 2.0, 5.0]))
    twist = rng.normal(size=(n,m,6)) * scale
    rhs = rng.normal(size=(p,q,6)) * (0.75 if index % 3 else 1e-8)
    points = rng.normal(size=(c,d,3)) * rng.choice([0.1, 1.0, 10.0])
    weights = rng.uniform(0.05, 2.0, size=(n,m))
    if index % 7 == 0 and n > 1 and m > 1:
        weights[0,0] = 0.0
    if index % 23 == 0 and n * m > 1:
        weights[...] = 0.01
        weights[0,0] = 1.0
    lines = [f"{n} {m} {p} {q} {c} {d}"]
    lines += emit_rows(twist) + emit_rows(rhs) + emit_rows(points)
    lines += [f"{weights[i,j]:.17e}" for j in range(m) for i in range(n)]
    proc = subprocess.run([str(exe)], input="\n".join(lines)+"\n", text=True,
                          capture_output=True, check=False)
    if proc.returncode:
        raise AssertionError(f"Fortran case {index}: {proc.stderr[-1500:]}")
    values = np.fromstring(proc.stdout, sep=" ")
    off = 0
    output = []
    for dims in [(n,m), (max(n,p),max(m,q)), (max(n,c),max(m,d)),
                 (max(n,c),max(m,d)), (n,m), (n,m), (n,m), (n,m)]:
        if len(output) in (2,3):
            data, off = read_vectors(values,off,dims)
        else:
            data, off = read_matrices(values,off,dims)
        output.append(data)
    means = []
    for dims in [(1,1), (1,m), (1,n), (n,m)]:
        data, off = read_matrices(values,off,dims)
        means.append(data)
    if off != len(values):
        raise AssertionError(f"Fortran output length {len(values)} parsed {off}")
    tfa = RigidTransform.from_exp_coords(twist)
    tfb = RigidTransform.from_exp_coords(rhs)
    reference = [tfa.as_matrix(), (tfa*tfb).as_matrix(), tfa.apply(points),
                 tfa.apply(points,inverse=True), tfa.inv().as_matrix(), tfa.as_matrix(), tfa.as_matrix(), tfa.as_matrix()]
    for i,(out,ref) in enumerate(zip(output,reference)):
        compare(out,ref,f"case {index} operation {i}",state)
    compare(means[0][0,0],tfa.mean(weights=weights).as_matrix(),f"case {index} mean",state)
    compare(means[1][0],tfa.mean(axis=0,weights=weights).as_matrix(),f"case {index} axis0",state)
    compare(means[2][0],tfa.mean(axis=1,weights=weights).as_matrix(),f"case {index} axis1",state)
    # SciPy's generic RigidTransform.pow can differ from exp/log for fractional
    # powers; compare the mathematically specified exp/log operation instead.
    power_reference = RigidTransform.from_exp_coords(tfa.as_exp_coords()*0.5).as_matrix()
    compare(means[3],power_reference,f"case {index} power",state)
    state["cases"] += 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    parser.add_argument("--seed", type=int, default=20261103)
    parser.add_argument("--cases", type=int, default=140)
    args = parser.parse_args()
    modes = ["checked", "optimized"] if args.mode == "both" else [args.mode]
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix="scifort-rigid-grid-") as td:
            exe = compile_driver(mode, Path(td))
            rng = np.random.default_rng(args.seed)
            state = {"max_rel":0.0,"cases":0,"values":0}
            for i in range(args.cases):
                run_case(exe,rng,state,i)
            print(f"{mode}: PASS {state['cases']} cases, {state['values']} values, "
                  f"max scale-relative discrepancy {state['max_rel']:.3e}", flush=True)


if __name__ == "__main__":
    main()
