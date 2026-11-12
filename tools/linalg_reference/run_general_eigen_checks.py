#!/usr/bin/env python3
"""Compile/test SciFort's nonsymmetric eig/eigvals against SciPy.

Compare spectra as unordered multisets and check normalized left/right
residuals. Eigenvectors are not compared entrywise because phases, scaling,
and eigenspaces are nonunique. No LAPACK library is linked to the Fortran code.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_matrix_functions.f90",
           "src/linalg_general_eigen.f90"]


def generate_cases(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = sorted(set([0, 1, 2, 3, 4, 5, 8, 12, max_order]))
    for k in range(count):
        n = sizes[k % len(sizes)]
        cmplx = (k // len(sizes)) % 2 == 1
        a = rng.normal(size=(n, n))
        if cmplx:
            a = a + 1j * rng.normal(size=(n, n))
        family = (k // (2 * len(sizes))) % 13
        if family == 1:
            a = np.triu(a)
        elif family == 2:
            a = np.tril(a)
        elif family == 3:
            a = np.diag(np.diag(a))
        elif family == 4:
            a = (a + a.conj().T) / 2
        elif family == 5:
            a = np.triu(a, -1)
        elif family == 6:
            a = np.eye(n, dtype=a.dtype) * (1.5 + 0j if cmplx else 1.5)
        elif family == 7:
            a = np.triu(a, 1)
        elif family == 8:
            a *= 1e120
        elif family == 9:
            a *= 1e-120
        elif family == 10:
            a = a - a.conj().T
        elif family == 11:
            if n:
                a = np.eye(n, dtype=a.dtype)
                if n > 1:
                    a += 1e-6 * np.diag(np.ones(n - 1), 1)
        elif family == 12:
            if n > 1:
                a[-1, :] *= 1e-6
        yield np.asarray(a)


def compile_driver(work: Path, mode: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / ("general_eig_" + mode)
    cmd = ["gfortran", "-std=f2018", *flags, *(str(ROOT / p) for p in SOURCES),
           str(ROOT / "tools/linalg_reference/driver_general_eigen.f90"), "-o", str(exe)]
    proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(f"Fortran compile failed:\n{proc.stderr}")
    return exe


def run_native(work: Path, mode: str):
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / ("general_eig_native_" + mode)
    cmd = ["gfortran", "-std=f2018", *flags, *(str(ROOT / p) for p in SOURCES),
           str(ROOT / "test/test_linalg_general_eigen.f90"), "-o", str(exe)]
    proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(f"Native test compilation failed:\n{proc.stderr}")
    proc = subprocess.run([str(exe)], cwd=work, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(f"Native test execution failed:\n{proc.stdout}\n{proc.stderr}")
    print(f"{mode}: {proc.stdout.strip()}")


def compare(exe: Path, cases: list[np.ndarray], tol: float):
    inp = [str(len(cases))]
    for a in cases:
        inp.append(f"{2 if np.iscomplexobj(a) else 1} {a.shape[0]}")
        inp.extend(f"{z.real:.17e} {z.imag:.17e}" for z in a.ravel(order="F"))
    proc = subprocess.run([str(exe)], input="\n".join(inp) + "\n",
                          capture_output=True, text=True, timeout=180,
                          env={**os.environ, "OPENBLAS_NUM_THREADS":"1", "OMP_NUM_THREADS":"1"})
    if proc.returncode:
        raise RuntimeError(f"Fortran driver failed: {proc.stderr[-5000:]}")
    lines = iter(proc.stdout.splitlines())
    worst = {"eigenvalue": 0., "right": 0., "left": 0.}
    for k, a in enumerate(cases):
        n = a.shape[0]
        statuses = [int(s) for s in next(lines).split()]
        if statuses != [0, 0]:
            raise AssertionError(f"case {k}: n={n}, statuses={statuses}")
        ww = np.asarray([[float(x) for x in next(lines).split()] for _ in range(n)])
        eig = ww[:, 0] + 1j * ww[:, 1] if n else np.empty(0, dtype=complex)
        eigonly = ww[:, 2] + 1j * ww[:, 3] if n else np.empty(0, dtype=complex)
        vv = np.asarray([[float(x) for x in next(lines).split()] for _ in range(n*n)])
        vr = (vv[:, 0] + 1j*vv[:, 1]).reshape((n,n),order="F") if n else np.empty((0,0),complex)
        vl = (vv[:, 2] + 1j*vv[:, 3]).reshape((n,n),order="F") if n else np.empty((0,0),complex)
        if n == 0:
            continue
        ref = linalg.eigvals(a)
        scale = max(np.max(np.abs(a)), 1e-300)
        norm_a = max(float(linalg.norm(a / scale)), 1.)
        # SciPy's eigenvalue ordering and eigenvector phase are arbitrary.
        assignment = np.abs(eig[:,None]/scale - ref[None,:]/scale)
        row, col = linear_sum_assignment(assignment)
        spectral = max(float(np.max(assignment[row, col])),
                       float(np.max(np.abs((eig-eigonly)/scale)))) / norm_a
        right = float(np.max(np.abs((a/scale) @ vr - vr * (eig/scale)[None,:]))) / norm_a
        left = float(np.max(np.abs((a/scale).conj().T @ vl - vl*(eig/scale).conj()[None,:]))) / norm_a
        norms = np.r_[np.linalg.norm(vr,axis=0), np.linalg.norm(vl,axis=0)]
        if not np.all(np.isfinite(norms)) or np.max(np.abs(norms-1)) > 1e-10:
            raise AssertionError(f"case {k}: vectors not normalized")
        for key, value in [("eigenvalue",spectral), ("right",right), ("left",left)]:
            worst[key] = max(worst[key], value)
            if not np.isfinite(value) or value > tol:
                raise AssertionError(f"case {k}: n={n}, kind={'complex' if np.iscomplexobj(a) else 'real'}, "
                                     f"{key}={value:.3e} exceeds {tol:.1e}")
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["checked", "optimized", "both"], default="both")
    parser.add_argument("--cases", type=int, default=180)
    parser.add_argument("--max-order", type=int, default=20)
    parser.add_argument("--seed", type=int, default=8108)
    parser.add_argument("--tolerance", type=float, default=5e-8)
    args = parser.parse_args()
    cases = list(generate_cases(args.cases,args.max_order,args.seed))
    with tempfile.TemporaryDirectory(prefix="scifort_general_eig_") as temp:
        for mode in (["checked", "optimized"] if args.mode == "both" else [args.mode]):
            run_native(Path(temp),mode)
            exe = compile_driver(Path(temp),mode)
            worst = compare(exe,cases,args.tolerance)
            print(f"{mode}: {len(cases)} SciPy comparisons passed; worst=" + ", ".join(
                f"{name} {value:.3e}" for name,value in worst.items()))

if __name__ == "__main__":
    main()
