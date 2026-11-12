#!/usr/bin/env python3
"""Verify banded automatic full-spectrum dispatch against SciPy and explicit methods."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_banded.f90",
           "src/linalg_hessenberg.f90", "src/linalg_banded_tridiagonal.f90",
           "src/linalg_banded_eigen.f90"]


def compile_driver(work: Path, mode: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / f"banded_auto_{mode}"
    command = ["gfortran", "-std=f2018", *flags,
               *[str(ROOT / src) for src in SOURCES],
               str(ROOT / "tools/linalg_reference/driver_banded_auto.f90"),
               "-o", str(exe)]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"Fortran compilation failed:\n{result.stderr[-4000:]}")
    return exe


def generate_cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    shapes = [(96, 3), (160, 4), (191, 1), (192, 2), (192, 8),
              (192, 12), (256, 0), (256, 3), (256, 8), (256, 20),
              (320, 1), (320, 4), (320, 12), (320, 32)]
    for it in range(count):
        n, kd = shapes[it % len(shapes)]
        typ = 1 + (it // len(shapes) + it) % 4
        complex_case = typ >= 3
        lower = typ in (1, 3)
        a = np.zeros((n, n), dtype=complex if complex_case else float)
        np.fill_diagonal(a, rng.normal(size=n))
        for k in range(1, kd+1):
            entries = 0.14*rng.normal(size=n-k)
            if complex_case:
                entries = entries + 0.14j*rng.normal(size=n-k)
            ii = np.arange(n-k)
            a[ii+k, ii] = entries
            a[ii, ii+k] = np.conj(entries)
        family = (it // 4) % 8
        if family == 1:
            a += 4*np.eye(n)
        elif family == 2:
            a *= 1e90
        elif family == 3:
            a *= 1e-90
        elif family == 4:
            a = np.diag(np.ones(n)).astype(a.dtype)
        elif family == 5:
            a = np.diag(np.linspace(-2., 2., n)).astype(a.dtype)
            if kd:
                idx = np.arange(n-1)
                a[idx+1, idx] = 1e-4
                a[idx, idx+1] = 1e-4
        elif family == 6 and kd:
            for j in range(n-1):
                if j % 8 == 7:
                    a[j+1:min(n,j+kd+1), j] = 0
                    a[j, j+1:min(n,j+kd+1)] = 0
        elif family == 7:
            a *= 1e-140
        ab = np.zeros((kd+1, n), dtype=a.dtype)
        for j in range(n):
            for k in range(min(kd,n-j-1)+1):
                if lower:
                    ab[k,j] = a[j+k,j]
                else:
                    ab[kd-k,j+k] = a[j,j+k]
        yield typ, n, kd, a, ab


def check(exe: Path, cases, tol: float):
    text = [str(len(cases))]
    for typ, n, kd, a, ab in cases:
        text.append(f"{typ} {n} {kd}")
        text.extend(f"{v.real:.17e} {v.imag:.17e}" for v in ab.ravel(order="F"))
    result = subprocess.run([str(exe)], input="\n".join(text)+"\n", capture_output=True,
                            text=True, timeout=180,
                            env={**os.environ, "OPENBLAS_NUM_THREADS":"1", "OMP_NUM_THREADS":"1"})
    if result.returncode:
        raise RuntimeError(f"Fortran execution failed: {result.stderr[-3500:]}")
    lines = iter(result.stdout.splitlines())
    worst = 0.
    selected_bulge = selected_dense = 0
    for it, (typ, n, kd, a, ab) in enumerate(cases):
        statuses = [int(x) for x in next(lines).split()]
        if statuses != [0,0,0,0]:
            raise AssertionError(f"case {it} n={n} kd={kd}: statuses {statuses}")
        outputs = np.array([[float(v) for v in next(lines).split()] for _ in range(n)])
        wa, wb, wd, we = outputs.T
        use_bulge = n >= 192 and kd <= min(8, max(1, n//16))
        chosen = wb if use_bulge else wd
        if use_bulge:
            selected_bulge += 1
        else:
            selected_dense += 1
        # The automatic result should be bit-for-bit identical to the chosen
        # explicit implementation unless the documented fallback was used.
        if not np.array_equal(wa, chosen) or not np.array_equal(we, chosen):
            raise AssertionError(f"case {it}: automatic dispatch did not match chosen method")
        wref = linalg.eig_banded(ab, lower=typ in (1,3), eigvals_only=True)
        scale = max(float(np.max(np.abs(a))), 1e-300)
        normalized = max(1., float(linalg.norm(a/scale, 'fro')))
        err = float(np.max(np.abs((wa-wref)/scale)))/normalized
        worst = max(worst,err)
        if not np.isfinite(err) or err > tol:
            raise AssertionError(f"case {it} n={n} kd={kd}, error={err:.3e} > {tol}")
    return worst, selected_bulge, selected_dense


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["checked", "optimized", "both"], default="both")
    ap.add_argument("--cases", type=int, default=56)
    ap.add_argument("--seed", type=int, default=116)
    ap.add_argument("--tolerance", type=float, default=1e-8)
    args = ap.parse_args()
    cases = list(generate_cases(args.cases, args.seed))
    modes = ("checked", "optimized") if args.mode == "both" else (args.mode,)
    with tempfile.TemporaryDirectory(prefix="scifort_banded_auto_") as dirname:
        for mode in modes:
            exe = compile_driver(Path(dirname), mode)
            worst, bulge, dense = check(exe, cases, args.tolerance)
            print(f"{mode}: {len(cases)} SciPy comparisons passed; "
                  f"auto bulge={bulge}, dense={dense}, max scaled error={worst:.3e}")


if __name__ == "__main__":
    main()
