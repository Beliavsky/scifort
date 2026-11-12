#!/usr/bin/env python3
"""Test O(n*bandwidth) eigenvalue selection using banded Sturm counts.

This exercises real/complex upper/lower packed bands and compares selected
spectra to scipy.linalg.eigvalsh. Nonunique eigenvectors are not involved.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_banded.f90", "src/linalg_hessenberg.f90",
           "src/linalg_banded_tridiagonal.f90", "src/linalg_banded_eigen.f90"]


def compile_program(work: Path, mode: str, name: str):
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    executable = work / (name + "_banded_sturm_" + mode)
    sources = [str(ROOT / name) for name in SOURCES]
    program = ("test/test_linalg_banded_sturm.f90" if name == "native"
               else "tools/linalg_reference/driver_banded_sturm.f90")
    cmd = ["gfortran", "-std=f2018", *flags, *sources,
           str(ROOT / program),
           "-o", str(executable)]
    p = subprocess.run(cmd, cwd=work, capture_output=True, text=True, check=False)
    if p.returncode:
        raise RuntimeError("Fortran compilation failed:\n" + p.stderr)
    return executable


def make_cases(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = sorted({1, 2, 3, 4, 7, 12, 24, 48, max_order})
    for index in range(count):
        n = int(rng.choice(sizes))
        kd = int(rng.choice([0, 1, 2, 3, 5, max(0, n-1)]))
        kd = min(kd, n-1)
        select_mode = (index // 4) % 4
        if n > 48:
            kd = min(kd,4)
            if select_mode == 0:
                select_mode = 1
        kind = 1+index % 4  # real upper/lower, complex upper/lower
        lower = kind in (2, 4)
        complex_case = kind in (3, 4)
        a = np.zeros((n, n), dtype=complex if complex_case else float)
        regime = (index // 16) % 10
        d = rng.normal(size=n)
        if regime == 0:
            d[:] = 0
        elif regime == 1:
            d[:] = 2
        elif regime == 2:
            d *= 1e-150
        elif regime == 3:
            d *= 1e150
        elif regime == 4:
            d[:] = np.arange(n) - n/2
        elif regime == 5:
            d[:] = 5
        elif regime == 6:
            d[:] = 1e-10
        a[np.arange(n),np.arange(n)] = d
        for j in range(n):
            for k in range(1, min(kd, n-1-j)+1):
                if regime in (0, 1) or (regime == 5 and (j % 6) == 5):
                    val = 0.
                else:
                    val = rng.normal() + (1j * rng.normal() if complex_case else 0)
                    if regime == 2:
                        val *= 1e-150
                    elif regime == 3:
                        val *= 1e150
                    elif regime in (6,7):
                        val *= 1e-11
                a[j+k,j] = val
                a[j,j+k] = np.conj(val)
        # SciPy packed convention, with unused slots deliberately left zero.
        ab = np.zeros((kd+1,n),dtype=complex)
        for j in range(n):
            for k in range(min(kd,n-1-j)+1):
                if lower:
                    ab[k,j] = a[j+k,j]
                else:
                    ab[kd-k,j+k] = a[j,j+k]
        wr = linalg.eigvalsh(a,check_finite=True)
        start = 1
        stop = n
        low = -1.0
        high = 1.0
        if select_mode in (1,3):
            start = int(rng.integers(1,n+1))
            stop = (min(n,start+5) if n > 48
                    else int(rng.integers(start,n+1)))
            expected = wr[start-1:stop]
        elif select_mode == 2:
            if n == 1:
                low = float(wr[0]-max(abs(wr[0]),1.0))
                high = float(wr[0]+max(abs(wr[0]),1.0))
            else:
                loidx = int(rng.integers(0,n))
                hiidx = (min(n,loidx+6) if n > 48
                         else int(rng.integers(loidx+1,n+1)))
                span = max(1.0,float(np.max(abs(wr))))
                low = float(wr[loidx-1]+0.5*(wr[loidx]-wr[loidx-1])) if loidx else float(wr[0]-span)
                high = float(wr[hiidx-1]+0.5*(wr[hiidx]-wr[hiidx-1])) if hiidx<n else float(wr[-1]+span)
            if low >= high:
                low, high = float(wr[0]-span), float(wr[-1]+span)
            expected = wr[(wr > low) & (wr <= high)]
        else:
            expected = wr
        capacity = len(expected) + (2 if select_mode==2 else 0)
        yield (kind,n,kd,select_mode,start,stop,low,high,capacity,ab,wr,expected)


def check_cases(executable: Path, cases, tolerance: float):
    lines = [str(len(cases))]
    for kind,n,kd,mode,first,last,low,high,capacity,ab,wr,expected in cases:
        lines.append(f"{kind} {n} {kd} {mode} {first} {last} {low:.17e} {high:.17e} {capacity}")
        for v in ab.T.flat:
            lines.append(f"{v.real:.17e} {v.imag:.17e}")
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    proc = subprocess.run([str(executable)], input="\n".join(lines)+"\n", text=True,
                          capture_output=True, timeout=130, env=env)
    if proc.returncode:
        raise RuntimeError(f"Fortran driver exited {proc.returncode}:\n{proc.stderr}\n{proc.stdout[-2000:]}")
    parts = iter(proc.stdout.splitlines())
    worst = 0.
    for idx,(_,n,kd,mode,_,_,_,_,capacity,ab,wr,expected) in enumerate(cases):
        status, count = (int(x) for x in next(parts).split())
        values = np.array([float(next(parts)) for _ in range(capacity)])
        if status or count != len(expected):
            raise AssertionError(f"case {idx} n={n} kd={kd} mode={mode}: status={status}, count={count}, expected={len(expected)}")
        if np.any(values[count:] != 0.):
            raise AssertionError(f"case {idx} did not zero the unused output")
        if count:
            denom = max(1e-300,float(np.max(np.abs(wr))))
            diff = float(np.max(abs(values[:count]-expected))/denom)
            worst = max(worst,diff)
            if not np.isfinite(diff) or diff > tolerance:
                raise AssertionError(f"case {idx}: error={diff:.4g}, tol={tolerance}, n={n}, kd={kd}, mode={mode}, got={values[:count]}, expected={expected}")
    return worst


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["checked","optimized","both"],default="both")
    p.add_argument("--cases",type=int,default=120)
    p.add_argument("--max-order",type=int,default=64)
    p.add_argument("--seed",type=int,default=9021)
    p.add_argument("--tol",type=float,default=1.e-10)
    args = p.parse_args()
    cases = list(make_cases(args.cases,args.max_order,args.seed))
    modes = ["checked","optimized"] if args.mode=="both" else [args.mode]
    with tempfile.TemporaryDirectory(prefix="scifort_band_sturm_") as wd:
        work = Path(wd)
        for mode in modes:
            native = compile_program(work,mode,"native")
            run = subprocess.run([str(native)], cwd=work, capture_output=True,
                                 text=True, timeout=90)
            if run.returncode:
                raise RuntimeError(f"{mode} native tests failed:\n{run.stdout}{run.stderr}")
            exe = compile_program(work,mode,"driver")
            worst = check_cases(exe,cases,args.tol)
            print(f"{mode}: native passed; {len(cases)} SciPy eigenvalue comparisons passed; "
                  f"worst normalized eigenvalue error={worst:.3e}")


if __name__ == "__main__":
    main()
