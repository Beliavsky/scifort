#!/usr/bin/env python3
"""Compare selected SciFort tridiagonal eigenvalues/eigenvectors against SciPy.

Use ``--mode both`` for separate runtime-checked and optimized GNU Fortran builds.
Value intervals are half-open (low, high]; index intervals are inclusive and
one-based on the SciFort side (zero-based in SciPy).
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


def compile_program(work: Path, name: str, mode: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace",
              "-ffpe-trap=invalid,zero,overflow"] if mode == "checked" else ["-O3"])
    exe = work / f"{name}_{mode}"
    source = ("test/test_linalg_tridiagonal_selection.f90" if name == "native" else
              "tools/linalg_reference/driver_tridiagonal_selection.f90")
    cmd = ["gfortran", "-std=f2018", *flags,
           *[str(ROOT / s) for s in SOURCES], str(ROOT / source), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"{mode} {name} compilation failed:\n{p.stderr}")
    return exe


def make_cases(total: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = sorted(set([1, 2, 3, 4, 5, 8, 12, 20, 32, max_order]))
    for case in range(total):
        n = int(rng.choice(sizes))
        d = rng.normal(size=n)
        e = rng.normal(size=max(0, n-1))
        family = (case//4) % 12
        if family == 0:
            d[:] = 0
            e[:] = 0
        elif family == 1:
            d[:] = 1
            e[:] = 0
        elif family == 2:
            d[:] = np.arange(n) - n/2
            e[:] = 0
        elif family == 3:
            d[:] = 2
            e[:] = -1
        elif family == 4:
            d[:] = 0
            e[:] = 1
        elif family == 5:
            d *= 1e150
            e *= 1e150
        elif family == 6:
            d *= 1e-150
            e *= 1e-150
        elif family == 7:
            d[:] = 2
            e[:] = 1e-7
        elif family == 8:
            d *= 1e-8
            e *= 1e-8
        elif family == 9:
            d *= 1e8
            e *= 1e8
        elif family == 10:
            d[:] = rng.choice([-2., -1., 0., 1., 2.], size=n)
            e[:] = 0
        else:
            d[:] = 2
            e[:] = 0
            if n > 1:
                e[0] = 1e-6
        mode = 1+case%4
        full = linalg.eigvalsh_tridiagonal(d,e,check_finite=True)
        if mode % 2 == 1:
            first = int(rng.integers(1,n+1))
            last = int(rng.integers(first,n+1))
            low,high = 0.,1.
            expected = full[first-1:last]
            selected_reference = linalg.eigvalsh_tridiagonal(
                d,e,select='i',select_range=(first-1,last-1))
        else:
            first,last = 1,n
            if case % 10 == 1:
                low = float(full[-1]+max(1.,abs(full[-1]))*1e-5)
                high = float(low+max(1.,abs(low))*1e-5)
            elif case % 10 == 3:
                low = float(full[0]-max(1.,abs(full[0]))*1e-5)
                high = float(full[-1]+max(1.,abs(full[-1]))*1e-5)
            else:
                # Mid-gaps avoid ambiguous inclusions at rounded boundaries.
                ka = int(rng.integers(0,n+1))
                kb = int(rng.integers(ka,n+1))
                scale = max(1e-300,float(max(abs(full))))
                low = float((full[ka-1]+full[ka])/2 if 0 < ka < n and
                            full[ka-1] < full[ka] else
                            full[0]-max(scale,1e-300)*1e-4 if ka == 0 else
                            full[-1]+max(scale,1e-300)*1e-4 if ka == n else
                            full[ka-1])
                high = float((full[kb-1]+full[kb])/2 if 0 < kb < n and
                             full[kb-1] < full[kb] else
                             full[0]-max(scale,1e-300)*1e-4 if kb == 0 else
                             full[-1]+max(scale,1e-300)*1e-4 if kb == n else
                             full[kb-1])
                if not low < high:
                    low = float(full[0]-max(scale,1e-300)*1e-4)
                    high = float(full[-1]+max(scale,1e-300)*1e-4)
            expected = full[(full>low) & (full<=high)]
            selected_reference = linalg.eigvalsh_tridiagonal(
                d,e,select='v',select_range=(low,high))
            if len(expected) != len(selected_reference):
                # Numerical interval endpoints can be ambiguous if eigenvalues
                # lie within one rounding unit of the specified boundary.
                expected = selected_reference
        yield mode, n, first, last, low, high, d, e, expected, selected_reference


def run_comparisons(exe: Path, cases: list, tol: float):
    data = [str(len(cases))]
    for mode,n,first,last,low,high,d,e,_,_ in cases:
        data.append(f"{mode} {n} {first} {last} {low:.17e} {high:.17e}")
        data.append(" ".join(f"{v:.17e}" for v in d))
        if n>1:
            data.append(" ".join(f"{v:.17e}" for v in e))
    env = {**os.environ, "OPENBLAS_NUM_THREADS":"1", "OMP_NUM_THREADS":"1"}
    p = subprocess.run([str(exe)], input="\n".join(data)+"\n", text=True,
                       capture_output=True, timeout=240, env=env)
    if p.returncode:
        raise RuntimeError(f"Fortran driver failed:\n{p.stdout[-1000:]}\n{p.stderr[-2000:]}")
    lines = iter(p.stdout.splitlines())
    max_eigen_error, max_residual, max_orth, max_ref_gap = 0.,0.,0.,0.
    for index,(mode,n,first,last,low,high,d,e,ref,selected_ref) in enumerate(cases):
        st,m = map(int,next(lines).split())
        if st != 0 or m != len(selected_ref):
            raise AssertionError(f"case {index}: status={st}, count={m}, "
                                 f"SciPy selected count={len(selected_ref)}, "
                                 f"mode={mode}, n={n}, interval=({low:.5g},{high:.5g}]")
        found = np.array([float(next(lines)) for _ in range(m)])
        scale = max(float(np.max(np.abs(d))), float(np.max(np.abs(e))) if n>1 else 0.,
                    1.e-300)
        eig_error = max(abs(found-ref),default=0)/scale
        max_eigen_error = max(max_eigen_error,eig_error)
        ref_gap = max(abs(selected_ref-ref),default=0)/scale
        max_ref_gap = max(max_ref_gap,ref_gap)
        if eig_error > tol or not np.all(np.isfinite(found)):
            raise AssertionError(f"case {index}: eigen error={eig_error:.3e}, "
                                 f"mode={mode}, n={n}, d={d}, e={e}, "
                                 f"result={found}, expected={ref}")
        if mode >= 3 and m > 0:
            z = np.array([float(next(lines)) for _ in range(n*m)]).reshape((n,m),order='F')
            az = d[:,None]*z
            az[1:] += e[:,None]*z[:-1]
            az[:-1] += e[:,None]*z[1:]
            residual = float(linalg.norm(az-z*found,ord='fro')) / (scale*max(1.,m**0.5))
            orth = float(linalg.norm(z.T@z-np.eye(m),ord='fro'))/m
            max_residual = max(max_residual,residual)
            max_orth = max(max_orth,orth)
            if residual > tol or orth > tol:
                raise AssertionError(f"case {index}, n={n}: residual={residual}, orth={orth}")
    return max_eigen_error,max_residual,max_orth,max_ref_gap


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode",choices=["both","checked","optimized"],default="both")
    parser.add_argument("--cases",type=int,default=350)
    parser.add_argument("--max-order",type=int,default=64)
    parser.add_argument("--seed",type=int,default=20261108)
    parser.add_argument("--tol",type=float,default=5e-12)
    args = parser.parse_args()
    if args.cases < 1 or args.max_order < 1 or args.tol <= 0:
        parser.error("--cases, --max-order and --tol must be positive")
    cases = list(make_cases(args.cases,args.max_order,args.seed))
    modes = ["checked","optimized"] if args.mode == "both" else [args.mode]
    with tempfile.TemporaryDirectory(prefix="scifort_tri_select_") as folder:
        work = Path(folder)
        for mode in modes:
            native = compile_program(work,"native",mode)
            p = subprocess.run([str(native)],cwd=work,capture_output=True,text=True,timeout=90)
            if p.returncode:
                raise RuntimeError(f"{mode} native tests failed:\n{p.stdout}\n{p.stderr}")
            driver = compile_program(work,"driver",mode)
            worst = run_comparisons(driver,cases,args.tol)
            print(f"{mode}: native passed, {len(cases)} SciPy comparisons passed; "
                  f"eigenvalue={worst[0]:.3e}, residual={worst[1]:.3e}, "
                  f"orthogonality={worst[2]:.3e}, "
                  f"SciPy selected-vs-full gap={worst[3]:.3e}")


if __name__ == "__main__":
    main()
