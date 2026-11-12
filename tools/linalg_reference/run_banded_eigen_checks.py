#!/usr/bin/env python3
"""Validate SciFort symmetric tridiagonal/Hermitian band eigenpairs versus SciPy.

The comparison checks eigenvalues, eigen-equation residuals, and orthogonality;
degenerate eigenspaces can have different but equally valid bases.
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


def compile_program(work: Path, name: str, mode: str):
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    executable = work / (name + "_" + mode)
    cmd = ["gfortran", "-std=f2018", *flags, *[str(ROOT / s) for s in SOURCES],
           str(ROOT / ("test/test_linalg_banded_eigen.f90" if name == "native"
                       else "tools/linalg_reference/driver_banded_eigen.f90")),
           "-o", str(executable)]
    result = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError("Fortran compilation failed: " + result.stderr)
    return executable


def examples(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = sorted(set([1, 2, 3, 4, 5, 7, 11, 16, 24, max_order]))
    for it in range(count):
        n = int(rng.choice(sizes))
        typ = 1 + (it % 5)
        family = (it // 5) % 11
        bw = min(n - 1, int(rng.choice([0, 1, 2, 3, 6, n - 1])))
        if typ == 1:
            d = rng.standard_normal(n)
            e = rng.standard_normal(n-1)
            if family == 0:
                d[:] = 0
                e[:] = 0
            elif family == 1:
                d[:] = 1
                e[:] = 0
            elif family == 2:
                d[:] = 1.0
                e[:] = 1.e-7
            elif family == 3:
                d[:] = 0
                e[:] = 1.0
            elif family == 4:
                d[:] = np.arange(n) - n/2
            elif family == 5:
                d *= 1.e150
                e *= 1.e150
            elif family == 6:
                d *= 1.e-150
                e *= 1.e-150
            elif family == 7:
                d[:] = 2
                e[:] = 1
            elif family == 8:
                d = np.linspace(-3, 3, n)
                e *= 1.e-12
            elif family == 9:
                d = np.ones(n)
                e[:] = 0
                if n > 1:
                    e[0] = 1.e-6
            else:
                d *= 10.0**(-rng.uniform(0, 8, size=n))
                e *= 10.0**(-rng.uniform(0, 8, size=n-1))
            a = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
            yield typ, n, 0, np.concatenate((d, e)).astype(complex), a, (d, e)
            continue
        complex_case = typ >= 4
        lower = typ in (3, 5)
        dtype = complex if complex_case else float
        a = np.zeros((n, n), dtype=dtype)
        for j in range(n):
            a[j, j] = rng.standard_normal()
            for i in range(j+1, min(n,j+bw+1)):
                value = (rng.standard_normal() + (1j*rng.standard_normal()
                         if complex_case else 0))
                a[i, j] = value
                a[j, i] = np.conj(value)
        if family == 0:
            a *= 0
        elif family == 1:
            a = np.eye(n, dtype=dtype)
        elif family == 2:
            a *= 0
            np.fill_diagonal(a, rng.choice([-1., 1.], size=n))
        elif family == 3:
            a += np.eye(n) * (4*bw+4)
        elif family == 4:
            a *= 1.e150
        elif family == 5:
            a *= 1.e-150
        elif family == 6:
            a *= 1.e-9
            a += np.eye(n) * 10.0
        elif family == 7:
            a *= 1.e-12
            a += np.eye(n)
        elif family == 8:
            a = a.real.astype(dtype)  # complex matrix with real entries
        elif family == 9:
            a *= 10.0**rng.uniform(-2, 2)
        else:
            a = (a*0) + np.eye(n, dtype=dtype)
            if bw > 0:
                for j in range(n-1):
                    a[j+1,j] = 1.e-6j if complex_case else 1.e-6
                    a[j,j+1] = np.conj(a[j+1,j])
        ab = np.zeros((bw+1,n), dtype=dtype)
        for j in range(n):
            for k in range(min(bw,n-j-1)+1 if lower else min(bw,j)+1):
                if lower:
                    ab[k,j] = a[j+k,j]
                else:
                    ab[bw-k,j] = a[j-k,j]
        yield typ, n, bw, ab.ravel(order='F').astype(complex), a, ab


def run_comparisons(executable: Path, cases, tol: float):
    data = [str(len(cases))]
    for typ,n,bw,flat,a,payload in cases:
        data.append(f"{typ} {n} {bw}")
        data.extend(f"{z.real:.17e} {z.imag:.17e}" for z in flat)
    p = subprocess.run([str(executable)], input="\n".join(data)+"\n", capture_output=True,
                       text=True, timeout=240, env={**os.environ, "OPENBLAS_NUM_THREADS": "1",
                                                    "OMP_NUM_THREADS": "1"})
    if p.returncode:
        raise RuntimeError("Fortran batch failed:\n"+p.stderr[-2500:])
    lines = iter(p.stdout.splitlines())
    worst = dict(eigenvalue=0., residual=0., orthogonality=0.)
    for idx,(typ,n,bw,flat,a,payload) in enumerate(cases):
        status = int(next(lines))
        if status:
            raise AssertionError(f"case {idx}, type={typ}, n={n}: status {status}")
        w = np.array([float(next(lines)) for _ in range(n)])
        v = np.array([complex(*map(float,next(lines).split())) for _ in range(n*n)]
                     ).reshape((n,n),order="F")
        if typ == 1:
            wr, vr = linalg.eigh_tridiagonal(*payload)
        else:
            wr, vr = linalg.eig_banded(payload, lower=typ in (3,5))
        scale = max(1.e-300,float(np.max(abs(a))))
        aa,ww,wref = a/scale,w/scale,wr/scale
        denom = max(1.,float(linalg.norm(aa,"fro")))
        ev_error = np.max(abs(ww-wref))/denom
        res = linalg.norm(aa@v-v*ww[np.newaxis,:],"fro")/denom
        ortho = linalg.norm(v.conj().T@v-np.eye(n),"fro")/n
        for key,value in (("eigenvalue",ev_error),("residual",res),
                          ("orthogonality",ortho)):
            worst[key] = max(worst[key],float(value))
            if not np.isfinite(value) or value > tol:
                raise AssertionError(f"case {idx}, type={typ}, n={n}, {key}={value:.3g}, "
                                     f"tolerance={tol:g}")
    return worst


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["both","checked","optimized"],default="both")
    parser.add_argument("--cases",type=int,default=300)
    parser.add_argument("--seed",type=int,default=20261107)
    parser.add_argument("--max-order",type=int,default=48)
    parser.add_argument("--tol",type=float,default=5.e-12)
    args=parser.parse_args()
    if args.cases<1 or args.max_order<1 or args.tol<=0:
        parser.error("positive --cases, --max-order and --tol required")
    cases=list(examples(args.cases,args.max_order,args.seed))
    modes=["checked","optimized"] if args.mode=="both" else [args.mode]
    with tempfile.TemporaryDirectory(prefix="scifort_banded_eig_") as path:
        work=Path(path)
        for mode in modes:
            native=compile_program(work,"native",mode)
            run=subprocess.run([str(native)],cwd=work,capture_output=True,text=True,timeout=60)
            if run.returncode:
                raise RuntimeError(f"{mode} native tests failed:\n"+run.stdout+run.stderr)
            driver=compile_program(work,"driver",mode)
            worst=run_comparisons(driver,cases,args.tol)
            print(f"{mode}: native passed; {len(cases)} SciPy comparisons passed; " +
                  ", ".join(f"{k}={v:.3e}" for k,v in worst.items()))


if __name__=="__main__":
    main()
