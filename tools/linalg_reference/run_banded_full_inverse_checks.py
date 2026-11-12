#!/usr/bin/env python3
"""Check opt-in full banded inverse iteration against SciPy and residuals."""
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
           "src/linalg_hessenberg.f90", "src/linalg_banded_tridiagonal.f90", "src/linalg_banded_eigen.f90"]


def compile_driver(work: Path, mode: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / f"banded_full_inverse_{mode}"
    cmd = ["gfortran", "-std=f2018", *flags, *[str(ROOT / f) for f in SOURCES],
           str(ROOT / "tools/linalg_reference/driver_banded_full_inverse.f90"), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"compilation failed: {p.stderr[-4000:]}")
    return exe


def cases(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = sorted({2, 3, 5, 8, 12, min(20, max_order), max_order})
    for it in range(count):
        n = int(rng.choice(sizes))
        kd = int(rng.choice([0, 1, 2, 3, min(5, n-1)]))
        kd = min(kd, n-1)
        typ = 1 + it % 4
        low = typ in (1, 3)
        complex_case = typ >= 3
        a = np.zeros((n,n),dtype=complex if complex_case else float)
        np.fill_diagonal(a, rng.normal(size=n))
        for j in range(n):
            for k in range(1,min(kd,n-j-1)+1):
                val = 0.3*rng.normal()+(0.3j*rng.normal() if complex_case else 0.)
                a[j+k,j]=val
                a[j,j+k]=np.conj(val)
        family=(it//4)%9
        if family == 0:
            a += 5*np.eye(n)
        elif family == 1:
            a = np.diag(rng.choice([-2., 0., 0., 1., 3.], size=n)).astype(a.dtype)
        elif family == 2:
            a *= 1e60
        elif family == 3:
            a *= 1e-60
        elif family == 4:
            a = np.diag(np.linspace(-2.,2.,n)).astype(a.dtype)
            for j in range(n-1):
                if kd:
                    a[j+1,j]=0.005
                    a[j,j+1]=0.005
        elif family == 5:
            a = np.diag(np.ones(n)).astype(a.dtype)
            for j in range(n-1):
                if kd:
                    a[j+1,j]=-0.15
                    a[j,j+1]=-0.15
        elif family == 6:
            if kd:
                for j in range(n):
                    for k in range(1,min(kd,n-j-1)+1):
                        if j//4 != (j+k)//4:
                            a[j+k,j]=a[j,j+k]=0
        elif family == 7:
            a = np.diag(np.ones(n)).astype(a.dtype)
            if kd:
                for j in range(n-1):
                    a[j+1,j]=1e-5
                    a[j,j+1]=1e-5
        else:
            a += np.eye(n)
        ab = np.zeros((kd+1,n),dtype=a.dtype)
        for j in range(n):
            for k in range(min(kd, n-j-1)+1):
                if low:
                    ab[k,j]=a[j+k,j]
                else:
                    ab[kd-k,j+k]=a[j,j+k]
        yield typ, n, kd, a, ab, family


def run(exe: Path, items, threshold: float):
    payload = [str(len(items))]
    for typ,n,kd,a,ab,family in items:
        payload.append(f"{typ} {n} {kd}")
        payload.extend(f"{v.real:.17e} {v.imag:.17e}" for v in ab.ravel(order="F"))
    p=subprocess.run([str(exe)],input="\n".join(payload)+"\n",capture_output=True,
                     text=True, timeout=180,env={**os.environ,"OPENBLAS_NUM_THREADS":"1",
                                                 "OMP_NUM_THREADS":"1"})
    if p.returncode:
        raise RuntimeError(f"Fortran execution error: {p.stderr[-3000:]}")
    lines=iter(p.stdout.splitlines())
    worst={"eig":0., "residual":0., "orth":0.}
    failures=[]
    for idx,(typ,n,kd,a,ab,family) in enumerate(items):
        status,found=map(int,next(lines).split())
        if status or found != n:
            failures.append((idx,typ,n,kd,family,status,found))
            continue
        w=np.array([float(next(lines)) for j in range(n)])
        z=np.array([complex(*map(float,next(lines).split())) for j in range(n*n)])
        z=z.reshape((n,n),order='F')
        wr=linalg.eig_banded(ab,lower=typ in (1,3),eigvals_only=True)
        scale=max(float(np.max(abs(a))),1e-300)
        norm=max(1.,float(linalg.norm(a/scale,'fro')))
        eig=float(np.max(abs((w-wr)/scale)))/norm
        res=float(linalg.norm((a/scale)@z-z*(w/scale)[None,:],'fro'))/norm
        orth=float(linalg.norm(z.conj().T@z-np.eye(n),'fro'))/max(1.,float(n))
        for k,v in (("eig",eig),("residual",res),("orth",orth)):
            worst[k]=max(worst[k],v)
            if not np.isfinite(v) or v > threshold:
                raise AssertionError(f"case {idx} n={n} kd={kd} type={typ} family={family}: "
                                     f"{k}={v:.3e} > {threshold:g}")
    if failures:
        raise AssertionError(f"failed {len(failures)}/{len(items)} cases: {failures[:15]}")
    return worst


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["checked","optimized","both"],default="both")
    ap.add_argument("--cases",type=int,default=80)
    ap.add_argument("--max-order",type=int,default=32)
    ap.add_argument("--seed",type=int,default=101)
    ap.add_argument("--tolerance",type=float,default=2e-8)
    args=ap.parse_args()
    items=list(cases(args.cases,args.max_order,args.seed))
    modes=["checked","optimized"] if args.mode=="both" else [args.mode]
    with tempfile.TemporaryDirectory(prefix="scifort_banded_full_") as temp:
        for mode in modes:
            work=Path(temp)
            exe=compile_driver(work,mode)
            worst=run(exe,items,args.tolerance)
            print(f"{mode}: {len(items)} comparisons passed; {worst=}")

if __name__ == "__main__":
    main()
