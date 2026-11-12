#!/usr/bin/env python3
"""Compare selected banded eigenpairs with SciPy, including large sparse bands."""
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


def compile_driver(work: Path, mode: str, source: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / f"{Path(source).stem}_{mode}"
    cmd = ["gfortran", "-std=f2018", *flags,
           *[str(ROOT / name) for name in SOURCES], str(ROOT / source), "-o", str(exe)]
    result = subprocess.run(cmd, cwd=work, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"{mode} compilation failed:\n{result.stderr}")
    return exe


def cases(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = sorted({2, 3, 4, 7, 12, 20, 36, min(64,max_order), max_order})
    for it in range(count):
        n = int(rng.choice(sizes))
        kd = min(n-1, int(rng.choice([0, 1, 2, 3, 4, 6])))
        typ = 1 + (it % 4)
        low = typ in (1, 3)
        complex_case = typ >= 3
        family = (it//4)%10
        a = np.zeros((n,n),dtype=complex if complex_case else float)
        np.fill_diagonal(a, rng.normal(size=n))
        for j in range(n):
            for k in range(1,min(kd,n-j-1)+1):
                value = (0.4*rng.normal() + (0.3j*rng.normal() if complex_case else 0))
                a[j+k,j]=value
                a[j,j+k]=np.conj(value)
        if family == 0:
            a += 3.*np.eye(n)
        elif family == 1:
            a *= 0.
            np.fill_diagonal(a, 2.)
            if kd:
                for j in range(n-1):
                    value=(-0.4j if complex_case else -0.4)
                    a[j+1,j]=value
                    a[j,j+1]=np.conj(value)
        elif family == 2:
            a = np.diag(rng.choice([-1.,0.,0.,2.,3.],size=n)).astype(a.dtype)
        elif family == 3:
            a *= 1e70
        elif family == 4:
            a *= 1e-70
        elif family == 5:
            a = np.diag(np.linspace(-2.,2.,n)).astype(a.dtype)
            if kd:
                for j in range(n-1):
                    a[j+1,j]=0.02
                    a[j,j+1]=0.02
        elif family == 6:
            a = np.diag(np.ones(n)).astype(a.dtype)
            if kd:
                for j in range(n-1):
                    a[j+1,j]=1e-5
                    a[j,j+1]=1e-5
        elif family == 7:
            # Repeated blocks, no coupling across some partitions.
            if kd:
                for j in range(n):
                    for k in range(1,min(kd,n-j-1)+1):
                        if j//4 != (j+k)//4:
                            a[j+k,j]=0
                            a[j,j+k]=0
        elif family == 8:
            a += 10.*np.eye(n)
        else:
            a += np.eye(n)
        # Choose one through four consecutive eigenvectors, often near the edges.
        m = min(n, int(rng.integers(1,5)))
        if it % 5 == 0:
            first=1
        elif it%5 == 1:
            first=n-m+1
        else:
            first=int(rng.integers(1,n-m+2))
        last=first+m-1
        ab=np.zeros((kd+1,n),dtype=a.dtype)
        for j in range(n):
            if low:
                for k in range(min(kd,n-j-1)+1):
                    ab[k,j]=a[j+k,j]
            else:
                for k in range(min(kd,j)+1):
                    ab[kd-k,j]=a[j-k,j]
        yield typ,n,kd,first,last,ab,a


def run(exe: Path, items, threshold: float):
    payload=[str(len(items))]
    for typ,n,kd,first,last,ab,a in items:
        payload.append(f"{typ} {n} {kd} {first} {last}")
        payload.extend(f"{v.real:.17e} {v.imag:.17e}" for v in ab.ravel(order="F"))
    proc=subprocess.run([str(exe)], input="\n".join(payload)+"\n", text=True,
                        capture_output=True, timeout=150,
                        env={**os.environ,"OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1"})
    if proc.returncode:
        raise RuntimeError(f"Fortran driver failed:\n{proc.stderr[-2000:]}")
    lines=iter(proc.stdout.splitlines())
    worst={"eig":0.,"residual":0.,"orth":0.}
    for idx,(typ,n,kd,first,last,ab,a) in enumerate(items):
        status,found=map(int,next(lines).split())
        m=last-first+1
        if status or found!=m:
            raise AssertionError(f"case {idx} n={n} kd={kd} family={(idx//4)%10} "
                                 f"type={typ} selection={first}:{last} status={status},found={found}")
        w=np.array([float(next(lines)) for _ in range(m)])
        z=np.array([complex(*map(float,next(lines).split())) for _ in range(n*m)])
        z=z.reshape((n,m),order='F')
        wr = linalg.eig_banded(ab,lower=typ in (1,3),eigvals_only=True)
        scale=max(float(np.max(abs(a))), 1e-300)
        norm=max(1.,float(linalg.norm(a/scale,'fro')))
        eigerr=float(np.max(abs((w-wr[first-1:last])/scale))/norm)
        res=float(linalg.norm((a/scale)@z-z*(w/scale)[None,:],'fro'))/norm
        ortho=float(linalg.norm(z.conj().T@z-np.eye(m),'fro'))/max(1,m)
        for key,value in (("eig",eigerr),("residual",res),("orth",ortho)):
            worst[key]=max(value,worst[key])
            if not np.isfinite(value) or value>threshold:
                raise AssertionError(f"case {idx}: {key}={value:.3g}, n={n} kd={kd} "
                                     f"type={typ} family={(idx//4)%10} limit={threshold:g}")
    return worst


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode",choices=["both","checked","optimized"],default="both")
    parser.add_argument("--cases",type=int,default=160)
    parser.add_argument("--max-order",type=int,default=80)
    parser.add_argument("--seed",type=int,default=15151)
    parser.add_argument("--tol",type=float,default=2e-8)
    args=parser.parse_args()
    tests=list(cases(args.cases,args.max_order,args.seed))
    with tempfile.TemporaryDirectory(prefix="scifort_band_vectors_") as folder:
        work=Path(folder)
        modes=["checked","optimized"] if args.mode=='both' else [args.mode]
        for mode in modes:
            native=compile_driver(work,mode,"test/test_linalg_banded_vectors.f90")
            subprocess.run([str(native)],check=True,capture_output=True,text=True,timeout=90)
            driver=compile_driver(work,mode,"tools/linalg_reference/driver_banded_vectors.f90")
            worst=run(driver,tests,args.tol)
            print(f"{mode}: native passed; {args.cases} SciPy comparisons passed; "
                  f"worst eigenvalue={worst['eig']:.3e}, residual={worst['residual']:.3e}, "
                  f"orthogonality={worst['orth']:.3e}")


if __name__=='__main__':
    main()
