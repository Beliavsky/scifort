#!/usr/bin/env python3
"""Selected tridiagonal eigenvectors: SciPy comparison and large-order test.

Independent reference checks focus on sorted eigenvalues, normalized residuals,
orthogonality, and projector agreement for resolved simple eigenvalues. For
multiple eigenvalues an orthonormal basis is not unique and direct comparison
of eigenvector entries is inappropriate.
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


def build(work: Path, source: str, mode: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / (Path(source).stem + "_" + mode)
    cmd = ["gfortran", "-std=f2018", *flags,
           *[str(ROOT / s) for s in SOURCES], str(ROOT / source), "-o", str(exe)]
    out = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if out.returncode:
        raise RuntimeError(f"{mode} build failed:\n{out.stderr}")
    return exe


def cases(count: int, seed: int, max_order: int):
    rng = np.random.default_rng(seed)
    possible_n = sorted(set([2, 3, 4, 5, 8, 16, 32, 64, 128, max_order]))
    for case in range(count):
        n = int(rng.choice(possible_n))
        d = rng.normal(size=n)
        e = rng.normal(size=n-1)
        family = case % 14
        if family == 0:  # a zero matrix, with highly non-unique eigenvectors
            d[:] = 0.
            e[:] = 0.
        elif family == 1:
            d[:] = 2.
            e[:] = -1.
        elif family == 2:
            d[:] = 2.
            e[:] = 1e-7
        elif family == 3:
            d[:] = 2.
            e[:] = 1e-11
        elif family == 4:
            d[:] = rng.choice([-2., 0., 1., 3.], size=n)
            e[:] = 0.
        elif family == 5:
            d[:] = 1.
            e[:] = 0.
            if n > 2:
                e[0] = 0.25
        elif family == 6:
            e[::3] = 0.  # independently repeated block spectra
        elif family == 7:
            d[:] = 1e12
            e *= 0.1
        elif family == 8:
            d *= 1e150
            e *= 1e150
        elif family == 9:
            d *= 1e-150
            e *= 1e-150
        elif family == 10:
            d[:] = 2.
            e[:] = 1e-5
            e[::2] = 0.
        elif family == 11:
            d[:] = np.linspace(-1, 1, n)
            e[:] = 1e-9
        elif family == 12:
            d[:] = 0.
            e[:] = 1.
        else:
            d[:] = 0.
            e[:] = 0.
            e[::2] = 1.

        first = int(rng.integers(1, n+1))
        k = int(rng.integers(1, min(7,n-first+2)))
        last = min(n, first+k-1)
        yield n, first, last, d, e


def check(exe: Path, inputs: list, tolerance: float):
    payload = [str(len(inputs))]
    for n, first, last, d, e in inputs:
        payload.append(f"{n} {first} {last}")
        payload.append(" ".join(f"{x:.17e}" for x in d))
        payload.append(" ".join(f"{x:.17e}" for x in e))
    proc = subprocess.run([str(exe)], input="\n".join(payload)+"\n",
                          text=True, capture_output=True, timeout=180,
                          env={**os.environ,"OPENBLAS_NUM_THREADS":"1", "OMP_NUM_THREADS":"1"})
    if proc.returncode:
        raise RuntimeError(f"driver failed: {proc.stderr[-1200:]}\n{proc.stdout[-1200:]}")
    lines = iter(proc.stdout.splitlines())
    worst_eig, worst_res, worst_orth, worst_angle = 0.,0.,0.,0.
    for index, (n, first, last, d, e) in enumerate(inputs):
        st, count = map(int, next(lines).split())
        m = last-first+1
        if st != 0 or count != m:
            raise AssertionError(f"case {index}: n={n}, st={st}, count={count}, family={index%14}")
        vals = np.array([float(next(lines)) for _ in range(m)])
        z = np.array([float(next(lines)) for _ in range(n*m)]).reshape((n,m),order="F")
        ref, zr = linalg.eigh_tridiagonal(d,e, select="i",
                                         select_range=(first-1,last-1))
        scale = max(np.max(abs(d)),np.max(abs(e)),1e-300)
        eig_err = np.max(abs(ref/scale-vals/scale))
        a_z = (d/scale)[:,None]*z
        a_z[1:] += (e/scale)[:,None]*z[:-1]
        a_z[:-1] += (e/scale)[:,None]*z[1:]
        residual = linalg.norm(a_z-z*(vals/scale),ord="fro") / np.sqrt(m)
        orth = linalg.norm(z.T@z-np.eye(m),ord="fro") / m
        if not (np.all(np.isfinite(vals)) and np.all(np.isfinite(z))):
            raise AssertionError(f"case {index}: nonfinite values")
        worst_eig=max(worst_eig,float(eig_err))
        worst_res=max(worst_res,float(residual))
        worst_orth=max(worst_orth,float(orth))
        if eig_err > tolerance or residual > tolerance or orth > tolerance:
            raise AssertionError(f"case {index} n={n} family={index%14}: "
                                 f"eig={eig_err:.3g}, res={residual:.3g}, orth={orth:.3g}")
        # Direct comparison only for well-resolved interior simple eigenvalues.
        full=linalg.eigvalsh_tridiagonal(d,e)
        ext=np.r_[np.inf, np.diff(full), np.inf]
        isolated = all(min(ext[first-1+i:first+1+i]) > 1e-4*scale for i in range(m))
        if isolated:
            angles=np.abs(np.sum(z*zr,axis=0))
            worst_angle=max(worst_angle,float(max(1-angles)))
            if np.min(angles) < 1-1e-7:
                raise AssertionError(f"case {index}: wrong isolated eigenvector")
    return worst_eig, worst_res, worst_orth, worst_angle


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode",choices=["checked","optimized","both"],default="both")
    parser.add_argument("--cases",type=int,default=280)
    parser.add_argument("--seed",type=int,default=1441)
    parser.add_argument("--max-order",type=int,default=128)
    parser.add_argument("--tol",type=float,default=8e-12)
    args=parser.parse_args()
    if args.cases < 1 or args.max_order < 2:
        parser.error("--cases >= 1 and --max-order >= 2 required")
    inputs=list(cases(args.cases,args.seed,args.max_order))
    modes=["checked","optimized"] if args.mode=="both" else [args.mode]
    with tempfile.TemporaryDirectory(prefix="scifort_tri_vec_") as td:
        work=Path(td)
        for mode in modes:
            native=build(work,"test/test_linalg_tridiagonal_vectors.f90",mode)
            subprocess.run([str(native)],cwd=work,check=True,stdout=subprocess.PIPE,
                           text=True,timeout=45)
            driver=build(work,"tools/linalg_reference/driver_tridiagonal_vectors.f90",mode)
            worst=check(driver, inputs, args.tol)
            print(f"{mode}: native passed; {len(inputs)} SciPy cases passed; "
                  f"eigen={worst[0]:.3e}, residual={worst[1]:.3e}, "
                  f"orth={worst[2]:.3e}, simple-vector error={worst[3]:.3e}")


if __name__ == "__main__":
    main()
