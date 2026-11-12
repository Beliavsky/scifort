#!/usr/bin/env python3
"""Independent randomized SciPy comparisons for SciFort packed-band solvers."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]


def pack_general(a: np.ndarray, lower: int, upper: int) -> np.ndarray:
    n = a.shape[0]
    packed = np.zeros((lower+upper+1,n), dtype=a.dtype)
    for j in range(n):
        for i in range(max(0,j-upper),min(n,j+lower+1)):
            packed[upper+i-j,j] = a[i,j]
    return packed


def pack_positive(a: np.ndarray, bw: int, lower: bool) -> np.ndarray:
    n = a.shape[0]
    packed = np.zeros((bw+1,n), dtype=a.dtype)
    for j in range(n):
        if lower:
            for i in range(j,min(n,j+bw+1)):
                packed[i-j,j] = a[i,j]
        else:
            for i in range(max(0,j-bw),j+1):
                packed[bw+i-j,j] = a[i,j]
    return packed


def build_cases(count: int, seed: int, max_order: int, stress: bool = False):
    rng = np.random.default_rng(seed)
    cases = []
    for q in range(count):
        n = int(rng.integers(1,max_order+1))
        # Some broad bands and some narrow bands, including zero bandwidth.
        bwmax = n-1 if q % 13 == 0 else min(n-1, 6)
        lower = int(rng.integers(0,bwmax+1))
        upper = int(rng.integers(0,bwmax+1))
        dtype = 2 if q % 2 else 1
        cmplx = dtype == 2
        nrhs = [1, 3, 2][q % 3]
        op = q % 3 + 1
        rhs = rng.normal(size=(n,nrhs))
        if cmplx:
            rhs = rhs + 1j*rng.normal(size=(n,nrhs))
        if op == 1:
            a = rng.normal(size=(n,n))
            if cmplx:
                a = a + 1j*rng.normal(size=(n,n))
            for j in range(n):
                for i in range(n):
                    if i-j > lower or j-i > upper:
                        a[i,j] = 0
            if stress and q % 5 == 3 and lower:
                # Favor a row swap on the first and later pivot columns.
                for i in range(n):
                    a[i,i] = 0.001
                a += np.eye(n,dtype=a.dtype)*0.001
            else:
                for i in range(n):
                    a[i,i] += sum(abs(a[i,:])) + 1.0
            if np.linalg.cond(a) > 1.e12:
                continue
            packed = pack_general(a,lower,upper)
            try:
                expected = linalg.solve_banded((lower,upper),packed,rhs)
            except np.linalg.LinAlgError:
                continue
            flag = 0
        else:
            bw = lower if q % 2 else upper
            lowerflag = bool(q % 2)
            l = np.zeros((n,n),dtype=np.complex128 if cmplx else np.float64)
            for j in range(n):
                l[j,j] = 1.5+rng.uniform(0,1)
                for i in range(j+1,min(n,j+bw+1)):
                    l[i,j] = 0.2*rng.normal() + (0.2j*rng.normal() if cmplx else 0)
            a = l @ l.conj().T
            packed = pack_positive(a,bw,lowerflag)
            lower = bw if lowerflag else 0
            upper = 0 if lowerflag else bw
            expected = (linalg.solveh_banded(packed,rhs,lower=lowerflag) if op==2 else
                        linalg.cho_solve_banded((linalg.cholesky_banded(packed,lower=lowerflag),lowerflag),rhs))
            flag = int(lowerflag)
        cases.append((op,dtype,lower,upper,n,nrhs,flag,packed,rhs,expected,a))
    return cases


def run(exe: Path, cases: list, tolerance: float):
    lines = [str(len(cases))]
    for op,dtype,lower,upper,n,nrhs,flag,packed,rhs,expected,a in cases:
        lines.append(f'{op} {dtype} {lower} {upper} {n} {nrhs} {flag}')
        for z in packed.reshape(-1,order='F'):
            lines.append(f'{z.real:.17e} {z.imag:.17e}')
        for z in rhs.reshape(-1,order='F'):
            lines.append(f'{z.real:.17e} {z.imag:.17e}')
    completed = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,capture_output=True)
    if completed.returncode:
        print(completed.stderr[-2000:])
        raise RuntimeError('Fortran driver failed')
    out = iter(completed.stdout.splitlines())
    worst = 0.0
    max_residual = 0.0
    for q,(op,dtype,lower,upper,n,nrhs,flag,packed,rhs,expected,a) in enumerate(cases):
        status = int(next(out))
        if status:
            raise AssertionError(f'Case {q} op={op} dtype={dtype} n={n} l={lower} u={upper} returned {status}')
        v = np.array([list(map(float,next(out).split())) for _ in range(n*nrhs)])
        x = (v[:,0]+1j*v[:,1]).reshape((n,nrhs),order='F')
        scale = max(float(np.linalg.norm(expected)),1.e-300)
        err = float(np.linalg.norm(x-expected)/scale)
        residual = float(np.linalg.norm(a@x-rhs) / max(float(np.linalg.norm(a)*np.linalg.norm(x) + np.linalg.norm(rhs)),1.e-300))
        worst=max(worst,err)
        max_residual=max(max_residual,residual)
        if not np.isfinite(err) or err > tolerance or residual > tolerance:
            raise AssertionError(f'Case {q} op={op} kind={dtype} n={n} l={lower} u={upper} err={err:.2g} resid={residual:.2g}')
    print(f'{exe.name}: {len(cases)} comparisons passed; maximum relative error={worst:.3e}; maximum residual={max_residual:.3e}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',required=True,type=Path)
    p.add_argument('--cases',type=int,default=400)
    p.add_argument('--seed',type=int,default=407)
    p.add_argument('--max-order',type=int,default=50)
    p.add_argument('--stress',action='store_true')
    p.add_argument('--tolerance',type=float,default=3e-10)
    a=p.parse_args()
    cases=build_cases(a.cases,a.seed,a.max_order,a.stress)
    run(a.executable.resolve(),cases,a.tolerance)


if __name__=='__main__':
    main()
