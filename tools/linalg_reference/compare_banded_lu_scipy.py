#!/usr/bin/env python3
"""Generate independent real/complex, pivoted N/T/C band-LU comparisons."""
from __future__ import annotations

from pathlib import Path
import subprocess
import numpy as np
from scipy import linalg
from compare_banded_scipy import pack_general


def build_cases(count: int, seed: int, max_order: int, stress: bool) -> list:
    rng = np.random.default_rng(seed)
    cases = []
    q = 0
    while len(cases) < count:
        n = int(rng.integers(1,max_order+1))
        if stress and q % 47 == 0:
            n = [80,120,200,350][(q//47)%4]
        maxbw = min(n-1, 8) if q % 13 else n-1
        if n > 70:
            maxbw = min(n-1, 5)
        lower = int(rng.integers(0,maxbw+1))
        upper = int(rng.integers(0,maxbw+1))
        dtype = np.float64 if q % 2 == 0 else np.complex128
        a = np.zeros((n,n),dtype=dtype)
        for j in range(n):
            for i in range(max(0,j-upper),min(n,j+lower+1)):
                a[i,j] = rng.normal() + (1j*rng.normal() if dtype==np.complex128 else 0)
        # Very large random band matrices may be exponentially ill-conditioned
        # even when nonsingular; use a diagonally dominant large-order family.
        if n > 70:
            for i in range(n):
                a[i,i] += 1.5 + np.sum(abs(a[i,:]))
            if lower and q % 7 == 0:
                a[0,0] = 1.e-3
                a[1,0] = 3.0
        elif q % 5 == 0:
            for i in range(n):
                a[i,i] += 1.2 + np.sum(abs(a[i,:]))
        elif q % 5 == 1 and lower:
            a[0,0] = 1e-4
            a[1,0] = 4.0
        elif q % 5 == 2 and lower:
            a[0,0] = 1e-7
            a[min(n-1,lower),0] = 6.0
        elif q % 5 == 3:
            for i in range(n):
                a[i,i] += 0.1
        if n <= 70 and np.linalg.cond(a) > 1.e7:
            q += 1
            continue
        # Two independent batches of multiple right-hand sides reuse the same LU.
        nrhs = 1 + q % 3
        b = [rng.normal(size=(n,nrhs)) for _ in range(2)]
        if dtype == np.complex128:
            b = [v + 1j*rng.normal(size=v.shape) for v in b]
        try:
            expected = []
            for op in ('N','T','C'):
                aop = a if op=='N' else a.T if op=='T' else a.conj().T
                expected.append([linalg.solve(aop, rhs, assume_a='gen', check_finite=True) for rhs in b])
        except np.linalg.LinAlgError:
            q += 1
            continue
        cases.append((a,pack_general(a,lower,upper),lower,upper,b,expected))
        q += 1
    return cases


def run(exe: Path, cases: list, tol: float = 3.e-10):
    lines=[str(len(cases))]
    for a, ab, lower, upper, bs, expected in cases:
        n=a.shape[0]
        nrhs=bs[0].shape[1]
        lines.append(f'{1 if a.dtype == np.float64 else 2} {lower} {upper} {n} {nrhs}')
        for z in ab.reshape(-1,order='F'):
            lines.append(f'{z.real:.17e} {z.imag:.17e}')
        for b in bs:
            for z in b.reshape(-1,order='F'):
                lines.append(f'{z.real:.17e} {z.imag:.17e}')
    completed=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,capture_output=True)
    if completed.returncode:
        raise RuntimeError(f'Fortran executable exited {completed.returncode}: {completed.stderr[-3000:]}')
    output=iter(completed.stdout.splitlines())
    max_err=0.0
    max_residual=0.0
    for index,(a,ab,lower,upper,bs,expected) in enumerate(cases):
        n=a.shape[0]
        nrhs=bs[0].shape[1]
        status=int(next(output))
        if status:
            raise AssertionError(f'case={index} LU failed n={n}, bands={lower},{upper}, status={status}')
        for op in range(3):
            aop = a if op==0 else a.T if op==1 else a.conj().T
            for rep in range(2):
                status=int(next(output))
                if status:
                    raise AssertionError(f'case={index} LU solve failed op={op}, rep={rep}, status={status}')
                vals=np.array([list(map(float,next(output).split())) for _ in range(n*nrhs)])
                x=(vals[:,0]+1j*vals[:,1]).reshape((n,nrhs),order='F')
                wanted=expected[op][rep]
                err=float(np.linalg.norm(x-wanted)/max(float(np.linalg.norm(wanted)),1.e-300))
                residual=float(np.linalg.norm(aop@x-bs[rep])/
                    max(float(np.linalg.norm(aop)*np.linalg.norm(x)+np.linalg.norm(bs[rep])),1.e-300))
                max_err=max(err,max_err)
                max_residual=max(residual,max_residual)
                if not np.isfinite(err) or err>tol or residual>tol:
                    raise AssertionError(f'case={index} n={n} lower={lower} upper={upper} '
                                         f'op={op} rhs={rep} relerr={err:.2e} residual={residual:.2e}')
    print(f'{exe.name}: {len(cases)*6} SciPy comparisons passed (from {len(cases)} LU factors); '
          f'max relative error={max_err:.3e}, max backward residual={max_residual:.3e}')
