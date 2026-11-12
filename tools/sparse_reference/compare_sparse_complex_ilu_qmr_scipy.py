#!/usr/bin/env python3
"""Randomized SciPy/NumPy checks for complex BiCG/QMR and complex ILU(0).

Usage:
    python compare_sparse_complex_ilu_qmr_scipy.py DRIVER --per-solver 100

SciPy is only a reference dependency; the Fortran implementation has none.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import subprocess
import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sla


@dataclass
class Case:
    matrix: np.ndarray
    rows: np.ndarray
    cols: np.ndarray
    values: np.ndarray
    b: np.ndarray
    x0: np.ndarray
    solver: int
    fmt: int
    mode: int
    cap: int


def complex_random(rng, shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def case_random(rng, solver, n, fmt, mode):
    if solver == 5 and n > 4 and rng.random() < 0.5:
        # Tridiagonal: ILU(0) is exact when no fill is generated.
        mask = np.eye(n, dtype=bool) | np.eye(n, k=1, dtype=bool) | np.eye(n, k=-1, dtype=bool)
    else:
        mask = rng.random((n,n)) < 0.16
        mask |= np.eye(n, dtype=bool)
    mat = complex_random(rng, (n, n)) * mask
    diag = 3.5 + n*0.7 + np.arange(n)*0.03 + (rng.normal(size=n)*0.15)*1j
    mat[np.diag_indices(n)] = diag
    rows, cols = np.nonzero(mat)
    vals = mat[rows,cols].copy()
    # Random duplicated entries plus explicit zeros test canonical construction.
    if len(vals) > 0:
        idx = rng.integers(len(vals), size=min(len(vals),5))
        rr, cc = rows[idx], cols[idx]
        other = vals[idx]*0.4
        vals[idx] *= 0.6
        rows = np.r_[rows,rr]
        cols = np.r_[cols,cc]
        vals = np.r_[vals,other]
    zcount = int(rng.integers(0,3))
    rows = np.r_[rows, rng.integers(n,size=zcount)]
    cols = np.r_[cols, rng.integers(n,size=zcount)]
    vals = np.r_[vals, np.zeros(zcount,complex)]
    permutation = rng.permutation(len(vals))
    rows,cols,vals = rows[permutation], cols[permutation], vals[permutation]
    A = sparse.coo_array((vals,(rows,cols)),shape=(n,n)).toarray()
    truth = complex_random(rng,n)
    b = A @ truth if solver != 5 else complex_random(rng,n)
    x0 = complex_random(rng,n)*0.2 if rng.random() < 0.65 else np.zeros(n,complex)
    return Case(A,rows,cols,vals,b,x0,solver,fmt,mode,max(40,n*8))


def ilu0_dense(a,rows,cols):
    n = a.shape[0]
    mask = np.zeros((n,n),dtype=bool)
    mask[rows,cols] = True
    work = a.copy()
    lower = np.eye(n, dtype=complex)
    for i in range(n):
        for j in range(i):
            if not mask[i,j]:
                continue
            factor = work[i,j] / work[j,j]
            lower[i,j] = factor
            for k in range(j+1,n):
                if mask[j,k] and mask[i,k]:
                    work[i,k] -= factor * work[j,k]
        # Store only the upper half in U.
    upper = np.triu(work)
    return lower @ upper


def pack(case):
    n = case.matrix.shape[0]
    lines = [f'{n} {len(case.values)} {case.solver} {case.fmt} {case.mode} {case.cap}']
    for i,j,z in zip(case.rows,case.cols,case.values):
        lines.append(f'{i+1} {j+1} {z.real:.17e} {z.imag:.17e}')
    for vector in (case.b,case.x0):
        for z in vector:
            lines.append(f'{z.real:.17e} {z.imag:.17e}')
    return '\n'.join(lines)+'\n'


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('executable')
    ap.add_argument('--per-solver',type=int,default=80)
    ap.add_argument('--seed',type=int,default=1728)
    args=ap.parse_args()
    rng=np.random.default_rng(args.seed)
    cases=[]
    for solver in [1,2,3,4,5]:
        for j in range(args.per_solver):
            n=int(rng.integers(1,18))
            cases.append(case_random(rng,solver,n,1+j%2,(j//2)%3 if solver!=5 else 2))
    proc=subprocess.run([args.executable],input=''.join(pack(c) for c in cases),
                        capture_output=True,text=True,check=True,timeout=160)
    lines=iter(proc.stdout.splitlines())
    worst=0.0
    for ind,c in enumerate(cases):
        n=c.matrix.shape[0]
        try:
            info_s,iterations_s,resid_s=next(lines).split()
            x=np.array([complex(*map(float,next(lines).split())) for _ in range(n)])
        except (ValueError,StopIteration) as exc:
            raise AssertionError(f'case {ind} malformed output') from exc
        info,iterations,resid=int(info_s),int(iterations_s),float(resid_s)
        if c.solver==5:
            mat=ilu0_dense(c.matrix,c.rows,c.cols)
            expected=np.linalg.solve(mat.conj().T if c.fmt==2 else mat,c.b)
            reference_resid=0.0
        else:
            expected=np.linalg.solve(c.matrix,c.b)
            reference_resid=np.linalg.norm(c.b-c.matrix @ x)
            if c.mode != 2:
                a=sparse.csr_array(c.matrix)
                method={1:sla.bicg,2:sla.qmr,3:sla.bicgstab,4:sla.gmres}[c.solver]
                kwargs=dict(x0=c.x0,rtol=1e-10,atol=1e-12,maxiter=c.cap)
                if c.solver==4:
                    kwargs['restart']=min(20,n)
                scipy_solution, code=method(a,c.b,**kwargs)
                if code!=0 or not np.allclose(scipy_solution,expected,atol=5e-9,rtol=5e-9):
                    raise AssertionError(f'SciPy reference failed at {ind}: {code}')
        err=np.linalg.norm(x-expected)
        worst=max(worst,err)
        if info!=0 or not np.allclose(x,expected,atol=3e-9,rtol=5e-9):
            raise AssertionError(f'case {ind}, n={n}, solver={c.solver}, fmt={c.fmt}, mode={c.mode}: '
                                 f'info={info}, steps={iterations}, error={err:.5g}, residual={resid:.5g}')
        if c.solver!=5 and not np.isclose(resid,reference_resid,atol=2e-10,rtol=2e-7):
            raise AssertionError(f'case {ind}: actual residual={reference_resid:.5g}, reported={resid:.5g}')
    if next(lines,None) is not None:
        raise AssertionError('extra Fortran output')
    print(f'PASS {len(cases)} randomized SciPy/NumPy comparisons, seed={args.seed}, max solution error {worst:.3e}')


if __name__=='__main__':
    main()
