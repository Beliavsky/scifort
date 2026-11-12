#!/usr/bin/env python3
"""Offline randomized validation for Fortran funm_multiply_krylov.

Usage: python compare_funm_krylov_scipy.py /path/to/funm_krylov_driver [--seed N] [--cases N]

Compares with dense SciPy matrix functions and SciPy's own Krylov action.
The Fortran driver also cross-checks the matrix-free and CSR overloads.
"""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class Case:
    matrix: np.ndarray
    rhs: np.ndarray
    t: float
    m: int
    max_restarts: int
    hermitian: bool
    function: int


def function_factory(mode: int):
    if mode == 0:
        return la.expm
    if mode == 1:
        return lambda x: x @ x @ x - x
    if mode == 2:
        return lambda x: np.eye(x.shape[0]) + x + 0.1 * (x @ x)
    if mode == 3:
        return la.sinm
    raise ValueError(mode)


def generate(seed: int, num: int) -> list[Case]:
    rng = np.random.default_rng(seed)
    tests = []
    for k in range(num):
        n = int(rng.integers(2, 23))
        hermitian = (k % 3 == 0)
        mat = np.zeros((n, n))
        if k % 7 == 0:
            mat[np.diag_indices(n)] = rng.normal(-0.5, 0.3, n)
        elif k % 7 == 1:
            mat[np.diag_indices(n)] = -0.4 + 0.1 * rng.normal(size=n)
            for j in range(n-1):
                mat[j,j+1] = 0.4 * rng.normal()
                mat[j+1,j] = mat[j,j+1] if hermitian else 0.4*rng.normal()
        else:
            mask = rng.uniform(size=(n,n)) < rng.uniform(0.08,0.22)
            mat[mask] = rng.normal(size=int(np.count_nonzero(mask)))
            mat[np.diag_indices(n)] += rng.normal(-0.5,0.25,n)
            if hermitian:
                mat = 0.5*(mat+mat.T)
        # Make a nontrivial but numerically controlled action.
        mat *= 1.1 / max(1.0,np.linalg.norm(mat,ord=1))
        b = rng.normal(size=n)
        if k % 19 == 0:
            b[:] = 0.0
        function = k % 4
        # Polynomial cases need enough Krylov vectors to accurately propagate.
        m = int(rng.integers(2, min(8,n)+1))
        if k % 29 == 0:
            m = n + 11     # oversized Krylov basis must behave as a full n-dimensional basis
        elif k % 5 == 0:
            m = n
        t = float(rng.uniform(-0.9,0.9))
        max_restarts = 20
        tests.append(Case(mat,b,t,m,max_restarts,hermitian,function))
    return tests


def check(driver: str, seed: int = 3381, cases: int = 250):
    inputs = []
    test_cases = generate(seed,cases)
    for case in test_cases:
        n=case.matrix.shape[0]
        inputs.append(f"{n} {case.m} {case.max_restarts} {int(case.hermitian)} {case.function} {case.t:.17g}")
        inputs.append(' '.join(f'{v:.17g}' for v in case.matrix.flatten(order='F')))
        inputs.append(' '.join(f'{v:.17g}' for v in case.rhs))
    result = subprocess.run([driver],input='\n'.join(inputs)+'\n',text=True,
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
    if result.returncode:
        raise RuntimeError(f'Fortran driver exited {result.returncode}: {result.stderr[-4000:]}')
    lines = result.stdout.splitlines()
    if len(lines) != 2*cases:
        raise RuntimeError(f'Expected {2*cases} lines, got {len(lines)}: {result.stdout[-1000:]}')
    worst_abs=0.0
    worst_rel=0.0
    worst_scipy=0.0
    max_mf=0.0
    unconverged=[]
    scipy_cases=0
    scipy_breakdown_bugs=0
    converged_cases=0
    num_good=0
    for k,case in enumerate(test_cases):
        n=case.matrix.shape[0]
        status,cycles,products,discrepancy = lines[2*k].split()
        status,cycles,products=int(status),int(cycles),int(products)
        discrepancy=float(discrepancy)
        got=np.fromstring(lines[2*k+1],sep=' ')
        if got.size != n:
            raise RuntimeError(f'Bad result shape for case {k}, n={n} got={got.shape}')
        max_mf=max(max_mf, discrepancy)
        func=function_factory(case.function)
        expected=np.asarray(func(case.t*case.matrix)@case.rhs, dtype=float)
        # SciPy 1.17 has a known incorrect slice on happy breakdown in a
        # restarted cycle; retain the independent dense reference even then.
        try:
            scipy = spla.funm_multiply_krylov(func,sp.csr_array(case.matrix),case.rhs,
                                             t=case.t,atol=0.0,rtol=1e-6,
                                             assume_a='hermitian' if case.hermitian else 'general',
                                             restart_every_m=case.m,max_restarts=case.max_restarts)
            scipy_cases += 1
            worst_scipy = max(worst_scipy, np.linalg.norm(scipy-expected)/max(1.0,np.linalg.norm(expected)))
        except ValueError as exc:
            if 'not aligned' not in str(exc):
                raise
            scipy_breakdown_bugs += 1
        err=np.linalg.norm(got-expected)
        relative=err/max(1.0,np.linalg.norm(expected))
        worst_abs=max(worst_abs,err)
        worst_rel=max(worst_rel,relative)
        if status == 0:
            converged_cases += 1
        if status > 0:
            unconverged.append((k,status,cycles,err,case.m,n))
        if status < 0 or not np.all(np.isfinite(got)) or discrepancy > 2e-10:
            raise AssertionError(f'case {k}: status={status} matrix-free diff={discrepancy}')
        if relative > (1e-7 if status == 0 else 2e-5):
            raise AssertionError(f'case {k}: relerr={relative:.3e} '
                                 f'status={status} m={case.m}, n={n}, kind={case.function}, symmetric={case.hermitian}')
        num_good+=1
    print(f'SciFort funm_multiply_krylov: {num_good}/{cases} dense SciPy comparisons PASS')
    print(f'Matrix-free/CSR max discrepancy: {max_mf:.3e}')
    print(f'Largest absolute error: {worst_abs:.3e}; scaled error: {worst_rel:.3e}')
    print(f'Largest scaled SciPy Krylov error on {scipy_cases} callable cases: {worst_scipy:.3e}')
    print(f'SciPy 1.17 restarted-breakdown exception cases: {scipy_breakdown_bugs}')
    print(f'Converged cases {converged_cases}; exceeded cycle budget in {len(unconverged)} cases')
    if unconverged:
        print('First unconverged:',unconverged[:5])


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('driver')
    p.add_argument('--seed',type=int,default=3381)
    p.add_argument('--cases',type=int,default=250)
    args=p.parse_args()
    check(args.driver,args.seed,args.cases)
