#!/usr/bin/env python3
"""Cross-check complex sparse LU/spsolve against SciPy and dense residuals.

The Fortran driver processes a batch of text-formatted CSR/CSC cases.  All
numerical expectations are computed independently with SciPy SuperLU and
NumPy.  The default samples include unsorted, duplicate COO coordinates,
explicit zeros, row swaps, multiple right-hand sides, and N/T/H solves.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess

import numpy as np
from scipy.sparse import coo_array
from scipy.sparse.linalg import splu, spsolve


@dataclass
class Case:
    n: int
    rows: np.ndarray
    cols: np.ndarray
    values: np.ndarray
    b: np.ndarray
    fmt: int
    api: int
    trans: str
    pivot_tol: float = 0.0


def generate(seed: int, count: int) -> list[Case]:
    rng = np.random.default_rng(seed)
    cases: list[Case] = []
    for case_id in range(count):
        n = int(rng.integers(1, 26))
        nrhs = int(rng.integers(1, 5))
        mask = rng.random((n, n)) < rng.uniform(0.10, 0.35)
        dense = (rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))) * mask
        dense[np.diag_indices(n)] += n * 0.8 + 2.0 + 1j * 0.2
        if case_id % 3 == 0 and n > 1:
            dense = dense[rng.permutation(n)]
        rr, cc = np.nonzero(dense)
        vv = dense[rr, cc]
        # Split some stored elements into duplicate COO contributions.
        if len(vv) and case_id % 2 == 0:
            dup_idx = rng.choice(len(vv),size=max(1,len(vv)//4),replace=False)
            extra_rows = rr[dup_idx]
            extra_cols = cc[dup_idx]
            extra_values = vv[dup_idx]*0.35
            vv = vv.copy()
            vv[dup_idx] -= extra_values
            rr = np.concatenate([rr,extra_rows])
            cc = np.concatenate([cc,extra_cols])
            vv = np.concatenate([vv,extra_values])
        # Retain some explicit zeros not present in the matrix.
        if case_id % 4 == 0:
            rr = np.append(rr,int(rng.integers(0,n)))
            cc = np.append(cc,int(rng.integers(0,n)))
            vv = np.append(vv,0.0+0.0j)
        perm = rng.permutation(len(vv))
        rr,cc,vv = rr[perm],cc[perm],vv[perm]
        b = rng.normal(size=(n,nrhs))+1j*rng.normal(size=(n,nrhs))
        operation = 'N' if case_id % 4 == 0 else ('T' if case_id % 4 == 1 else 'H')
        api = 1 if operation == 'N' and case_id % 2 == 0 else 0
        fmt = (case_id % 2) + 1
        cases.append(Case(n,rr,cc,vv,b,fmt,api,operation))
    return cases


def run(driver: Path, cases: list[Case]) -> None:
    lines: list[str] = []
    for case in cases:
        codes = {'N':1,'T':2,'H':3}
        lines.append(f'{case.n} {len(case.values)} {case.b.shape[1]} '
                     f'{case.fmt} {case.api} {codes[case.trans]} {case.pivot_tol:.17g}')
        for r,c,v in zip(case.rows,case.cols,case.values):
            lines.append(f'{r+1} {c+1} {v.real:.17g} {v.imag:.17g}')
        for z in case.b.reshape(-1,order='F'):
            lines.append(f'{z.real:.17g} {z.imag:.17g}')
    output = subprocess.run([str(driver.resolve())],input='\n'.join(lines)+'\n',
                            capture_output=True,text=True,check=True)
    output_lines = iter(output.stdout.splitlines())
    biggest_error = 0.0
    for case_id,case in enumerate(cases):
        try:
            status = int(next(output_lines))
        except (ValueError, StopIteration) as exc:
            raise AssertionError(f'Case {case_id}: invalid driver output') from exc
        if status != 0:
            raise AssertionError(f'Case {case_id}: factorization status={status}')
        try:
            received = np.array([complex(*map(float,next(output_lines).split()))
                                 for _ in range(case.n*case.b.shape[1])])
        except StopIteration as exc:
            raise AssertionError(f'Case {case_id}: incomplete output') from exc
        received = received.reshape(case.b.shape,order='F')
        a = coo_array((case.values,(case.rows,case.cols)),shape=(case.n,case.n)).tocsc()
        a.sum_duplicates()
        if case.api == 1:
            expected = np.asarray(spsolve(a,case.b)).reshape(case.b.shape)
        else:
            expected = splu(a,permc_spec='NATURAL').solve(case.b,trans=case.trans)
        err = float(np.max(np.abs(received-expected)))
        res = (a.toarray() if case.trans == 'N' else
               a.toarray().T if case.trans == 'T' else a.toarray().conj().T)
        resid = float(np.max(np.abs(res@received-case.b)))
        scale = max(1.0,float(np.max(np.abs(expected))))
        if err > 2.0e-10*scale or resid > 2.0e-10*max(1.0,float(np.max(np.abs(case.b)))):
            raise AssertionError(f'Case {case_id}, n={case.n}, trans={case.trans}, '
                                 f'fmt={case.fmt}, API={case.api}: err={err:.3e}, resid={resid:.3e}')
        biggest_error = max(biggest_error,err)
    if next(output_lines,None) is not None:
        raise AssertionError('Unexpected extra output from driver')
    print(f'PASS {len(cases)} complex direct sparse SciPy cases; '
          f'max absolute solution difference {biggest_error:.3e}',flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver',type=Path)
    parser.add_argument('--count',type=int,default=250)
    parser.add_argument('--seed',type=int,default=2667)
    args = parser.parse_args()
    if args.count <= 0:
        parser.error('--count must be positive')
    run(args.driver,generate(args.seed,args.count))


if __name__ == '__main__':
    main()
