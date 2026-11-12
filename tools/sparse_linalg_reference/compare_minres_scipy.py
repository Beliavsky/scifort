#!/usr/bin/env python3
"""Check SciFort MINRES and true residual against SciPy on symmetric CSR matrices."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import sparse
from scipy.sparse import linalg

ROOT = Path(__file__).resolve().parents[2]


def create_cases(seed, count):
    rng = np.random.default_rng(seed)
    stream = [str(count)]
    cases = []
    for case in range(count):
        n = int(rng.integers(2, 46))
        kind = case % 5
        density = float(rng.uniform(.08, .4))
        upper = np.triu((rng.random((n, n)) < density) * rng.uniform(-.45, .45, (n, n)), 1)
        dense = upper + upper.T
        if kind in (0, 1):
            dense[np.diag_indices(n)] = np.sum(np.abs(dense), axis=1) + rng.uniform(.9, 2.0, n)
        elif kind == 2:
            dense[np.diag_indices(n)] = -(np.sum(np.abs(dense), axis=1) + rng.uniform(.9, 2.0, n))
        elif kind == 3:
            dense[np.diag_indices(n)] = np.where(np.arange(n) % 2 == 0, 1, -1) * (
                np.sum(np.abs(dense), axis=1) + rng.uniform(.9, 2.0, n))
        else:
            dense[np.diag_indices(n)] = np.where(rng.random(n) < .5, 1, -1) * (
                np.sum(np.abs(dense), axis=1) + rng.uniform(.9, 2.0, n))
        shift = float(rng.uniform(-.3, .3)) if case % 3 == 0 else 0.0
        use_jacobi = int(kind in (0, 1) and case % 2 == 0)
        xtrue = rng.normal(size=n)
        b = (dense - shift * np.eye(n)) @ xtrue
        x0 = rng.normal(size=n) if case % 7 == 0 else np.zeros(n)
        coo = sparse.coo_matrix(dense)
        cap = 8*n
        rtol = 1e-10
        stream.append(f'{n} {coo.nnz} {shift:.17g} {rtol:.17g} {cap} {use_jacobi}')
        stream.extend(f'{i+1} {j+1} {v:.17g}' for i, j, v in zip(coo.row, coo.col, coo.data))
        stream.extend(f'{t:.17g}' for t in b)
        stream.extend(f'{t:.17g}' for t in x0)
        cases.append((dense, b, x0, shift, use_jacobi, xtrue, rtol))
    return '\n'.join(stream) + '\n', cases


def verify(output, cases):
    lines = iter(output.splitlines())
    totals = {'solved': 0, 'shifted': 0, 'preconditioned': 0}
    for i, (dense,b,x0,shift,jacobi,xtrue,rtol) in enumerate(cases,1):
        try:
            assert next(lines).strip() == f'CASE {i}'
            tokens = next(lines).split()
            assert tokens[0] == 'MINRES'
            status, steps = map(int,tokens[1:3])
            residual = float(tokens[3])
            tokens = next(lines).split()
            assert tokens[0] == 'X'
            x = np.asarray([float(t) for t in tokens[1:]])
            assert len(x) == len(b)
            A = sparse.csr_matrix(dense)
            M = sparse.diags(1/np.diag(dense),format='csr') if jacobi else None
            ref, info = linalg.minres(A,b,x0=x0,shift=shift,rtol=rtol,maxiter=8*len(b),M=M)
            assert info == 0, f'SciPy status {info}'
            assert status == 0, f'SciFort status {status}, residual={residual:.4g}'
            actual = np.linalg.norm(b-(dense-shift*np.eye(len(b)))@x)
            np.testing.assert_allclose(residual,actual,rtol=1e-7,atol=1e-11)
            assert actual <= max(2e-11,rtol*np.linalg.norm(b)), (actual,np.linalg.norm(b))
            np.testing.assert_allclose(x,xtrue,rtol=2e-7,atol=2e-7)
            # SciPy 1.17 computes its initial residual as b-A@x0, omitting
            # the shift*x0 term. With BOTH shift and nonzero x0 its returned
            # result can have an incorrect residual. Compare the mathematical
            # solution instead for that combination.
            if not (shift and np.any(x0)):
                np.testing.assert_allclose(x,ref,rtol=2e-6,atol=2e-6)
            assert steps <= 8*len(b)
            totals['solved'] += 1
            totals['shifted'] += bool(shift)
            totals['preconditioned'] += bool(jacobi)
        except (AssertionError, ValueError, StopIteration) as e:
            raise AssertionError(f'case {i}: {e}') from e
    assert next(lines,None) is None
    return totals


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cases',type=int,default=240)
    p.add_argument('--seed',type=int,default=20261105)
    p.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    opts = p.parse_args()
    stdin,cases = create_cases(opts.seed,opts.cases)
    for mode in (('checked','optimized') if opts.mode == 'both' else (opts.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-minres-ref-') as d:
            tmp=Path(d)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cmd=['gfortran','-std=f2018',*flags,'-J',str(tmp),'-I',str(tmp)]
            objects=[]
            for label,path in [('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                               ('lsmr','src/optimize/lsmr.f90'),('sparse_linalg','src/sparse_linalg.f90')]:
                obj=tmp/(label+'.o')
                subprocess.run([*cmd,'-c',str(ROOT/path),'-o',str(obj)],cwd=tmp,check=True)
                objects.append(str(obj))
            exe=tmp/'minres_reference'
            subprocess.run([*cmd,str(ROOT/'tools/sparse_linalg_reference/sparse_minres_reference_driver.f90'),
                           *objects,'-o',str(exe)],cwd=tmp,check=True)
            result=subprocess.run([str(exe)],input=stdin,text=True,capture_output=True,check=True,cwd=tmp)
            counts=verify(result.stdout,cases)
            print(f'{mode}: {counts}, seed={opts.seed}',flush=True)


if __name__=='__main__':
    main()
