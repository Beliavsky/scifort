#!/usr/bin/env python3
"""Randomized tests for SciFort no-pivot ILUT versus SciPy SuperLU and dense LU.

Compare mathematical solutions, actual residuals, factor identity and triangular
structure, not identical fill/drop patterns (SciPy uses SuperLU's pivoting, area
rule and permutations). Original, independent Python numerical checks.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spilu as scipy_spilu

ROOT = Path(__file__).resolve().parents[2]


def make_cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    text = [str(count)]
    cases = []
    for case_id in range(count):
        n = int(rng.integers(3, 35))
        nrhs = 2
        density = float(rng.uniform(0.08, 0.27))
        a = rng.normal(scale=0.4, size=(n,n))
        a *= rng.random((n,n)) < density
        np.fill_diagonal(a,0.0)
        diag = np.sum(abs(a),axis=1) + rng.uniform(1.0,1.7,size=n)
        np.fill_diagonal(a,diag)
        mode = case_id % 5
        tol = (0.0,1e-4,0.003,0.03,0.12)[mode]
        fill = (100.0,8.0,3.0,1.5,1.0)[mode]
        shift = 0.05 if case_id % 7 == 0 else 0.0
        storage = case_id % 2
        b = rng.normal(size=(n,nrhs))
        original = sparse.coo_matrix(a)
        rows = original.row.tolist()
        cols = original.col.tolist()
        vals = original.data.tolist()
        rows.extend([0,0]); cols.extend([n-1,n-1]); vals.extend([0.25,-0.25])
        rows.append(0); cols.append(0); vals.append(0.0)
        text.append(f'{n} {len(rows)} {nrhs} {storage} {tol:.17g} {fill:.17g} {shift:.17g}')
        text.extend(f'{r+1} {c+1} {v:.17g}' for r,c,v in zip(rows,cols,vals))
        text.extend(f'{v:.17g}' for column in b.T for v in column)
        structural = sparse.coo_matrix((vals,(rows,cols)),shape=a.shape).tocsr()
        cases.append((a,b,tol,fill,shift,structural))
    return '\n'.join(text)+'\n', cases


def independent_ilut(a, tol: float, fill: float, shift: float) -> np.ndarray:
    """Pure Python sparse-dict ILUT reference with the documented row policies."""
    n = a.shape[0]
    factored = []
    for i in range(n):
        w = {int(c):float(v) for c,v in zip(
            a.indices[a.indptr[i]:a.indptr[i+1]],
            a.data[a.indptr[i]:a.indptr[i+1]])}
        scale = max((abs(v) for v in w.values()),default=0.)
        nlo = sum(c < i for c in w)
        nhi = sum(c > i for c in w)
        if shift:
            w[i] = w.get(i,0.)+shift
            if w[i] == 0.:
                del w[i]
        thresh = tol*scale
        for j in range(i):
            if j not in w:
                continue
            mult = w[j]/factored[j][j]
            if abs(mult) <= thresh:
                del w[j]
                continue
            w[j] = mult
            for col,value in sorted(factored[j].items()):
                if col <= j:
                    continue
                updated = w.get(col,0.)-mult*value
                if updated == 0.:
                    w.pop(col,None)
                else:
                    w[col] = updated
        if not w.get(i,0.):
            raise AssertionError(f'failed reference pivot on row {i}')
        for col in list(w):
            if col != i and abs(w[col]) <= thresh:
                del w[col]
        for cols, originally in (([c for c in w if c < i],nlo),
                                 ([c for c in w if c > i],nhi)):
            limit = min(n,int(min(float(n),fill*max(1,originally))))
            for col in sorted(cols,key=lambda c:(abs(w[c]),c))[:max(0,len(cols)-limit)]:
                del w[col]
        factored.append(w)
    result = np.zeros((n,n))
    for i,w in enumerate(factored):
        for col,value in w.items():
            result[i,col] = value
    return result


def parse_results(output: str, cases):
    lines = iter(output.splitlines())
    count_exact = 0
    count_dropped = 0
    count_scipy = 0
    count_gmres = 0
    peak_relative = 0.0
    for case_id, (a,b,tol,fill,shift,original_csr) in enumerate(cases,1):
        line = next(lines).split()
        assert line == ['CASE',str(case_id)], line
        line = next(lines).split()
        assert line == ['STATUS','0'], (case_id,line)
        nnz = int(next(lines).split()[1]); n = a.shape[0]
        fac = np.zeros_like(a)
        for _ in range(nnz):
            words = next(lines).split()
            assert words[0] == 'LU'
            fac[int(words[1])-1,int(words[2])-1] = float(words[3])
        assert np.all(np.isfinite(fac))
        independent = independent_ilut(original_csr,tol,fill,shift)
        assert np.allclose(fac,independent,atol=5e-12,rtol=5e-12),case_id
        lo = np.tril(fac,k=-1)+np.eye(n)
        hi = np.triu(fac)
        assert np.all(np.diag(hi) != 0)
        assert next(lines).split() == ['SOLVE','0']
        x = np.zeros_like(b)
        for col in range(b.shape[1]):
            words = next(lines).split()
            assert words[0] == 'X'
            x[:,col] = [float(v) for v in words[1:]]
        shifted = a+shift*np.eye(n)
        rhs_error = np.linalg.norm(shifted@x-b)/np.linalg.norm(b)
        factor_error = np.linalg.norm(lo@hi-shifted)/np.linalg.norm(shifted)
        exact_factor_solve = np.linalg.solve(lo@hi,b)
        assert np.allclose(x,exact_factor_solve,atol=5e-12,rtol=5e-12),case_id
        assert np.isfinite(rhs_error) and np.isfinite(factor_error)
        assert rhs_error < 0.60, (case_id,tol,fill,rhs_error)
        assert factor_error < 0.60, (case_id,tol,fill,factor_error)
        if tol == 0.0:
            assert rhs_error < 5e-12, (case_id,rhs_error)
            assert factor_error < 5e-12, (case_id,factor_error)
            count_exact += 1
        else:
            count_dropped += 1
        words = next(lines).split()
        assert words[0] == 'GMRES'
        info = int(words[1]); iterations = int(words[2]); residual = float(words[3])
        assert info == 0, (case_id,info,iterations)
        assert residual < 2e-9*np.linalg.norm(b[:,0]), (case_id,residual)
        count_gmres += 1
        peak_relative = max(peak_relative,rhs_error)
        try:
            ref = scipy_spilu(sparse.csc_matrix(shifted),drop_tol=tol,fill_factor=fill,
                permc_spec='NATURAL',diag_pivot_thresh=0.0,
                options={'Equil':False})
            ref_x = ref.solve(b)
            sci_resid = np.linalg.norm(shifted @ ref_x-b)/np.linalg.norm(b)
            assert np.isfinite(sci_resid)
            if tol == 0.0:
                assert sci_resid < 5e-12
            count_scipy += 1
        except RuntimeError:
            # SuperLU may reject a particular small fill budget; SciFort's
            # independent solution is still checked above.
            pass
    assert not list(lines)
    return count_exact,count_dropped,count_scipy,count_gmres,peak_relative


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=160)
    parser.add_argument('--seed',type=int,default=9841)
    args=parser.parse_args()
    input_text,cases=make_cases(args.cases,args.seed)
    for mode in ('checked','optimized') if args.mode=='both' else (args.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-spilu-ref-') as td:
            build=Path(td)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cmd=['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
            names=[('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                ('lsmr','src/optimize/lsmr.f90'),('sparse_linalg','src/sparse_linalg.f90')]
            for name,file in names:
                subprocess.run([*cmd,'-c',str(ROOT/file),'-o',str(build/(name+'.o'))],
                    check=True,cwd=build)
            exe=build/'spilu_reference'
            subprocess.run([*cmd,str(ROOT/'tools/sparse_linalg_reference/sparse_spilu_reference_driver.f90'),
                *[str(build/(name+'.o')) for name,_ in names],'-o',str(exe)],
                check=True,cwd=build)
            out=subprocess.run([str(exe)],input=input_text,encoding='utf8',
                capture_output=True,check=True,cwd=build).stdout
            exact,dropped,reference,gmres,peak=parse_results(out,cases)
            print(f'{mode}: {args.cases} ILUT cases passed, exact={exact}, dropped={dropped}, '
                f'SciPy SuperLU comparisons={reference}, GMRES={gmres}, '
                f'peak relative residual={peak:.3e}',flush=True)

if __name__=='__main__': main()
