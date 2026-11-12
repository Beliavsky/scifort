#!/usr/bin/env python3
"""Cross-check SciFort CSC conversion, diagonal and triangular operations with SciPy.

All numerical reference results are independently computed with SciPy 1.x.
Zero-size, rectangular, duplicate, stored-zero, and diagonal-offset cases are included.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
LABELS = ("A", "C", "R", "T", "B", "U", "L", "D", "S")


def diagonal_bounds(m: int, n: int, offset: int) -> tuple[int, int]:
    start = max(0, -offset)
    stop = min(m, n-offset)
    return start, max(start, stop)


def diagonal_len(m: int, n: int, offset: int) -> int:
    start, stop = diagonal_bounds(m, n, offset)
    return stop-start


def form_cases(rng: np.random.Generator, count: int):
    lines = [str(count)]
    expectations = []
    for case in range(count):
        m = int(rng.integers(0, 22))
        n = int(rng.integers(0, 22))
        if case % 17 == 0:
            m, n = 0, 0
        elif case % 17 == 1:
            m, n = 16, 0
        elif case % 17 == 2:
            m, n = 0, 14
        elif case % 17 == 3:
            m, n = 2, 20
        elif case % 17 == 4:
            m, n = 20, 2
        nnz = int(rng.integers(0, m*n//2+1))
        rows = rng.integers(0, m, size=nnz)
        cols = rng.integers(0, n, size=nnz)
        values = rng.integers(-12, 13, size=nnz).astype(float)/4
        if nnz:
            values[rng.random(nnz) < .15] = 0
        coo = sparse.coo_matrix((values,(rows,cols)),shape=(m,n))
        csr = coo.tocsr()
        csr.sum_duplicates()
        dense = csr.toarray()
        offsets = [-2, 0, 2]
        if case % 5 == 0:
            offsets = [0, -3, 1, 9]
        diag_ptr = [1]
        diags = []
        for off in offsets:
            d = rng.integers(-9, 10, size=diagonal_len(m,n,off)).astype(float)/4
            if d.size and rng.random() < .4:
                d[int(rng.integers(0,len(d)))] = 0.0
            diags.append(d)
            diag_ptr.append(diag_ptr[-1]+len(d))
        k_cut = int(rng.integers(-25,26))
        k_set = int(rng.integers(-25,26))
        k_diag = int(rng.integers(-25,26))
        d_set = rng.integers(-8,9,size=diagonal_len(m,n,k_set)).astype(float)/2
        if d_set.size and case % 3 == 0:
            d_set[0] = 0
        lines.append(f"{m} {n} {nnz}")
        lines.extend(f"{r+1} {c+1} {v:.17g}" for r,c,v in zip(rows,cols,values))
        lines.append(str(len(offsets)))
        lines.append(' '.join(map(str,offsets)))
        lines.append(' '.join(map(str,diag_ptr)))
        lines.extend(f"{v:.17g}" for d in diags for v in d)
        lines.append(f"{k_cut} {k_set} {k_diag}")
        lines.append(str(len(d_set)))
        lines.extend(f"{v:.17g}" for v in d_set)

        # SciPy diags can reject offset diagonals of zero length for some
        # dimension/offset combinations. Assemble the same matrix explicitly,
        # then also verify with scipy.sparse.diags for compatible shapes.
        diags_dense = np.zeros((m,n))
        for off,d in zip(offsets,diags):
            start, stop = diagonal_bounds(m,n,off)
            for i,v in zip(range(start,stop),d):
                diags_dense[i,i+off] = v
        if m>0 and n>0 and all(diagonal_len(m,n,off)>0 for off in offsets):
            scipy_diags = sparse.diags(diags,offsets,shape=(m,n),format='csr').toarray()
            np.testing.assert_allclose(diags_dense, scipy_diags)

        set_result = dense.copy()
        start, stop = diagonal_bounds(m,n,k_set)
        for i,v in zip(range(start,stop),d_set):
            set_result[i,i+k_set] = v
        csr_set = csr.copy()
        # SciPy cannot safely set an empty offset on all zero-size shapes.
        if m>0 and n>0 and d_set.size:
            csr_set.setdiag(d_set,k=k_set)
            np.testing.assert_allclose(set_result,csr_set.toarray())
        expected = {
            'A':dense,
            'C':csr.tocsc().toarray(),
            'R':csr.tocsc().tocsr().toarray(),
            'T':csr.T.tocsc().toarray(),
            'B':diags_dense,
            'U':sparse.triu(csr,k=k_cut,format='csr').toarray(),
            'L':sparse.tril(csr,k=k_cut,format='csr').toarray(),
            'D':np.array([csr.diagonal(k=k_diag)]),
            'S':set_result,
        }
        expectations.append(expected)
    return '\n'.join(lines)+'\n',expectations


def verify(output: str, expected_cases) -> None:
    lines = iter(output.splitlines())
    for case_id, expected in enumerate(expected_cases, start=1):
        headline = next(lines)
        assert headline == f'CASE {case_id}', (case_id, headline)
        for label in LABELS:
            header = next(lines).split()
            assert header[0] == label, (case_id,label,header)
            m,n = map(int,header[1:])
            result = np.empty((m,n))
            for r in range(m):
                line = next(lines).strip()
                vals = [float(v) for v in line.split()]
                assert len(vals) == n, (case_id,label,len(vals),n)
                result[r] = vals
            exp = expected[label]
            assert result.shape == exp.shape,(case_id,label,result.shape,exp.shape)
            try:
                np.testing.assert_allclose(result,exp,rtol=1e-12,atol=1e-12)
            except AssertionError as err:
                raise AssertionError(f'case={case_id} operation={label}: {err}') from err
    assert next(lines,None) is None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('both','checked','optimized'),default='both')
    parser.add_argument('--cases',type=int,default=350)
    parser.add_argument('--seed',type=int,default=20261104)
    args=parser.parse_args()
    cases,expected=form_cases(np.random.default_rng(args.seed),args.cases)
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-sparse-csc-') as name:
            build=Path(name)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc=['gfortran','-std=f2018','-ffree-line-length-none',*flags,'-J',str(build),'-I',str(build)]
            objects=[]
            for basename in ('kinds','sparse'):
                obj=build/(basename+'.o')
                subprocess.run([*cc,'-c',str(ROOT/'src'/(basename+'.f90')),'-o',str(obj)],cwd=build,check=True)
                objects.append(obj)
            exe=build/'sparse_csc_reference'
            subprocess.run([*cc,str(ROOT/'tools/sparse_reference/sparse_csc_reference_driver.f90'),
                            *(str(obj) for obj in objects),'-o',str(exe)],cwd=build,check=True)
            p=subprocess.run([str(exe)],input=cases,text=True,capture_output=True,cwd=build)
            if p.returncode:
                raise RuntimeError(f'{mode}: Fortran failure ({p.returncode})\nSTDOUT\n{p.stdout[-3000:]}\nSTDERR\n{p.stderr[-3000:]}')
            verify(p.stdout,expected)
            print(f'{mode}: {args.cases} randomized CSC and diagonal cases PASS '
                  f'({args.cases*len(LABELS)} matrix comparisons)')

if __name__=='__main__':
    main()
