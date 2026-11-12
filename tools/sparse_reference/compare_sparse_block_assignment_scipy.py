#!/usr/bin/env python3
"""Randomized SciPy checks for canonical CSR block assembly and sparse scatter.

Missing blocks are `None` in SciPy and `csr_matrix()` in SciFort.  Explicit
row/column sizes in the Fortran reference test also handle all-None rows and
columns that SciPy cannot infer without materialized sparse zero blocks.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]


@dataclass
class Case:
    op: int
    initial: sparse.csr_matrix | None = None
    grid: list[list[sparse.csr_matrix | None]] | None = None
    heights: list[int] | None = None
    widths: list[int] | None = None
    coord_updates: list[tuple[int, int, float]] | None = None
    rows: list[int] | None = None
    cols: list[int] | None = None
    values: np.ndarray | float | sparse.csr_matrix | None = None
    one: tuple[int, int, float] | None = None


def random_csr(rng: np.random.Generator, nr: int, nc: int) -> sparse.csr_matrix:
    count = int(rng.integers(0, nr*nc+5)) if nr*nc else 0
    rows = rng.integers(nr,size=count)
    cols = rng.integers(nc,size=count)
    vals = rng.integers(-9,10,size=count).astype(float)/4
    if count:
        vals[:min(3,count)] = 0.0
    return sparse.coo_matrix((vals,(rows,cols)),shape=(nr,nc)).tocsr()


def generate(rng: np.random.Generator, n: int) -> list[Case]:
    all_cases = []
    for index in range(n):
        op = index % 6 + 1
        if op == 1:
            nr = int(rng.integers(1,5))
            nc = int(rng.integers(1,5))
            heights = list(map(int,rng.integers(0,6,size=nr)))
            widths = list(map(int,rng.integers(0,6,size=nc)))
            grid = []
            for i in range(nr):
                grid.append([random_csr(rng,heights[i],widths[j]) if rng.random()<0.55
                             else None for j in range(nc)])
            all_cases.append(Case(op,grid=grid,heights=heights,widths=widths))
            continue
        m = int(rng.integers(0,9))
        k = int(rng.integers(0,9))
        matrix = random_csr(rng,m,k)
        if op == 2:
            nupdates = int(rng.integers(0,34)) if m*k else 0
            updates = []
            for _ in range(nupdates):
                r = int(rng.integers(m))
                c = int(rng.integers(k))
                v = float(rng.integers(-15,16))/4
                updates.append((r,c,v))
            all_cases.append(Case(op,initial=matrix,coord_updates=updates))
        elif op in (3,4,5):
            nr = int(rng.integers(0,9)) if m else 0
            nc = int(rng.integers(0,9)) if k else 0
            rows = list(map(int,rng.integers(m,size=nr))) if m else []
            cols = list(map(int,rng.integers(k,size=nc))) if k else []
            if op == 3:
                values = rng.integers(-9,10,size=(nr,nc)).astype(float)/4
            elif op == 4:
                values = float(rng.integers(-9,10))/4
            else:
                values = random_csr(rng,nr,nc)
            all_cases.append(Case(op,initial=matrix,rows=rows,cols=cols,values=values))
        else:
            if m*k:
                r,c = int(rng.integers(m)),int(rng.integers(k))
                one = (r,c,float(rng.integers(-9,10))/4)
            else:
                one = (-1,0,1.0)
            all_cases.append(Case(op,initial=matrix,one=one))
    return all_cases


def serialize_csr(buf: StringIO, a: sparse.csr_matrix) -> None:
    coo = a.tocoo()
    buf.write(f'{a.shape[0]} {a.shape[1]} {len(coo.data)}\n')
    for r,c,v in zip(coo.row,coo.col,coo.data):
        buf.write(f'{r+1} {c+1} {float(v):.17g}\n')


def serialize(cases: list[Case]) -> str:
    buf = StringIO()
    buf.write(f'{len(cases)}\n')
    for case in cases:
        buf.write(f'{case.op}\n')
        if case.op == 1:
            assert case.grid is not None and case.heights is not None and case.widths is not None
            nr,nc = len(case.heights),len(case.widths)
            buf.write(f'{nr} {nc}\n')
            for h in case.heights: buf.write(f'{h}\n')
            for w in case.widths: buf.write(f'{w}\n')
            for j in range(nc):
                for i in range(nr):
                    block = case.grid[i][j]
                    buf.write(f'{int(block is not None)}\n')
                    if block is not None: serialize_csr(buf,block)
            continue
        assert case.initial is not None
        serialize_csr(buf,case.initial)
        if case.op == 2:
            updates = case.coord_updates or []
            buf.write(f'{len(updates)}\n')
            for r,c,v in updates:
                buf.write(f'{r+1} {c+1} {v:.17g}\n')
        elif case.op in (3,4,5):
            rows,cols = case.rows or [],case.cols or []
            buf.write(f'{len(rows)} {len(cols)}\n')
            for r in rows: buf.write(f'{r+1}\n')
            for c in cols: buf.write(f'{c+1}\n')
            if case.op == 3:
                values = np.asarray(case.values)
                for j in range(len(cols)):
                    for i in range(len(rows)):
                        buf.write(f'{float(values[i,j]):.17g}\n')
            elif case.op == 4:
                buf.write(f'{float(case.values):.17g}\n')
            else:
                assert sparse.issparse(case.values)
                serialize_csr(buf,case.values.tocsr())
        else:
            assert case.one is not None
            r,c,v = case.one
            buf.write(f'{r+1} {c+1} {v:.17g}\n')
    return buf.getvalue()


def expected(case: Case) -> np.ndarray | None:
    if case.op == 1:
        assert case.grid is not None and case.heights is not None and case.widths is not None
        # Use SciPy's general block assembler.  Rows/cols with no occupied
        # blocks require explicit empty sparse blocks to specify shape.
        filled = [[case.grid[i][j] if case.grid[i][j] is not None else
                   sparse.csr_matrix((case.heights[i],case.widths[j]))
                   for j in range(len(case.widths))] for i in range(len(case.heights))]
        return sparse.bmat(filled,format='csr').toarray()
    assert case.initial is not None
    result = case.initial.toarray()
    if case.op == 2:
        for r,c,v in case.coord_updates or []:
            result[r,c] = v
    elif case.op in (3,4,5):
        rows,cols = case.rows or [],case.cols or []
        if case.op == 3:
            vals = np.asarray(case.values)
        elif case.op == 4:
            vals = np.full((len(rows),len(cols)),float(case.values))
        else:
            assert sparse.issparse(case.values)
            vals = case.values.toarray()
        # Explicit Fortran array-element assignment order is column-major.
        for j,c in enumerate(cols):
            for i,r in enumerate(rows):
                result[r,c] = vals[i,j]
    elif case.op == 6:
        assert case.one is not None
        r,c,v = case.one
        if r<0 or c<0 or r>=result.shape[0] or c>=result.shape[1]:
            return None  # expected invalid, unchanged CSR
        result[r,c] = v
    return result


def compare(exe: Path, cases: list[Case]) -> None:
    proc = subprocess.run([str(exe)],input=serialize(cases),text=True,
                          capture_output=True,check=True)
    lines = iter(proc.stdout.splitlines())
    for i,case in enumerate(cases):
        info,m,n,nnz = map(int,next(lines).split())
        truth = expected(case)
        if truth is None:
            if info != 1: raise AssertionError(f'case {i}: expected invalid coordinate, got {info}')
            continue
        if info != 0: raise AssertionError(f'case {i}, op {case.op}: status {info}')
        r,c,values = [],[],[]
        for _ in range(nnz):
            rr,cc,val = next(lines).split()
            r.append(int(rr)-1); c.append(int(cc)-1); values.append(float(val))
        coords = list(zip(r,c))
        if coords != sorted(coords) or len(coords) != len(set(coords)):
            raise AssertionError(f'case {i}: noncanonical coordinates {coords}')
        got = sparse.coo_matrix((values,(r,c)),shape=(m,n)).toarray()
        if truth.shape != got.shape or not np.allclose(truth,got,rtol=0,atol=2e-13):
            raise AssertionError(f'case {i}, op {case.op}: shapes {truth.shape},{got.shape}\n{truth}\n{got}')
        if case.op == 1:
            assert case.grid is not None
            exact_nnz = sum(x.nnz for row in case.grid for x in row if x is not None)
            if nnz != exact_nnz:
                raise AssertionError(f'case {i}: stored zeros lost in block assembly: {nnz} vs {exact_nnz}')
    if next(lines,None) is not None:
        raise AssertionError('extra Fortran output')


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    p.add_argument('--cases',type=int,default=420)
    p.add_argument('--seed',type=int,default=20261105)
    args = p.parse_args()
    cases = generate(np.random.default_rng(args.seed),args.cases)
    with tempfile.TemporaryDirectory(prefix='scifort-block-assignment-') as td:
        build = Path(td)
        for mode in (('checked','optimized') if args.mode == 'both' else (args.mode,)):
            flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc = ['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
            subprocess.run([*cc,'-c',str(ROOT/'src/kinds.f90'),'-o',str(build/'kinds.o')],check=True)
            subprocess.run([*cc,'-c',str(ROOT/'src/sparse.f90'),'-o',str(build/'sparse.o')],check=True)
            exe = build/'sparse_block_assignment_reference'
            subprocess.run([*cc,str(ROOT/'tools/sparse_reference/sparse_block_assignment_reference_driver.f90'),
                            str(build/'kinds.o'),str(build/'sparse.o'),'-o',str(exe)],check=True)
            compare(exe,cases)
            print(f'{mode}: {len(cases)} SciPy/NumPy block/assignment comparisons passed (seed {args.seed})',flush=True)


if __name__ == '__main__':
    main()
