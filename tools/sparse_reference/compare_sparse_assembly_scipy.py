#!/usr/bin/env python3
"""Randomized CSR assembly/indexing comparisons with SciPy (Fortran indices 1-based).

Run both native builds; compare true sparse numerical values and canonical
storage. All nonzero input coordinates, including deliberately stored zeros,
are passed as COO records to the Fortran reference driver.
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
    mats: list[sparse.csr_matrix]
    rows: list[int] | None = None
    cols: list[int] | None = None
    bounds: tuple[int, int, int, int, int, int] | None = None


def make_matrix(rng: np.random.Generator, m: int, n: int) -> sparse.csr_matrix:
    count = int(rng.integers(0, m*n + 7)) if m*n else 0
    rows = rng.integers(m, size=count)
    cols = rng.integers(n, size=count)
    data = rng.integers(-8, 9, size=count).astype(float) / 4
    if len(data) > 2:
        data[:2] = 0.0  # preserve a few explicit zero coordinates
    return sparse.coo_matrix((data, (rows, cols)), shape=(m, n)).tocsr()


def generate(rng: np.random.Generator, cases: int) -> list[Case]:
    result: list[Case] = []
    for i in range(cases):
        op = i % 6 + 1
        if op <= 3:
            count = int(rng.integers(1, 5)) if op != 3 else int(rng.integers(0, 5))
            same = int(rng.integers(0, 9))
            matrices: list[sparse.csr_matrix] = []
            for j in range(count):
                m = same if op == 1 else int(rng.integers(0, 9))
                n = same if op == 2 else int(rng.integers(0, 9))
                matrices.append(make_matrix(rng, m, n))
            result.append(Case(op, matrices))
        elif op == 4:
            dims = rng.integers(0, 8, size=4)
            result.append(Case(op, [make_matrix(rng, *map(int, dims[:2])),
                                    make_matrix(rng, *map(int, dims[2:]))]))
        elif op == 5:
            m = int(rng.integers(0, 10))
            n = int(rng.integers(0, 10))
            rows = list(map(int, rng.integers(m, size=int(rng.integers(0, 12))))) if m else []
            cols = list(map(int, rng.integers(n, size=int(rng.integers(0, 12))))) if n else []
            result.append(Case(op,[make_matrix(rng,m,n)],rows,cols))
        else:
            m = int(rng.integers(0, 10))
            n = int(rng.integers(0, 10))
            rs = int(rng.choice([-3,-2,-1,1,2,3]))
            cs = int(rng.choice([-3,-2,-1,1,2,3]))
            r0,r1 = valid_slice_bounds(rng,m,rs)
            c0,c1 = valid_slice_bounds(rng,n,cs)
            result.append(Case(op,[make_matrix(rng,m,n)],bounds=(r0,r1,c0,c1,rs,cs)))
    return result


def valid_slice_bounds(rng: np.random.Generator, n: int, step: int) -> tuple[int,int]:
    if step > 0:
        start = int(rng.integers(1,n+2))
        stop = int(rng.integers(0,n+1))
    else:
        start = int(rng.integers(0,n+1))
        stop = int(rng.integers(1,n+2))
    return start,stop


def slice_positions(first: int, last: int, step: int) -> list[int]:
    return list(range(first, last+(1 if step > 0 else -1), step))


def expected(case: Case) -> sparse.csr_matrix:
    mats = case.mats
    if case.op == 1:
        return sparse.hstack(mats, format='csr')
    if case.op == 2:
        return sparse.vstack(mats, format='csr')
    if case.op == 3:
        return sparse.block_diag(mats,format='csr') if mats else sparse.csr_matrix((0, 0))
    if case.op == 4:
        return sparse.kron(*mats,format='csr')
    if case.op == 5:
        return mats[0][case.rows, :][:, case.cols].tocsr()
    r0,r1,c0,c1,rs,cs = case.bounds
    rows = [r-1 for r in slice_positions(r0,r1,rs)]
    cols = [c-1 for c in slice_positions(c0,c1,cs)]
    return mats[0][rows, :][:, cols].tocsr()


def prepare_input(cases: list[Case]) -> str:
    buf = StringIO()
    buf.write(f'{len(cases)}\n')
    for case in cases:
        buf.write(f'{case.op}\n')
        if case.op <= 3:
            buf.write(f'{len(case.mats)}\n')
        for a in case.mats:
            coo = a.tocoo()
            buf.write(f'{a.shape[0]} {a.shape[1]} {coo.nnz}\n')
            for r,c,v in zip(coo.row,coo.col,coo.data):
                buf.write(f'{r+1} {c+1} {float(v):.17g}\n')
        if case.op == 5:
            rows = case.rows or []
            cols = case.cols or []
            buf.write(f'{len(rows)} {len(cols)}\n')
            for r in rows:
                buf.write(f'{r+1}\n')
            for c in cols:
                buf.write(f'{c+1}\n')
        if case.op == 6:
            buf.write(' '.join(str(x) for x in case.bounds)+'\n')
    return buf.getvalue()


def compare(exe: Path, cases: list[Case]) -> None:
    proc = subprocess.run([str(exe)], input=prepare_input(cases), text=True,
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
    lines = iter(proc.stdout.splitlines())
    for i,case in enumerate(cases):
        try:
            info,m,n,nnz = map(int,next(lines).split())
        except (StopIteration,ValueError) as exc:
            raise AssertionError(f'bad header on case {i}: {case.op}') from exc
        if info != 0:
            raise AssertionError(f'Fortran status {info} on case {i}: {case}')
        row,col,vals = [],[],[]
        for _ in range(nnz):
            r,c,v = next(lines).split()
            row.append(int(r)-1); col.append(int(c)-1); vals.append(float(v))
        out = sparse.coo_matrix((vals,(row,col)),shape=(m,n)).tocsr()
        truth = expected(case)
        if (m,n) != truth.shape:
            raise AssertionError(f'shape mismatch case {i}: got {(m,n)}, expected {truth.shape}')
        diff = out.toarray()-truth.toarray()
        scale = max(1.0,float(np.max(abs(truth.data),initial=0.0)))
        if np.max(abs(diff),initial=0.0) > 2e-12*scale:
            raise AssertionError(f'value mismatch case {i}: op {case.op}\ngot {out.toarray()}\nexpected {truth.toarray()}')
        if row:
            coords = list(zip(row,col))
            if len(coords) != len(set(coords)) or coords != sorted(coords):
                raise AssertionError(f'noncanonical output case {i}: {coords}')
    if next(lines,None) is not None:
        raise AssertionError('extra Fortran output records')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=420)
    parser.add_argument('--seed',type=int,default=13579)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    cases = generate(rng,args.cases)
    modes = ('checked','optimized') if args.mode=='both' else (args.mode,)
    with tempfile.TemporaryDirectory(prefix='scifort-sparse-assembly-') as td:
        build = Path(td)
        for mode in modes:
            flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc = ['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
            subprocess.run([*cc,'-c',str(ROOT/'src/kinds.f90'),'-o',str(build/'kinds.o')],check=True)
            subprocess.run([*cc,'-c',str(ROOT/'src/sparse.f90'),'-o',str(build/'sparse.o')],check=True)
            exe = build/'assembly_ref'
            subprocess.run([*cc,str(ROOT/'tools/sparse_reference/sparse_assembly_reference_driver.f90'),
                            str(build/'kinds.o'),str(build/'sparse.o'),'-o',str(exe)],check=True)
            compare(exe,cases)
            print(f'{mode}: {len(cases)} randomized assembly/indexing comparisons passed (seed {args.seed})',flush=True)

if __name__ == '__main__':
    main()
