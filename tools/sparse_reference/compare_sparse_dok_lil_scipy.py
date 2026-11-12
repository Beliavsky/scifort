#!/usr/bin/env python3
"""Independent randomized DOK/LIL comparisons against SciPy 1.17+.

Builds a dedicated Fortran driver in both bounds-checked and optimized modes.
The input includes duplicate COO entries, explicit zeros, and hundreds of
randomized assignment, addition, deletion, and lookup sequences.  Both
numerical values and stored-entry counts are checked at every step.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from numpy.testing import assert_allclose
from scipy.sparse import coo_array, dok_array, lil_array

ROOT = Path(__file__).resolve().parents[2]


@dataclass
class Case:
    m: int
    n: int
    row: np.ndarray
    col: np.ndarray
    val: np.ndarray
    x: np.ndarray
    ops: list[tuple[int, int, int, float]]


def make_cases(ncases: int, seed: int) -> list[Case]:
    rng = np.random.default_rng(seed)
    cases = []
    shapes = [(0,0), (0,7), (7,0), (1,1), (1,9), (9,1),
              (2,3), (8,8), (3,11), (11,2)]
    for it in range(ncases):
        m, n = shapes[it % len(shapes)] if it < 2*len(shapes) else (
            int(rng.integers(0,13)), int(rng.integers(0,13)))
        nnz = int(rng.integers(0,41)) if m and n else 0
        row = rng.integers(0,m,size=nnz).astype(int) if nnz else np.array([],dtype=int)
        col = rng.integers(0,n,size=nnz).astype(int) if nnz else np.array([],dtype=int)
        val = rng.integers(-5,6,size=nnz).astype(float)
        if nnz > 0:
            val[::4] = 0.0
            row[::5] = row[0]
            col[::5] = col[0]
        x = rng.integers(-4,5,size=n).astype(float)
        operations = []
        if m and n:
            for _ in range(int(rng.integers(12,60))):
                i = int(rng.integers(0,m))
                j = int(rng.integers(0,n))
                op = int(rng.integers(1,5))
                v = float(rng.integers(-4,5))
                if rng.integers(0,4) == 0:
                    v = 0.0
                operations.append((op,i,j,v))
        cases.append(Case(m,n,row,col,val,x,operations))
    return cases


def case_input(cases: list[Case]) -> str:
    lines = [str(len(cases))]
    for c in cases:
        lines.append(f'{c.m} {c.n} {len(c.val)} {len(c.ops)}')
        lines.extend(f'{i+1} {j+1} {v:.17g}'
                     for i,j,v in zip(c.row,c.col,c.val))
        if c.n:
            lines.append(' '.join(f'{v:.17g}' for v in c.x))
        lines.extend(f'{op} {i+1} {j+1} {v:.17g}' for op,i,j,v in c.ops)
    return '\n'.join(lines)+'\n'


def compile_driver(build: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    command = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    objects = []
    for name in ('kinds','sparse','sparse_coo','sparse_dok_lil'):
        source = ROOT/'src'/f'{name}.f90'
        out = build/f'{name}.o'
        subprocess.run([*command,'-c',str(source),'-o',str(out)],cwd=build,check=True)
        objects.append(str(out))
    binary = build/'dok_lil_reference'
    subprocess.run([*command, str(ROOT/'tools/sparse_reference/sparse_dok_lil_reference_driver.f90'),
                    *objects,'-o',str(binary)],cwd=build,check=True)
    return binary


def compare(mode: str, cases: list[Case]) -> None:
    with tempfile.TemporaryDirectory(prefix=f'scifort-dok-lil-{mode}-') as dirname:
        exe = compile_driver(Path(dirname), mode)
        completed = subprocess.run([str(exe)],input=case_input(cases),
                                   text=True,capture_output=True,check=True)
        lines = iter(completed.stdout.splitlines())
        for icase,c in enumerate(cases,1):
            initial = coo_array((c.val,(c.row,c.col)),shape=(c.m,c.n))
            matrices = [dok_array(initial),lil_array(initial)]
            for k,(op,i,j,v) in enumerate(c.ops,1):
                if op == 1:
                    for a in matrices:
                        a[i,j] = v
                elif op == 2:
                    for a in matrices:
                        a[i,j] = float(a[i,j])+v
                elif op == 3:
                    for a in matrices:
                        a[i,j] = 0.0
                tag, got_icase, got_step, vd, vl, nd, nl = next(lines).split()
                assert (tag,int(got_icase),int(got_step)) == ('EDIT',icase,k)
                assert_allclose([float(vd),float(vl)],
                                [float(a[i,j]) for a in matrices],atol=2e-12,rtol=0)
                assert (int(nd),int(nl)) == tuple(a.nnz for a in matrices), (
                    mode,icase,k,'nnz',int(nd),int(nl),[a.nnz for a in matrices])
            parts = next(lines).split()
            assert parts[0] == 'CASE'
            assert tuple(map(int,parts[1:])) == (
                icase,c.m,c.n,matrices[0].nnz,matrices[1].nnz), (mode,icase,parts)
            dense = [a.toarray() for a in matrices]
            product = [a@c.x for a in matrices]
            for i in range(c.m):
                for j in range(c.n):
                    tag,ri,cj,vd,vl = next(lines).split()
                    assert (tag,int(ri),int(cj)) == ('D',i+1,j+1)
                    assert_allclose([float(vd),float(vl)],
                                    [d[i,j] for d in dense],atol=2e-12,rtol=0)
                tag,ri,vd,vl = next(lines).split()
                assert (tag,int(ri)) == ('Y',i+1)
                assert_allclose([float(vd),float(vl)],
                                [z[i] for z in product],atol=1e-10,rtol=0)
            for name, a, suff in [('DOK',matrices[0],'D'),('LIL',matrices[1],'L')]:
                csr = a.tocsr()
                csr.sort_indices()
                for k in range(a.nnz):
                    tag, idx, value = next(lines).split()
                    assert tag == f'CSR{suff}'
                    assert int(idx) == csr.indices[k]+1, (mode,icase,name,'index',k)
                    assert_allclose(float(value),csr.data[k],atol=1e-11,rtol=0)
                for k in range(c.m+1):
                    tag, ptr = next(lines).split()
                    assert tag == f'PTR{suff}'
                    assert int(ptr) == csr.indptr[k]+1, (mode,icase,name,'ptr',k)
        assert next(lines,None) is None, f'Unexpected extra driver output for {mode}'
        print(f'{mode}: {len(cases)} randomized SciPy DOK/LIL cases: PASS', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('checked','optimized','both'),default='both')
    parser.add_argument('--trials',type=int,default=250)
    parser.add_argument('--seed',type=int,default=20261105)
    args = parser.parse_args()
    cases = make_cases(args.trials,args.seed)
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        compare(mode,cases)


if __name__ == '__main__':
    main()
