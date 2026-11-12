#!/usr/bin/env python3
"""Compare native CSR deletion, reshape and sparse scatter against SciPy 1.17+.

Each test uses one-based Fortran indices and zero-based NumPy reference indices.
Run with --mode both and --cases 250 to cover both build configurations.
"""
from __future__ import annotations

import argparse
import io
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]


def random_sparse(rng: np.random.Generator, m: int, n: int) -> sparse.csr_matrix:
    nz = int(rng.integers(0, m*n + 5)) if m*n else 0
    r = rng.integers(0,m,size=nz)
    c = rng.integers(0,n,size=nz)
    x = rng.integers(-12,13,size=nz).astype(float) / 4
    if nz:
        x[:min(3,nz)] = 0.0
    return sparse.coo_matrix((x,(r,c)),shape=(m,n)).tocsr()


def serialize(buf: io.StringIO, matrix: sparse.csr_matrix) -> None:
    coo = matrix.tocoo()
    for r,c,v in zip(coo.row,coo.col,coo.data):
        buf.write(f'{r+1} {c+1} {float(v):.17g}\n')


def generate(rng: np.random.Generator, cases: int):
    buf = io.StringIO()
    expected = []
    for index in range(cases):
        op = index % 5 + 1
        m = int(rng.integers(0, 13))
        n = int(rng.integers(0, 13))
        a = random_sparse(rng,m,n)
        dense = a.toarray()
        rr: list[int] = []
        cc: list[int] = []
        src = sparse.csr_matrix((0,0))
        newrows, newcols = 0, 0
        if op == 1:
            take = int(rng.integers(0, 30)) if m*n else 0
            rr = rng.integers(m,size=take).tolist() if take else []
            cc = rng.integers(n,size=take).tolist() if take else []
            for r,c in zip(rr,cc):
                dense[r,c] = 0
        elif op == 2:
            take = int(rng.integers(0, 2*m+1)) if m else 0
            rr = rng.integers(m,size=take).tolist() if take else []
            dense = np.delete(dense,np.unique(np.asarray(rr,dtype=int)),axis=0)
        elif op == 3:
            take = int(rng.integers(0, 2*n+1)) if n else 0
            cc = rng.integers(n,size=take).tolist() if take else []
            dense = np.delete(dense,np.unique(np.asarray(cc,dtype=int)),axis=1)
        elif op == 4:
            count = m*n
            if count == 0:
                options = [(0,0),(0,3),(4,0),(0,8)]
                newrows,newcols = options[int(rng.integers(0,len(options)))]
            else:
                divisors = [p for p in range(1,count+1) if count%p==0]
                newrows = divisors[int(rng.integers(0,len(divisors)))]
                newcols = count//newrows
            dense = a.reshape((newrows,newcols),order='C').toarray()
        elif op == 5:
            take_r = int(rng.integers(0,15)) if m else 0
            take_c = int(rng.integers(0,15)) if n else 0
            rr = rng.integers(m,size=take_r).tolist() if take_r else []
            cc = rng.integers(n,size=take_c).tolist() if take_c else []
            src = random_sparse(rng,take_r,take_c)
            values = src.toarray()
            # Match existing SciFort column-major scatter order.
            for j,c in enumerate(cc):
                for i,r in enumerate(rr):
                    dense[r,c] = values[i,j]
        src_nz = src.nnz if op == 5 else 0
        buf.write(f'{op} {m} {n} {a.nnz} {len(rr)} {len(cc)} '
                  f'{newrows} {newcols} {src_nz}\n')
        serialize(buf,a)
        for r in rr:
            buf.write(f'{r+1}\n')
        for c in cc:
            buf.write(f'{c+1}\n')
        if op == 5:
            serialize(buf,src)
        expected.append((dense,op))
    return buf.getvalue(),expected


def run(mode: str, cases: int, seed: int):
    rng = np.random.default_rng(seed)
    input_text, expected = generate(rng,cases)
    with tempfile.TemporaryDirectory(prefix=f'scifort-delete-{mode}-') as temp:
        work = Path(temp)
        flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
        cmd = ['gfortran','-std=f2018',*flags,'-J',str(work),'-I',str(work)]
        for source,basename in [('src/kinds.f90','kinds'),('src/sparse.f90','sparse')]:
            subprocess.run([*cmd,'-c',str(ROOT/source),'-o',str(work/(basename+'.o'))],cwd=work,check=True)
        exe=work/'compare_deletion'
        subprocess.run([*cmd,str(ROOT/'tools/sparse_reference/sparse_deletion_reference_driver.f90'),
                        str(work/'kinds.o'),str(work/'sparse.o'),'-o',str(exe)],cwd=work,check=True)
        output = subprocess.run([str(exe)],input=input_text,text=True,stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,check=True,cwd=work).stdout.splitlines()
    cursor = 0
    checked = {op:0 for op in range(1,6)}
    for index,(reference,op) in enumerate(expected):
        head = output[cursor].split()
        cursor += 1
        if head[0] != 'RESULT' or len(head) != 5:
            raise AssertionError(f'bad header at case {index}: {head}')
        info,m,n,nnz = map(int,head[1:])
        if info != 0:
            raise AssertionError(f'failure status {info} at case {index}, op {op}')
        if (m,n) != reference.shape:
            raise AssertionError(f'shape mismatch at case {index}, op {op}: {(m,n)} vs {reference.shape}')
        got = np.zeros((m,n))
        for _ in range(nnz):
            r,c,value = output[cursor].split()
            cursor += 1
            got[int(r)-1,int(c)-1] = float(value)
        np.testing.assert_allclose(got,reference,rtol=0,atol=1e-12,
                                   err_msg=f'case={index}, op={op}')
        checked[op] += 1
    assert cursor == len(output), f'unexpected extra output: {output[cursor:cursor+3]}'
    print(f'{mode}: {cases} randomized SciPy/NumPy checks passed; by operation {checked}',flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    p.add_argument('--cases',type=int,default=250)
    p.add_argument('--seed',type=int,default=20261105)
    args=p.parse_args()
    for mode in ('checked','optimized') if args.mode=='both' else (args.mode,):
        run(mode,args.cases,args.seed)

if __name__ == '__main__':
    main()
