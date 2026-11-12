#!/usr/bin/env python3
"""Randomized SciFort sparse direct-solver comparisons with SciPy and NumPy.

Includes independent LU reconstruction (P*A=L*U), direct solves and
triangular solves.  No external solver source code is used by SciFort.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sla

ROOT = Path(__file__).resolve().parents[2]


def make_cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    text = [str(count)]
    cases = []
    for case_id in range(count):
        kind = case_id % 6
        n = int(rng.integers(1, 27))
        nrhs = int(rng.integers(1, 5))
        density = float(rng.uniform(0.1, 0.6))
        a = rng.normal(scale=0.5, size=(n, n))
        a *= rng.random((n,n)) < density
        if kind < 2:
            for i in range(n):
                a[i,i] += np.sum(np.abs(a[i,:])) + .75
            if kind == 1 and n > 2:
                # Force row pivoting without losing nonsingularity.
                a[[0,1],:] = a[[1,0],:]
        else:
            a = np.tril(a) if kind in (2,4) else np.triu(a)
            if kind in (2,3):
                np.fill_diagonal(a,rng.choice([-1.,1.],n)*rng.uniform(.7,2.5,size=n))
            else:
                np.fill_diagonal(a,0.0)
        a0 = a.copy()
        original = sparse.coo_matrix(a)
        rows = original.row.tolist()
        columns = original.col.tolist()
        values = original.data.tolist()
        # Repeated coordinate on a structural zero, canceled in COO summation.
        rows.extend([0,0])
        columns.extend([n-1,n-1])
        values.extend([.25,-.25])
        # Explicit zero with repeated coordinate of a nonzero entry.
        if n > 1:
            rows.append(0)
            columns.append(0)
            values.append(0.0)
        exact = rng.normal(size=(n,nrhs))
        effective = a0.copy()
        if kind >= 4:
            np.fill_diagonal(effective,1.0)
        b = effective @ exact
        text.append(f'{n} {len(rows)} {nrhs} {kind} {case_id % 2}')
        text.extend(f'{r+1} {c+1} {v:.17g}' for r,c,v in zip(rows,columns,values))
        text.extend(f'{v:.17g}' for col in b.T for v in col)
        cases.append((kind,a0,b,exact))
    return '\n'.join(text)+'\n',cases


def read_values(iterator, label, n):
    line=next(iterator).split()
    assert line and line[0] == label, (label,line[:5])
    result=np.asarray([float(v) for v in line[1:]])
    assert result.shape == (n,),(label,result.shape,n)
    return result


def validate(output: str, cases):
    lines=iter(output.splitlines())
    count_lu=count_tri=pivoted=0
    for case_no,(kind,a,b,exact) in enumerate(cases,1):
        n,nrhs=b.shape
        try:
            assert next(lines) == f'CASE {case_no}'
            if kind < 2:
                name,status=next(lines).split()
                assert name == 'LU' and int(status) == 0,(name,status)
                perm=read_values(lines,'PERM',n).astype(int)-1
                assert set(perm) == set(range(n))
                if np.any(perm != np.arange(n)):
                    pivoted+=1
                word,nnz=next(lines).split()
                assert word == 'NNZ'
                nnz=int(nnz)
                rr=[];cc=[];vv=[]
                for _ in range(nnz):
                    r,c,v=next(lines).split()
                    rr.append(int(r)-1)
                    cc.append(int(c)-1)
                    vv.append(float(v))
                combined=sparse.coo_matrix((vv,(rr,cc)),shape=(n,n)).toarray()
                lower=np.tril(combined,-1)+np.eye(n)
                upper=np.triu(combined)
                np.testing.assert_allclose(lower@upper,a[perm,:],rtol=1.e-11,atol=5.e-11)
                name,status=next(lines).split()
                assert name=='FSOLVE' and int(status)==0
                x=np.column_stack([read_values(lines,'X',n) for _ in range(nrhs)])
                name,status=next(lines).split()
                assert name=='SOLVE' and int(status)==0
                y=np.column_stack([read_values(lines,'Y',n) for _ in range(nrhs)])
                assert np.max(np.abs(a@x-b)) < 4.e-10
                np.testing.assert_allclose(x,exact,rtol=2.e-10,atol=2.e-10)
                np.testing.assert_allclose(y,x,rtol=2.e-12,atol=2.e-12)
                # SciPy reference uses SuperLU, not our row elimination algorithm.
                reference=sla.spsolve(sparse.csc_matrix(a),b[:,0])
                np.testing.assert_allclose(x[:,0],reference,rtol=1.e-10,atol=1.e-10)
                count_lu+=1
            else:
                name,status=next(lines).split()
                assert name == 'TRI' and int(status) == 0,(name,status)
                x=np.column_stack([read_values(lines,'X',n) for _ in range(nrhs)])
                lower=kind in (2,4)
                unit=kind >= 4
                scipy_a=sparse.csr_matrix(a)
                for j in range(nrhs):
                    reference=sla.spsolve_triangular(scipy_a,b[:,j],lower=lower,unit_diagonal=unit)
                    np.testing.assert_allclose(x[:,j],reference,rtol=1.e-11,atol=1.e-11)
                np.testing.assert_allclose(x,exact,rtol=1.e-11,atol=1.e-11)
                count_tri+=1
        except (AssertionError,ValueError,StopIteration) as exc:
            raise AssertionError(f'case {case_no}, shape {a.shape}, kind {kind}: {exc}') from exc
    assert next(lines,None) is None,'unexpected output after final case'
    return count_lu,count_tri,pivoted


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=180)
    parser.add_argument('--seed',type=int,default=20261105)
    args=parser.parse_args()
    data,cases=make_cases(args.cases,args.seed)
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-direct-') as tmp:
            build=Path(tmp)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc=['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
            objects=[]
            for name,source in (('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                                ('lsmr','src/optimize/lsmr.f90'),('sparse_linalg','src/sparse_linalg.f90')):
                obj=build/f'{name}.o'
                subprocess.run([*cc,'-c',str(ROOT/source),'-o',str(obj)],cwd=build,check=True)
                objects.append(str(obj))
            exe=build/'reference_driver'
            subprocess.run([*cc,str(ROOT/'tools/sparse_linalg_reference/sparse_direct_reference_driver.f90'),
                            *objects,'-o',str(exe)],cwd=build,check=True)
            result=subprocess.run([str(exe)],input=data,text=True,capture_output=True,cwd=build)
            if result.returncode:
                raise RuntimeError(f'{mode}: {result.stderr[-3000:]}\n{result.stdout[-3000:]}')
            lu,tri,pivoted=validate(result.stdout,cases)
            print(f'{mode}: {args.cases} cases; {lu} LU, {tri} triangular; {pivoted} row-permuted: PASS',flush=True)


if __name__ == '__main__':
    main()
