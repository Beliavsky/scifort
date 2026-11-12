#!/usr/bin/env python3
"""Compare SciFort's LSQR and ILU(0) to SciPy and an independent ILU(0)."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import sparse
from scipy.sparse import linalg
from scipy import linalg as scipy_linalg

ROOT = Path(__file__).resolve().parents[2]


def create_cases(seed: int, count: int):
    rng = np.random.default_rng(seed)
    stream = [str(count)]
    cases = []
    for case_id in range(count):
        kind = case_id % 4
        n = int(rng.integers(2, 22))
        m = n + int(rng.integers(1, 8)) if kind == 0 else (
            max(1, n-int(rng.integers(1, min(n, 7)+1))) if kind == 1 else n)
        raw = sparse.random(m,n,density=float(rng.uniform(.06, .28)),format='coo',
                            random_state=rng,
                            data_rvs=lambda size: rng.uniform(-0.45,.45,size))
        dense = raw.toarray()
        if kind == 2 or kind == 3:
            # Strict diagonal dominance: ILU(0) is pivot-stable.
            for i in range(n):
                if i > 0: dense[i,i-1] = -.25
                if i < n-1: dense[i,i+1] = -.2
                dense[i,i] = 1.5 + np.sum(np.abs(dense[i,:]))
        else:
            for i in range(min(m,n)):
                dense[i,i] += 1.0 + np.sum(np.abs(dense[i,:]))
        a = sparse.csr_matrix(dense)
        xtrue = rng.normal(size=n)
        b = np.asarray(a@xtrue)
        if kind == 0 and case_id % 8 == 0:
            b = b + 0.03*rng.normal(size=m)
        damping = .15 if case_id % 7 == 0 else 0.0
        x0 = rng.normal(size=n) if case_id % 6 == 0 else np.zeros(n)
        coo = a.tocoo()
        rows = coo.row.tolist()
        cols = coo.col.tolist()
        data = coo.data.tolist()
        if kind > 1:
            # SciFort canonical CSR retains these zero structural entries.
            rows.append(0)
            cols.append(0)
            data.append(0.0)
        cap = 8*max(m,n)
        stream.append(f'{m} {n} {len(rows)} {damping:.17g} {cap}')
        stream.extend(f'{i+1} {j+1} {v:.17g}' for i,j,v in zip(rows,cols,data))
        stream.extend(f'{t:.17g}' for t in b)
        stream.extend(f'{t:.17g}' for t in x0)
        cases.append((a,b,damping,x0,kind))
    return '\n'.join(stream)+'\n', cases


def ilu0_python(a):
    """Independent no-fill Crout row updates on SciPy CSR's fixed sparsity."""
    a = a.tocsr(copy=True)
    a.sum_duplicates()
    a.sort_indices()
    n = a.shape[0]
    x = a.toarray()
    pattern = np.zeros((n,n),dtype=bool)
    for i in range(n):
        pattern[i,a.indices[a.indptr[i]:a.indptr[i+1]]] = True
    for i in range(n):
        for j in range(i):
            if not pattern[i,j]: continue
            if abs(x[j,j]) <= np.finfo(float).tiny:
                raise ValueError('zero pivot')
            x[i,j] /= x[j,j]
            for k in range(j+1,n):
                if pattern[j,k] and pattern[i,k]:
                    x[i,k] -= x[i,j]*x[j,k]
        if not pattern[i,i] or abs(x[i,i]) <= np.finfo(float).tiny:
            raise ValueError('missing or zero pivot')
    return x,pattern


def parse_vector(line, label, n):
    tokens = line.split()
    assert tokens[0] == label, (label,tokens[:4])
    values = np.asarray([float(v) for v in tokens[1:]])
    assert values.shape == (n,), (label,n,values.shape)
    return values


def validate(output, cases):
    lines = iter(output.splitlines())
    counts = {'LSQR': 0,'ILU0': 0, 'BiCGSTAB-ILU0': 0, 'GMRES-ILU0': 0}
    for case_no,(a,b,damp,x0,kind) in enumerate(cases,1):
        n = a.shape[1]
        try:
            assert next(lines) == f'CASE {case_no}'
            record = next(lines).split()
            assert record[0] == 'LSQR',record
            status,iterations = map(int,record[1:3])
            r1,r2,gradient = map(float,record[3:6])
            x = parse_vector(next(lines),'X',n)
            variance = parse_vector(next(lines),'VAR',n)
            scipy_out = linalg.lsqr(a,b,damp=damp,atol=1.e-12,btol=1.e-12,
                                   conlim=1.e14,iter_lim=8*max(a.shape),x0=x0,
                                   calc_var=True)
            expected = scipy_out[0]
            assert 0 <= status <= 3 or status == 7, (status,iterations)
            assert np.all(np.isfinite(variance)) and np.all(variance >= 0)
            np.testing.assert_allclose(r1,np.linalg.norm(b-a@x),rtol=1.e-8,atol=1.e-10)
            np.testing.assert_allclose(r2,np.hypot(np.linalg.norm(b-a@x),damp*np.linalg.norm(x-x0)),
                                       rtol=1.e-8,atol=1.e-10)
            g = a.T @ (b-a@x)-damp**2*(x-x0)
            np.testing.assert_allclose(gradient,np.linalg.norm(g),rtol=1.e-7,atol=1.e-10)
            np.testing.assert_allclose(x,expected,rtol=1.e-6,atol=1.e-6)
            # Variance estimates are iteration-sensitive in nearly rank-deficient
            # Lanczos recurrences. Positivity is checked here; diagonal systems
            # have exact analytic variance tests in the native regression.
            counts['LSQR'] += 1
            if kind >= 2:
                line = next(lines).split()
                assert line[0] == 'ILU' and int(line[1]) == 0, line
                lu_data = parse_vector(next(lines),'DATA',a.nnz)
                reference,pattern = ilu0_python(a)
                expected_data = reference[np.repeat(np.arange(n),np.diff(a.indptr)),a.indices]
                np.testing.assert_allclose(lu_data,expected_data,rtol=1.e-11,atol=1.e-11)
                # Compare factor solve against independent triangular L and U.
                line = next(lines).split()
                assert line[0] == 'ISOLVE' and int(line[1]) == 0, line
                ix = parse_vector(next(lines),'IX',n)
                L = np.tril(reference,-1)+np.eye(n)
                U = np.triu(reference)
                expected_ix = scipy_linalg.solve_triangular(U,scipy_linalg.solve_triangular(
                    L,b,lower=True,unit_diagonal=True),lower=False)
                np.testing.assert_allclose(ix,expected_ix,rtol=1.e-11,atol=1.e-11)
                counts['ILU0'] += 1
                for label,name in [('B','BiCGSTAB-ILU0'),('G','GMRES-ILU0')]:
                    line = next(lines).split()
                    assert line[0] == label
                    info,steps = map(int,line[1:3])
                    residual = float(line[3])
                    xx = parse_vector(next(lines),label+'X',n)
                    assert info == 0 and steps >= 0, (label,info,steps)
                    actual = np.linalg.norm(b-a@xx)
                    np.testing.assert_allclose(actual,residual,rtol=1.e-8,atol=1.e-10)
                    assert actual <= 1.e-11*np.linalg.norm(b)+1.e-10,(label,actual)
                    np.testing.assert_allclose(xx,np.linalg.solve(a.toarray(),b),rtol=1.e-8,atol=1.e-8)
                    counts[name] += 1
        except (AssertionError, ValueError, StopIteration) as exc:
            raise AssertionError(f'case {case_no}, kind={kind}, shape={a.shape}: {exc}') from exc
    assert next(lines,None) is None, 'trailing Fortran output'
    return counts


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cases',type=int,default=160)
    p.add_argument('--seed',type=int,default=20261104)
    p.add_argument('--mode',choices=('both','checked','optimized'),default='both')
    args = p.parse_args()
    stdin, cases = create_cases(args.seed,args.cases)
    modes = ('checked','optimized') if args.mode=='both' else (args.mode,)
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-lsqr-ilu-') as directory:
            root=Path(directory)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc=['gfortran','-std=f2018',*flags,'-J',str(root),'-I',str(root)]
            objects=[]
            for name,source in [('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                                ('lsmr','src/optimize/lsmr.f90'),('sparse_linalg','src/sparse_linalg.f90')]:
                path=root/f'{name}.o'
                subprocess.run([*cc,'-c',str(ROOT/source),'-o',str(path)],cwd=root,check=True)
                objects.append(str(path))
            exe=root/'reference_driver'
            subprocess.run([*cc,str(ROOT/'tools/sparse_linalg_reference/sparse_lsqr_ilu0_reference_driver.f90'),
                            *objects,'-o',str(exe)],cwd=root,check=True)
            result=subprocess.run([str(exe)],input=stdin,text=True,capture_output=True,cwd=root)
            if result.returncode:
                raise RuntimeError(f'{mode}: {result.stderr[-4000:]}\n{result.stdout[-4000:]}')
            counts=validate(result.stdout,cases)
            print(f'{mode}: {args.cases} graphs / matrices: {counts} PASS',flush=True)


if __name__=='__main__':
    main()
