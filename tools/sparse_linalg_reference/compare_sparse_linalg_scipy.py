#!/usr/bin/env python3
"""Randomized SciPy numerical comparisons for SciFort sparse iterative solvers."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse
from scipy.sparse import linalg

ROOT = Path(__file__).resolve().parents[2]


def cases_from_seed(seed: int, count: int):
    rng = np.random.default_rng(seed)
    lines = [str(count)]
    cases = []
    for case_no in range(count):
        kind = case_no % 3 + 1
        n = int(rng.integers(2, 18))
        m = n if kind != 3 else n + int(rng.integers(1, 10))
        # Build reproducible, diagonally well-conditioned random matrices.
        offdiag = sparse.random(m, n, density=float(rng.uniform(.08,.4)),
                                data_rvs=lambda count: rng.uniform(-.6,.6,count),
                                random_state=rng, format='csr')
        if kind == 1:
            offdiag = (offdiag + offdiag.T)*.5
            dense = offdiag.toarray()
            dense[np.diag_indices(n)] += 2.0 + np.sum(np.abs(dense), axis=1)
        else:
            dense = offdiag.toarray()
            dense[np.arange(n),np.arange(n)] += 2.0 + np.sum(np.abs(dense[:n]),axis=1)
        a = sparse.csr_matrix(dense)
        xtrue = rng.normal(size=n)
        b = a @ xtrue
        # Nonzero residual least squares cases stress stopping logic.
        if kind == 3 and case_no % 2 == 0:
            b = b + 0.04*rng.normal(size=m)
        # Store some explicit zeros and duplicate triplets to exercise COO canonicalization.
        coo = a.tocoo()
        rows = coo.row.tolist()
        cols = coo.col.tolist()
        values = coo.data.tolist()
        if n > 0:
            rows.append(0)
            cols.append(0)
            values.append(0.0)
        lines.append(f'{kind} {m} {n} {len(values)}')
        lines.extend(f'{i+1} {j+1} {v:.17g}' for i,j,v in zip(rows,cols,values))
        lines.extend(f'{v:.17g}' for v in b)
        cases.append((kind, a, b))
    return '\n'.join(lines)+'\n', cases


def reference(label, a, b):
    opts = dict(rtol=1.e-10, atol=0.0)
    m = a.shape[1]
    if label in ('C','CJ'):
        M = sparse.diags(1.0/a.diagonal(),format='csr') if label.endswith('J') else None
        return linalg.cg(a,b,**opts,maxiter=10*m,M=M)
    if label in ('B','BJ'):
        M = sparse.diags(1.0/a.diagonal(),format='csr') if label.endswith('J') else None
        return linalg.bicgstab(a,b,**opts,maxiter=10*m,M=M)
    if label in ('G','GJ'):
        M = sparse.diags(1.0/a.diagonal(),format='csr') if label.endswith('J') else None
        return linalg.gmres(a,b,**opts,restart=min(8,m),maxiter=100,M=M)
    kw = dict(atol=1.e-12,btol=1.e-12,maxiter=10*max(a.shape))
    if label == 'LD':
        kw['damp'] = 0.15
    result = linalg.lsmr(a,b,**kw)
    return result[0], result[1]


def verify(output: str, cases) -> dict[str,int]:
    lines = iter(output.splitlines())
    counts: dict[str,int] = {}
    for case_id,(kind,a,b) in enumerate(cases,1):
        line = next(lines,None)
        assert line == f'CASE {case_id}', (case_id,line)
        methods = (('C','CJ') if kind == 1 else ()) + (
                  ('B','BJ','G','GJ') if kind != 3 else ()) + (
                  ('L','LD') if kind == 3 else ('L',))
        for label in methods:
            tokens = next(lines).split()
            assert tokens[0] == label, (case_id,label,tokens)
            status = int(tokens[1])
            steps = int(tokens[2])
            reported_residual = float(tokens[3])
            x = np.array([float(s) for s in tokens[4:]])
            assert x.shape == (a.shape[1],), (case_id,label,x.shape)
            sx, si = reference(label,a,b)
            try:
                if label not in ('L','LD'):
                    assert status == 0 and si == 0, (status,si)
                    assert 0 <= steps <= 10*a.shape[1]*100
                    actual = float(np.linalg.norm(b-a@x))
                    assert actual <= 1.e-10*np.linalg.norm(b) + 1.e-12, actual
                else:
                    assert status >= 0, status
                    actual = float(np.linalg.norm(b-a@x))
                np.testing.assert_allclose(reported_residual,actual,atol=1.e-10,rtol=1.e-9)
                np.testing.assert_allclose(x,sx,rtol=1.e-7,atol=1.e-7)
            except (AssertionError, ValueError) as exc:
                raise AssertionError(f'case {case_id}, method {label}: {exc}') from exc
            counts[label] = counts.get(label,0) + 1
    assert next(lines,None) is None, 'extra Fortran output'
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('both','checked','optimized'),default='both')
    parser.add_argument('--cases',type=int,default=240)
    parser.add_argument('--seed',type=int,default=20261104)
    opts = parser.parse_args()
    stdin,cases = cases_from_seed(opts.seed,opts.cases)
    for mode in (('checked','optimized') if opts.mode == 'both' else (opts.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-linalg-reference-') as d:
            build=Path(d)
            flags=('-O0','-g','-fcheck=all','-fbacktrace') if mode=='checked' else ('-O3',)
            cc=['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
            objects=[]
            for name,file in [('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                              ('lsmr','src/optimize/lsmr.f90'),
                              ('sparse_linalg','src/sparse_linalg.f90')]:
                output=build/f'{name}.o'
                subprocess.run([*cc,'-c',str(ROOT/file),'-o',str(output)],cwd=build,check=True)
                objects.append(str(output))
            exe=build/'sparse_linalg_reference'
            subprocess.run([*cc,str(ROOT/'tools/sparse_linalg_reference/sparse_linalg_reference_driver.f90'),
                            *objects,'-o',str(exe)],cwd=build,check=True)
            result=subprocess.run([str(exe)],input=stdin,text=True,capture_output=True,cwd=build)
            if result.returncode:
                raise RuntimeError(f'{mode} Fortran failed: {result.stderr[-4000:]}\n{result.stdout[-4000:]}')
            counts=verify(result.stdout,cases)
            print(f'{mode}: {opts.cases} randomized matrices, {sum(counts.values())} '
                  f'SciPy solver comparisons PASS ({counts})',flush=True)

if __name__=='__main__':
    main()
