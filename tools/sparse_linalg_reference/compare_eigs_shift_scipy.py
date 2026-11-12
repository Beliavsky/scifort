#!/usr/bin/env python3
"""Independent numerical checks of complex/real shift-invert CSR eigs."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigs as scipy_eigs, ArpackNoConvergence

ROOT = Path(__file__).resolve().parents[2]


def build_cases(seed: int, count: int, smallest: bool = False):
    rng = np.random.default_rng(seed)
    data = [str(count)]
    cases = []
    for case in range(count):
        n = int(rng.integers(14, 33))
        a = np.diag(np.linspace(-4, 4, n))
        for i in range(n-1):
            a[i,i+1] = rng.uniform(-0.012,0.012)
            a[i+1,i] = rng.uniform(-0.012,0.012)
        if smallest:
            a = np.diag(np.linspace(2.0,9.0,n))
            a[0,0] = 0.13
            a[1,1] = -0.27
            for i in range(n-1):
                a[i,i+1] = rng.uniform(-0.01,0.01)
                a[i+1,i] = rng.uniform(-0.01,0.01)
            sigma = 0.0j
            k = 2
        elif case % 3 == 0:
            # Real shift between adjacent simple eigenvalues.
            middle = int(rng.integers(3,n-3))
            sigma = complex(float(a[middle,middle] +
                            0.22*(a[middle+1,middle+1]-a[middle,middle])))
            k = 2
        else:
            # Isolated complex conjugate eigenvalues of a real operator.
            re = float(rng.uniform(-1.5,1.5))
            im = float(rng.uniform(2,5))
            a[0:2,0:2] = [[re,-im],[im,re]]
            k = 1
            sigma = complex(re+.11, (1 if case%3 == 1 else -1)*(im-.13))
        tol = 5e-8
        cap = 110
        ncv = min(n,12)
        coo = sparse.coo_matrix(a)
        data.append(f'{n} {coo.nnz} {k} {sigma.real:.17g} {sigma.imag:.17g} {tol:.17g} {cap} {ncv} {'SM' if smallest else 'LM'}')
        data.extend(f'{r+1} {c+1} {v:.17g}' for r,c,v in zip(coo.row,coo.col,coo.data))
        cases.append((a,k,sigma,tol,smallest))
    return '\n'.join(data)+'\n', cases


def verify(output, cases):
    stream = iter(output.splitlines())
    totals = {'cases': 0, 'complex_shifts': 0, 'largest_residual': 0.0, 'max_iter': 0,
              'scipy_verified': 0}
    for case, (a,k,sigma,tol,smallest) in enumerate(cases,1):
        assert next(stream).strip() == f'CASE {case}'
        status, steps = map(int,next(stream).split()[1:])
        line=next(stream).split()
        assert line[0]=='LAMBDA'
        numbers=np.array(list(map(float,line[1:])))
        vals=numbers[::2]+1j*numbers[1::2]
        line=next(stream).split()
        assert line[0]=='RESIDUAL'
        errs=np.array(list(map(float,line[1:])))
        vectors=[]
        for j in range(k):
            line=next(stream).split()
            assert line[:2] == ['VECTOR',str(j+1)]
            numbers=np.array(list(map(float,line[2:])))
            vectors.append(numbers[::2]+1j*numbers[1::2])
        vectors=np.array(vectors).T
        assert status==0, f'case {case}: status={status}, steps={steps}, sigma={sigma}, n={len(a)}'
        spectrum=np.linalg.eigvals(a)
        expected=sorted(spectrum,key=lambda lam:abs(lam-sigma))[:k]
        for lam in vals:
            assert min(abs(lam-expected)) < 1e-5, f'case {case}: wrong value {lam}, {expected}'
        actual=np.linalg.norm(a@vectors-vectors*vals[None,:],axis=0)
        np.testing.assert_allclose(actual,errs,rtol=1e-4,atol=1e-10)
        assert max(actual)<tol*max(1,np.max(np.sum(abs(a),axis=1)))*1.1, (case,actual)
        np.testing.assert_allclose(np.linalg.norm(vectors,axis=0),1,atol=1e-7)
        # SciPy on complex CSR uses the full shifted operator.  Its real
        # CSR API defaults to an OPpart projection for nonreal sigma.
        if k<len(a)-2:
            try:
                reference=scipy_eigs(sparse.csr_matrix(a,dtype=complex),k=k,sigma=sigma,which='LM',
                                     tol=1e-9,maxiter=1000,return_eigenvectors=False)
                for lam in reference:
                    assert min(abs(lam-spectrum))<1e-5
                for lam in vals:
                    assert min(abs(lam-reference))<1e-5, (case,lam,reference)
                totals['scipy_verified'] += 1
            except ArpackNoConvergence:
                pass
        totals['cases']+=1
        totals['complex_shifts']+=int(sigma.imag!=0)
        totals['largest_residual']=max(totals['largest_residual'],float(max(actual)))
        totals['max_iter']=max(totals['max_iter'],steps)
    assert next(stream,None) is None
    return totals


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--cases',type=int,default=60)
    parser.add_argument('--seed',type=int,default=20261104)
    parser.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    parser.add_argument('--smallest',action='store_true',help='Test reduced-subspace SM shift-invert at zero')
    options=parser.parse_args()
    data,cases=build_cases(options.seed,options.cases,options.smallest)
    for mode in ('checked','optimized') if options.mode=='both' else (options.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-eigs-shift-') as directory:
            tmp=Path(directory)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc=['gfortran','-std=f2018',*flags,'-J',str(tmp),'-I',str(tmp)]
            objects=[]
            for name, source in [('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                                 ('lsmr','src/optimize/lsmr.f90'),('sparse_linalg','src/sparse_linalg.f90')]:
                obj=tmp/(name+'.o')
                subprocess.run([*cc,'-c',str(ROOT/source),'-o',str(obj)],cwd=tmp,check=True)
                objects.append(str(obj))
            exe=tmp/'ref'
            subprocess.run([*cc,str(ROOT/'tools/sparse_linalg_reference/sparse_eigs_shift_reference_driver.f90'),
                            *objects,'-o',str(exe)],cwd=tmp,check=True)
            run=subprocess.run([str(exe)],input=data,text=True,capture_output=True,cwd=tmp,check=True)
            print(mode,verify(run.stdout,cases),flush=True)


if __name__=='__main__':
    main()
