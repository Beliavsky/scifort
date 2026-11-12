#!/usr/bin/env python3
"""Compare SciFort symmetric Ritz eigenpairs to SciPy and dense NumPy eigensolvers."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]


def create_cases(seed, count):
    rng = np.random.default_rng(seed)
    stream = [str(count)]
    cases = []
    for ci in range(count):
        n = int(rng.integers(5, 37))
        k = int(rng.integers(1,min(n-1,3)+1))
        which = ('LM','LA','SA','SM')[ci % 4]
        diag = np.linspace(-6.0,5.0,n)
        if which == 'LM':
            diag = np.linspace(-8.5,5.0,n) if (ci//4)%2 else np.linspace(1.0,9.0,n)
        diag += rng.uniform(-.03,.03,size=n)
        a=np.diag(diag)
        # A tridiagonal perturbation gives eigenvectors with nontrivial
        # support while retaining well-separated symmetric Ritz eigenvalues.
        for i in range(n-1):
            a[i,i+1]=a[i+1,i]=float(rng.uniform(-.12,.12))
        # Some random symmetric off-diagonal structures test the small
        # Rayleigh--Ritz eigensolver beyond narrowly banded matrices.
        if ci % 7 == 0:
            upper = np.triu((rng.random((n, n)) < .25) * rng.uniform(-.2, .2, (n, n)), 1)
            a += upper + upper.T
        # Include reduced block subspaces, not just complete bases.
        nc = n if which == 'SM' or ci % 3 == 0 or ci % 7 == 0 else min(n, max(2*k+8,16))
        cap = 650 if nc < n else 10
        tol = 2.e-7 if nc < n else 1.e-10
        coo=sparse.coo_matrix(a)
        stream.append(f'{n} {coo.nnz} {k} {which} {tol:.17g} {cap} {nc}')
        stream.extend(f'{i+1} {j+1} {v:.17g}' for i,j,v in zip(coo.row,coo.col,coo.data))
        cases.append((a,k,which,tol,nc))
    return '\n'.join(stream)+'\n',cases


def expected_eigs(a,k,which):
    spectrum=np.linalg.eigvalsh(a)
    if which=='LA':
        vals=spectrum[-k:]
    elif which=='SA':
        vals=spectrum[:k]
    elif which=='LM':
        vals=spectrum[np.argsort(np.abs(spectrum))[-k:]]
    else:
        vals=spectrum[np.argsort(np.abs(spectrum))[:k]]
    return np.sort(vals)


def verify(output,cases):
    lines=iter(output.splitlines())
    counts={'cases':0,'reduced_subspaces':0,'scipy_near_ties':0,'LM':0,'LA':0,'SA':0,'SM':0}
    for ci,(a,k,which,tol,nc) in enumerate(cases,1):
        try:
            assert next(lines).strip()==f'CASE {ci}'
            record=next(lines).split()
            assert record[0]=='EIGSH'
            status,steps=map(int,record[1:3])
            record=next(lines).split()
            assert record[0]=='LAMBDA'
            values=np.array([float(v) for v in record[1:]])
            record=next(lines).split()
            assert record[0]=='RESIDUAL'
            residuals=np.array([float(v) for v in record[1:]])
            vecs=[]
            for j in range(k):
                record=next(lines).split()
                assert record[:2]==['VECTOR',str(j+1)]
                vecs.append([float(v) for v in record[2:]])
            v=np.array(vecs).T
            assert status==0, f'nonconvergence status={status} steps={steps} which={which}, ncv={nc}, n={len(a)}'
            assert values.shape==(k,) and v.shape==(len(a),k)
            reference=expected_eigs(a,k,which)
            np.testing.assert_allclose(values,reference,rtol=3e-5,atol=3e-5)
            # SciPy ARPACK agrees with the fully assembled self-adjoint matrix.
            scipy_vals=np.sort(eigsh(sparse.csr_matrix(a),k,which=which,return_eigenvectors=False,
                                    tol=1e-10,maxiter=1000))
            spectrum=np.linalg.eigvalsh(a)
            if which in ('LM','SM'):
                sorted_by_target=np.sort(np.abs(spectrum))
                gap=(sorted_by_target[k]-sorted_by_target[k-1]) if which=='SM' else (sorted_by_target[-k]-sorted_by_target[-k-1])
            else:
                gap=(spectrum[-k]-spectrum[-k-1]) if which=='LA' else (spectrum[k]-spectrum[k-1])
            if gap >= .01:
                np.testing.assert_allclose(values,scipy_vals,rtol=3e-5,atol=3e-5)
            else:
                # ARPACK can return the wrong member of a nearly tied
                # |lambda| pair; the independently computed dense spectrum
                # is the stricter oracle for exact rank selection.
                counts['scipy_near_ties']+=1
            actual=np.linalg.norm(a@v-v*values[None,:],axis=0)
            np.testing.assert_allclose(actual,residuals,rtol=1e-6,atol=1e-10)
            assert np.max(actual)<=tol*max(1.0,np.linalg.norm(a,np.inf)) * 1.01
            np.testing.assert_allclose(v.T@v,np.eye(k),rtol=1e-6,atol=1e-6)
            counts['cases']+=1
            counts['reduced_subspaces']+=int(nc<len(a))
            counts[which]+=1
        except (AssertionError,ValueError,StopIteration) as e:
            raise AssertionError(f'case {ci}, mode {which}, size {len(a)}, k={k}: {e}') from e
    assert next(lines,None) is None
    return counts


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--cases',type=int,default=100)
    p.add_argument('--seed',type=int,default=20261106)
    p.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    options=p.parse_args()
    data,cases=create_cases(options.seed,options.cases)
    for mode in ('checked','optimized') if options.mode=='both' else (options.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-eigsh-') as d:
            tmp=Path(d)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc=['gfortran','-std=f2018',*flags,'-J',str(tmp),'-I',str(tmp)]
            obj=[]
            for name,source in [('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                                ('lsmr','src/optimize/lsmr.f90'),('sparse_linalg','src/sparse_linalg.f90')]:
                output=tmp/(name+'.o')
                subprocess.run([*cc,'-c',str(ROOT/source),'-o',str(output)],cwd=tmp,check=True)
                obj.append(str(output))
            exe=tmp/'reference'
            subprocess.run([*cc,str(ROOT/'tools/sparse_linalg_reference/sparse_eigsh_reference_driver.f90'),
                            *obj,'-o',str(exe)],cwd=tmp,check=True)
            result=subprocess.run([str(exe)],input=data,text=True,capture_output=True,check=True,cwd=tmp)
            print(mode,verify(result.stdout,cases),flush=True)

if __name__=='__main__':
    main()
