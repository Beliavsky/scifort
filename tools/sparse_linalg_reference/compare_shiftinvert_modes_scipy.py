#!/usr/bin/env python3
"""Validate normal shifted eigsh LA, SA, SM transformed modes against SciPy."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, ArpackNoConvergence

ROOT = Path(__file__).resolve().parents[2]
MODES=('LA','SA','SM')


def cases_for(seed, count):
    rng=np.random.default_rng(seed)
    stream=[str(count)]
    samples=[]
    for t in range(count):
        mode=MODES[t%3]
        n=int(rng.integers(14, 43))
        k=int(rng.integers(1, 4))
        diag=np.linspace(-5.0,6.0,n)
        a=np.diag(diag)
        for i in range(n-1):
            val=float(rng.uniform(-.015,.015))
            a[i,i+1]=a[i+1,i]=val
        pivot=int(rng.integers(4,n-5))
        sigma=float(diag[pivot]+.27*(diag[pivot+1]-diag[pivot]))
        nc=n if mode=='SM' or t%3==0 else min(n,16)
        cap=140 if nc<n else 3
        tol=4e-7
        coo=sparse.coo_matrix(a)
        stream.append(f'{n} {coo.nnz} {k} {sigma:.17g} {tol:.17g} {cap} {nc} {mode}')
        stream.extend(f'{i+1} {j+1} {v:.17g}' for i,j,v in zip(coo.row,coo.col,coo.data))
        samples.append((a,k,sigma,mode,tol,nc))
    return '\n'.join(stream)+'\n',samples


def verify(output,samples):
    lines=iter(output.splitlines())
    counts={'cases':0,'reduced_subspaces':0,'max_residual':0.0,
            'scipy_converged':0,'scipy_nonconverged':0}
    for ic,(a,k,sigma,mode,tol,nc) in enumerate(samples,1):
        assert next(lines).strip()==f'CASE {ic}'
        status,its=map(int,next(lines).split()[1:])
        vals=np.array([float(x) for x in next(lines).split()[1:]])
        errors=np.array([float(x) for x in next(lines).split()[1:]])
        vecs=[]
        for j in range(k):
            parts=next(lines).split()
            assert parts[:2]==['VECTOR',str(j+1)]
            vecs.append([float(x) for x in parts[2:]])
        vecs=np.array(vecs).T
        assert status==0, f'case {ic} {mode} n={len(a)} ncv={nc}: status {status}, iterations {its}'
        spectrum=np.linalg.eigvalsh(a)
        transform=1/(spectrum-sigma)
        if mode=='LA':
            chosen=np.argsort(transform)[-k:]
        elif mode=='SA':
            chosen=np.argsort(transform)[:k]
        else:
            chosen=np.argsort(abs(transform))[:k]
        expected=np.sort(spectrum[chosen])
        np.testing.assert_allclose(vals,expected,atol=1e-4,rtol=1e-4)
        try:
            reference=np.sort(eigsh(sparse.csr_matrix(a),k=k,sigma=sigma,which=mode,
                                    tol=1e-11,maxiter=2000,return_eigenvectors=False))
            np.testing.assert_allclose(vals,reference,atol=1e-4,rtol=1e-4)
            counts['scipy_converged']+=1
        except ArpackNoConvergence:
            # ARPACK can struggle with which='SM' even with shift-invert.
            # The independent dense eigvalsh reference above still checks it.
            counts['scipy_nonconverged']+=1
        actual=np.linalg.norm(a@vecs-vecs*vals[None,:],axis=0)
        np.testing.assert_allclose(actual,errors,atol=1e-10,rtol=1e-4)
        assert max(actual)<=tol*max(1,max(abs(a).sum(axis=1)))*1.05
        counts['cases']+=1
        counts['reduced_subspaces']+=int(nc<len(a))
        counts['max_residual']=max(counts['max_residual'],float(max(actual)))
    assert next(lines,None) is None
    return counts


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--cases',type=int,default=75)
    p.add_argument('--seed',type=int,default=20261104)
    p.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    args=p.parse_args()
    data,samples=cases_for(args.seed,args.cases)
    for mode in ('checked','optimized') if args.mode=='both' else (args.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-eigsh-shift-modes-') as d:
            tmp=Path(d)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc=['gfortran','-std=f2018',*flags,'-J',str(tmp),'-I',str(tmp)]
            objects=[]
            for name,source in [('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                                ('lsmr','src/optimize/lsmr.f90'),('sparse_linalg','src/sparse_linalg.f90')]:
                obj=tmp/f'{name}.o'
                subprocess.run([*cc,'-c',str(ROOT/source),'-o',str(obj)],check=True,cwd=tmp)
                objects.append(str(obj))
            exe=tmp/'reference'
            subprocess.run([*cc,str(ROOT/'tools/sparse_linalg_reference/sparse_shiftinvert_reference_driver.f90'),
                            *objects,'-o',str(exe)],check=True,cwd=tmp)
            run=subprocess.run([str(exe)],input=data,text=True,capture_output=True,check=True,cwd=tmp)
            print(mode,verify(run.stdout,samples),flush=True)


if __name__=='__main__':
    main()
