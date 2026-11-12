#!/usr/bin/env python3
"""Cross-check symmetric interior shift-invert eigenpairs against SciPy."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

ROOT=Path(__file__).resolve().parents[2]

def generate(seed,count):
    rng=np.random.default_rng(seed)
    stream=[str(count)]
    cases=[]
    for case in range(count):
        n=int(rng.integers(8,45))
        k=int(rng.integers(1,4))
        diag=np.linspace(-6,7,n)+rng.uniform(-.08,.08,n)
        a=np.diag(diag)
        for i in range(n-1):
            a[i,i+1]=a[i+1,i]=float(rng.uniform(-.035,.035))
        if case%7==0:
            a[0,2]=a[2,0]=0.01
        # Shift between neighboring eigenvalues, avoiding near-singular solves.
        position=int(rng.integers(1,n-2))
        sigma=float(diag[position]+.26*(diag[position+1]-diag[position]))
        nc=n if case%3 else min(n,18)
        cap=90 if nc<n else 3
        tol=2e-7
        coo=sparse.coo_matrix(a)
        stream.append(f'{n} {coo.nnz} {k} {sigma:.17g} {tol:.17g} {cap} {nc} LM')
        stream.extend(f'{i+1} {j+1} {v:.17g}' for i,j,v in zip(coo.row,coo.col,coo.data))
        cases.append((a,k,sigma,tol,nc))
    return '\n'.join(stream)+'\n',cases

def verify(text,cases):
    lines=iter(text.splitlines())
    summary={'cases':0,'reduced_subspaces':0,'max_residual':0.0}
    for ic,(a,k,sigma,tol,nc) in enumerate(cases,1):
        assert next(lines).strip()==f'CASE {ic}'
        status,its=map(int,next(lines).split()[1:])
        values=np.array([float(v) for v in next(lines).split()[1:]])
        residuals=np.array([float(v) for v in next(lines).split()[1:]])
        vecs=[]
        for j in range(k):
            parts=next(lines).split()
            assert parts[:2]==['VECTOR',str(j+1)]
            vecs.append([float(v) for v in parts[2:]])
        vecs=np.array(vecs).T
        assert status==0, f'case {ic}: code {status}, n={len(a)} sigma={sigma}, k={k} nc={nc}'
        dense=np.linalg.eigvalsh(a)
        exact=np.sort(dense[np.argsort(abs(dense-sigma))[:k]])
        np.testing.assert_allclose(values,exact,atol=2e-5,rtol=2e-5)
        scipy_vals=np.sort(eigsh(sparse.csr_matrix(a),k=k,which='LM',sigma=sigma,
                                tol=1e-12,return_eigenvectors=False))
        np.testing.assert_allclose(values,scipy_vals,atol=2e-5,rtol=2e-5)
        measured=np.linalg.norm(a@vecs-vecs*values[None,:],axis=0)
        np.testing.assert_allclose(measured,residuals,atol=1e-10,rtol=1e-4)
        assert max(measured) <= tol*max(1,np.max(np.sum(abs(a),axis=1)))*1.01
        summary['max_residual']=max(summary['max_residual'],float(max(measured)))
        summary['cases']+=1
        summary['reduced_subspaces']+=int(nc<len(a))
    assert next(lines,None) is None
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--cases',type=int,default=90)
    p.add_argument('--seed',type=int,default=20261104)
    p.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    options=p.parse_args()
    data,cases=generate(options.seed,options.cases)
    for mode in ('checked','optimized') if options.mode=='both' else (options.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-shiftinvert-') as d:
            tmp=Path(d)
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            cc=['gfortran','-std=f2018',*flags,'-J',str(tmp),'-I',str(tmp)]
            objects=[]
            for name,file in [('kinds','src/kinds.f90'),('sparse','src/sparse.f90'),
                              ('lsmr','src/optimize/lsmr.f90'),('sparse_linalg','src/sparse_linalg.f90')]:
                obj=tmp/(name+'.o')
                subprocess.run([*cc,'-c',str(ROOT/file),'-o',str(obj)],check=True,cwd=tmp)
                objects.append(str(obj))
            exe=tmp/'reference'
            subprocess.run([*cc,str(ROOT/'tools/sparse_linalg_reference/sparse_shiftinvert_reference_driver.f90'),
                            *objects,'-o',str(exe)],check=True,cwd=tmp)
            run=subprocess.run([str(exe)],input=data,capture_output=True,text=True,check=True,cwd=tmp)
            print(mode,verify(run.stdout,cases),flush=True)
if __name__=='__main__':
    main()
