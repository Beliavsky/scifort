#!/usr/bin/env python3
"""Check complex nonsymmetric SciFort eigs against SciPy and NumPy eigenpairs."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigs, ArpackNoConvergence

ROOT=Path(__file__).resolve().parents[2]
MODES=('LM','SM','LR','SR','LI','SI')

def priority(vals,mode):
    vals=np.asarray(vals)
    return {'LM':lambda: np.abs(vals),'SM':lambda: -np.abs(vals),'LR':lambda: vals.real,
            'SR':lambda: -vals.real,'LI':lambda: vals.imag,'SI':lambda: -vals.imag}[mode]()

def cases_for(seed,count):
    rng=np.random.default_rng(seed)
    rows=[str(count)]
    cases=[]
    for t in range(count):
        n=int(rng.integers(8,31))
        k=int(rng.integers(1,4))
        mode=MODES[t%len(MODES)]
        if t%3==0:
            # Real 2x2 rotation blocks and real diagonals: exact imaginary eigenpairs.
            a=np.zeros((n,n))
            i=0
            while i<n:
                if i+1<n and rng.random()<.7:
                    re=float(rng.uniform(-5,5)); im=float(rng.uniform(.2,5))
                    a[i,i]=a[i+1,i+1]=re
                    a[i,i+1]=-im
                    a[i+1,i]=im
                    i+=2
                else:
                    a[i,i]=rng.uniform(-7,7)
                    i+=1
            # Small similarity transformations within a band.
            a += np.triu((rng.random((n,n))<.07)*rng.normal(scale=.15,size=(n,n)),1)
        elif t%3==1:
            # Nonsymmetric tridiagonal matrices with complex spectral components.
            a=np.diag(np.linspace(-5,5,n)+rng.uniform(-.1,.1,n))
            for i in range(n-1):
                a[i,i+1]=rng.uniform(-2,2)
                a[i+1,i]=rng.uniform(-2,2)
        else:
            a=rng.normal(size=(n,n)) * (rng.random((n,n))<.28)
            a += np.diag(rng.uniform(-3,3,n))
        nc=n
        cap=10
        tol=1e-8
        coo=sparse.coo_matrix(a)
        rows.append(f'{n} {coo.nnz} {k} {mode} {tol:.17g} {cap} {nc}')
        rows.extend(f'{i+1} {j+1} {v:.17g}' for i,j,v in zip(coo.row,coo.col,coo.data))
        cases.append((a,k,mode,nc,tol))
    return '\n'.join(rows)+'\n',cases

def reduced_cases(seed,count):
    """Moderately large sparse matrices with identifiable isolated target modes."""
    rng=np.random.default_rng(seed)
    stream=[str(count)]
    cases=[]
    for t in range(count):
        n=int(rng.integers(50,81))
        mode=('LM','LR','SR','LI','SI')[t%5]
        k=2
        a=np.diag(rng.uniform(-.45,.45,n))
        for i in range(n-1):
            a[i,i+1]=float(rng.uniform(-.025,.025))
            a[i+1,i]=float(rng.uniform(-.025,.025))
        if mode in ('LM','LR'):
            a[0,0]=8.0
            a[1,1]=7.0
        elif mode=='SR':
            a[0,0]=-8.0
            a[1,1]=-7.0
        else:
            a[:4,:4]=0.0
            for i,mag in ((0,8.),(2,7.)):
                a[i,i+1]=-mag
                a[i+1,i]=mag
            a[3,4]=.01
            a[4,3]=-.01
        nc=20
        cap=160
        tol=2e-7
        coo=sparse.coo_matrix(a)
        stream.append(f'{n} {coo.nnz} {k} {mode} {tol:.17g} {cap} {nc}')
        stream.extend(f'{i+1} {j+1} {v:.17g}' for i,j,v in zip(coo.row,coo.col,coo.data))
        cases.append((a,k,mode,nc,tol))
    return '\n'.join(stream)+'\n',cases

def verify(output,cases):
    it=iter(output.splitlines())
    counter={'cases':0,'max_residual':0.,'scipy_distinct_real_matrix_selections':0,'scipy_arpack_failures':0}
    for c,(a,k,mode,ncv,tol) in enumerate(cases,1):
        assert next(it).strip()==f'CASE {c}', f'case marker {c}'
        status,steps=map(int,next(it).split()[1:])
        line=next(it).split(); assert line[0]=='LAMBDA'
        seq=np.array(list(map(float,line[1:])))
        vals=seq[::2]+1j*seq[1::2]
        line=next(it).split(); assert line[0]=='RESIDUAL'
        residuals=np.array(list(map(float,line[1:])))
        cols=[]
        for j in range(k):
            line=next(it).split(); assert line[:2]==['VECTOR',str(j+1)]
            seq=np.array(list(map(float,line[2:])))
            cols.append(seq[::2]+1j*seq[1::2])
        vecs=np.array(cols).T
        assert status==0, f'case {c}: eigs status={status} n={a.shape[0]} k={k} mode={mode}'
        dense_vals=np.linalg.eigvals(a)
        np.testing.assert_allclose(np.linalg.norm(vecs,axis=0),1,atol=1e-7)
        actual=np.linalg.norm(a@vecs-vecs*vals[None,:],axis=0)
        np.testing.assert_allclose(residuals,actual,atol=2e-10,rtol=1e-4)
        assert max(actual)<tol*max(1,np.max(np.sum(abs(a),axis=1)))*2, (c,max(actual))
        ref_top=priority(dense_vals,mode)
        boundary=np.sort(ref_top)[-k]
        # For tied real or imag parts, multiple conjugate roots may legitimately
        # be returned; compare the selected spectral priority and global spectrum.
        assert np.min(priority(vals,mode))>=boundary-1e-6, (c,vals,boundary)
        assert all(min(abs(lam-dense_vals))<1e-6 for lam in vals)
        if ncv == a.shape[0] and k < a.shape[0]-2:
            try:
                reference=eigs(sparse.csr_matrix(a),k=k,which=mode,tol=1e-10,
                               maxiter=2000,return_eigenvectors=False)
                # ARPACK may select conjugate pairs rather than the literal
                # largest imaginary part on real nonsymmetric matrices.
                assert all(min(abs(lam-dense_vals)) < 1e-4 for lam in reference)
                if np.min(priority(reference,mode)) < boundary-5e-4:
                    counter['scipy_distinct_real_matrix_selections']+=1
            except ArpackNoConvergence:
                counter['scipy_arpack_failures']+=1
        counter['cases']+=1
        counter['max_residual']=max(counter['max_residual'],float(max(actual)))
    assert next(it,None) is None
    return counter

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--cases',type=int,default=120)
    p.add_argument('--seed',type=int,default=20261104)
    p.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    p.add_argument('--reduced',action='store_true',help='Test ncv<n sparse subspaces')
    args=p.parse_args()
    data,cases=(reduced_cases if args.reduced else cases_for)(args.seed,args.cases)
    for mode in ('checked','optimized') if args.mode=='both' else (args.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-eigs-') as d:
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
            subprocess.run([*cc,str(ROOT/'tools/sparse_linalg_reference/sparse_eigs_reference_driver.f90'),
                            *obj,'-o',str(exe)],cwd=tmp,check=True)
            result=subprocess.run([str(exe)],input=data,text=True,capture_output=True,check=True,cwd=tmp)
            print(mode,verify(result.stdout,cases),flush=True)
if __name__=='__main__':
    main()
