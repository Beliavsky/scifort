#!/usr/bin/env python3
"""Cross-check SciFort LOBPCG vs SciPy and independent dense eigenproblems.

Finite small matrices with diagonal SPD generalized B allow independent
numpy/scipy.linalg.eigh reference. Eigenvector orientation and repeated
subspace bases are not compared entrywise; true residuals/orthogonality are.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import warnings
import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import lobpcg as scipy_lobpcg

ROOT = Path(__file__).resolve().parents[2]


def cases(count: int, seed: int):
    rng=np.random.default_rng(seed)
    for case in range(count):
        mode=case % 6
        n=int(rng.integers(11, 40))
        k=int(rng.integers(1, min(4,n-2)))
        ny=1 if mode==3 else 0
        x0=rng.normal(size=(n,k))
        a=rng.normal(size=(n,n))*0.035
        mask=rng.random((n,n))<0.11
        a=np.tril(a*mask,-1)
        a=a+a.T
        d=rng.uniform(1.2,5.5,size=n)
        # Keep eigenvalues separated for these convergence tests.
        d += 0.1*np.arange(n)/n
        if mode in (1,5): d=-d
        np.fill_diagonal(a,d)
        b=np.diag(rng.uniform(0.9,1.7,size=n))
        if mode in (2,3):
            # Introduce off-diagonal mass terms while preserving SPD.
            cross=rng.normal(scale=0.008,size=(n,n))
            cross=np.tril(cross*(rng.random((n,n))<0.12),-1)
            b=b+cross+cross.T
        m=np.diag(1.0/np.abs(np.diag(a)))
        y=rng.normal(size=(n,ny)) if ny else np.empty((n,0))
        if mode==3:
            # The B-orthogonal constraint complement is null((B*Y)^T).
            z=linalg.null_space(y.T@b)
            exact=linalg.eigh(z.T@a@z,z.T@b@z,eigvals_only=True)
        elif mode==2:
            exact=linalg.eigh(a,b,eigvals_only=True)
        else:
            exact=np.linalg.eigvalsh(a)
        largest=mode in (0,2,4)
        expected=exact[-k:][::-1] if largest else exact[:k]
        yield dict(n=n,k=k,ny=ny,mode=mode,a=a,b=b,m=m,y=y,x0=x0,
                   expected=expected,largest=largest)


def native_input(data):
    lines=[str(len(data))]
    for c in data:
        lines.append(f"{c['n']} {c['k']} {c['ny']} {c['mode']}")
        for field in ('a','b','m','x0'):
            lines.extend(' '.join(f'{v:.17g}' for v in row) for row in c[field])
        if c['ny']:
            lines.extend(' '.join(f'{v:.17g}' for v in row) for row in c['y'])
    return '\n'.join(lines)+'\n'


def run(mode, count, seed):
    dataset=list(cases(count,seed))
    with tempfile.TemporaryDirectory(prefix='scifort-lobpcg-') as temp:
        d=Path(temp)
        flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
        cmd=['gfortran','-std=f2018',*flags,'-J',str(d),'-I',str(d)]
        objects=[]
        for src in ('src/kinds.f90','src/sparse.f90','src/optimize/lsmr.f90',
                    'src/sparse_linalg.f90'):
            out=d/(Path(src).stem+'.o')
            subprocess.run([*cmd,'-c',str(ROOT/src),'-o',str(out)],cwd=d,check=True,
                stdout=subprocess.DEVNULL)
            objects.append(str(out))
        exe=d/'driver'
        subprocess.run([*cmd,str(ROOT/'tools/sparse_linalg_reference/sparse_lobpcg_reference_driver.f90'),
                        *objects,'-o',str(exe)],cwd=d,check=True)
        process=subprocess.run([str(exe)],input=native_input(dataset),capture_output=True,
                               text=True,timeout=180,check=True)
        lines=process.stdout.splitlines()
        if len(lines)!=len(dataset):
            raise AssertionError(f'expected {len(dataset)} lines, got {len(lines)}: {process.stdout[:500]}')
        worst=0.0; successful=0; nonconverged=0; scipy_matches=0
        for case_index,(case,line) in enumerate(zip(dataset,lines)):
            n=case['n'];k=case['k'];mode=case['mode']
            try:
                values=np.fromstring(line,sep=' ')
            except ValueError as e:
                raise AssertionError(f'case {case_index}, mode {mode}, n {n}, line={line[:800]!r}') from e
            if len(values)!=2+2*k+2+n*k:
                raise AssertionError((mode,n,k,line))
            info,steps=map(int,values[:2]);w=values[2:2+k];r=values[2+k:2+2*k]
            orth,constr=values[2+2*k:4+2*k]
            vec=values[4+2*k:].reshape((k,n)).T
            if info<0:
                raise AssertionError(f'numerical breakdown: {mode=} {info=} {n=}')
            if np.max(r) < 1e-7 and info!=0:
                raise AssertionError('correct residual but not marked success')
            # Any result, even unsuccessful, must have truthful residuals and orthogonal eigenvectors.
            assert orth<1e-7 and constr<1e-7, (mode,n,orth,constr)
            amat=case['a']
            bmat=case['b'] if mode in (2,3) else np.eye(n)
            gram=vec.T@bmat@vec
            assert np.max(np.abs(gram-np.eye(k)))<1e-7, (mode,n,'independent B-orthogonality')
            residual=amat@vec-(bmat@vec)*w[None,:]
            if mode==3:
                y=case['y']
                residual=residual-(bmat@y)@np.linalg.solve(y.T@bmat@y,y.T@residual)
            verified=np.linalg.norm(residual,axis=0)
            assert np.max(np.abs(verified-r))<3e-9, (mode,n,verified,r)
            if info == 0:
                successful+=1
                err=np.max(np.abs(w-case['expected']))
                worst=max(worst,err)
                assert err<1e-5, (mode,n,w,case['expected'],r)
                assert max(r)<1.01e-7, (mode,n,r)
            else:
                nonconverged+=1
                if np.max(r)<1.0e-7:
                    raise AssertionError('false nonconvergence')
            # Cross-check a deterministic stratified subset with SciPy.
            # The full set above is always checked against dense eigh.
            if case_index % 5 != 0:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                opts=dict(largest=case['largest'],tol=1e-7,maxiter=160)
                if mode in (2,3): opts['B']=sparse.csr_matrix(case['b'])
                if mode==3: opts['Y']=case['y']
                if mode==4: opts['M']=sparse.csr_matrix(case['m'])
                if mode==5: opts['M']=sparse.csr_matrix(case['m'])
                try:
                    refw,refv=scipy_lobpcg(sparse.csr_matrix(case['a']),case['x0'].copy(),**opts)
                except (ValueError, NotImplementedError, np.linalg.LinAlgError):
                    continue
            if info==0 and np.max(np.abs(refw-case['expected']))<1e-5:
                assert np.max(np.abs(w-refw))<1e-5, (w,refw)
                scipy_matches+=1
        print(f'{mode}: {count} cases, converged={successful}, nonconverged={nonconverged}, '
              f'SciPy agreement={scipy_matches}, max eigenvalue error={worst:.3g}')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--mode',choices=('both','checked','optimized'),default='both')
    p.add_argument('--cases',type=int,default=90)
    p.add_argument('--seed',type=int,default=20261105)
    a=p.parse_args()
    for mode in ('checked','optimized') if a.mode=='both' else (a.mode,):
        run(mode,a.cases,a.seed)

if __name__=='__main__':
    main()
