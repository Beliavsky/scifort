#!/usr/bin/env python3
"""Numerical oracle for hard-deflated symmetric Lanczos with real shift-invert."""
from __future__ import annotations

import argparse
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]


def make_case(rng: np.random.Generator, ci: int):
    n = (16, 22, 28)[(ci // 4) % 3]
    k = 2
    mode = ('LM', 'LA', 'SA', 'SM')[ci % 4]
    diagonal = np.linspace(-3.5, 4.0, n) + rng.uniform(-0.03, 0.03, n)
    a = np.diag(diagonal)
    for j in range(n-1):
        a[j,j+1] = a[j+1,j] = rng.uniform(-0.04, 0.04)
    spectrum = np.linalg.eigvalsh(a)
    middle = n//2
    sigma = spectrum[middle-1]+0.31*(spectrum[middle]-spectrum[middle-1])
    p = n if mode == 'SM' else min(n, 16)
    return a, k, p, 60, mode, 1.0e-8, 1, sigma


def expected(a, k, mode, sigma):
    spectrum = np.linalg.eigvalsh(a)
    inverse = 1.0/(spectrum-sigma)
    if mode == 'LM':
        scores = -abs(inverse)
    elif mode == 'LA':
        scores = -inverse
    elif mode == 'SA':
        scores = inverse
    else:
        scores = abs(inverse)
    return np.sort(spectrum[np.argsort(scores)[:k]])


def run(mode, cases, seed):
    flags = (['-O0','-g','-fcheck=all','-fbacktrace']
             if mode == 'checked' else ['-O3'])
    with tempfile.TemporaryDirectory(prefix='scifort-locked-shift-') as tmp:
        build = Path(tmp)
        cc = ['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
        objs = []
        for name, rel in [('kinds','src/kinds.f90'),
                          ('sparse','src/sparse.f90'),
                          ('lsmr','src/optimize/lsmr.f90'),
                          ('sparse_linalg','src/sparse_linalg.f90')]:
            obj = build/(name+'.o')
            subprocess.run([*cc,'-c',str(ROOT/rel),'-o',str(obj)],cwd=build,check=True)
            objs.append(str(obj))
        executable = build/'locked_shift'
        driver = ROOT/'tools/sparse_linalg_reference/sparse_lanczos_locked_reference_driver.f90'
        subprocess.run([*cc,str(driver),*objs,'-o',str(executable)],cwd=build,check=True)
        rng = np.random.default_rng(seed)
        data = [make_case(rng, i) for i in range(cases)]
        lines = [str(len(data))]
        for a,k,p,cap,which,tol,b,sigma in data:
            lines.append(f'{len(a)} {k} {p} {cap} {which} {tol:.17g} {b} 1 {sigma:.17g}')
            lines.append(' '.join(f'{x:.17g}' for x in a.flatten(order='F')))
        result = subprocess.run([str(executable)],input='\n'.join(lines)+'\n',
                                text=True,capture_output=True,cwd=build,check=True)
        output = iter(result.stdout.splitlines())
        restarted = 0
        for i,(a,k,p,cap,which,tol,b,sigma) in enumerate(data, start=1):
            assert next(output).strip() == f'CASE {i}'
            info,cycle,locked,products = map(int,next(output).split())
            pairs = [np.fromstring(next(output),sep=' ') for _ in range(k)]
            eig = np.array([x[0] for x in pairs])
            errors = np.array([x[1] for x in pairs])
            vec = np.column_stack([x[2:] for x in pairs])
            assert info == 0 and locked == k, (i,which,info,cycle,errors)
            target = expected(a,k,which,sigma)
            np.testing.assert_allclose(eig,target,rtol=1e-6,atol=2e-6)
            actual_errors = np.linalg.norm(a@vec-vec*eig[None,:],axis=0)
            np.testing.assert_allclose(actual_errors,errors,rtol=1e-4,atol=1e-10)
            assert max(errors)<tol*max(1,np.linalg.norm(a,np.inf))*1.02
            np.testing.assert_allclose(vec.T@vec,np.eye(k),atol=2e-7)
            scipy_eigs = np.sort(eigsh(csr_matrix(a),k,which=which,sigma=sigma,
                                       tol=1e-10,maxiter=2000,return_eigenvectors=False))
            np.testing.assert_allclose(eig,scipy_eigs,atol=2e-6)
            restarted += cycle > 1
        assert next(output,None) is None
        print(f'{mode}: {cases} hard-deflated shifted comparisons passed; '
              f'{restarted} required multiple cycles',flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    parser.add_argument('--cases',type=int,default=36)
    parser.add_argument('--seed',type=int,default=20261104)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error('cases must be positive')
    for mode in (['checked','optimized'] if args.mode == 'both' else [args.mode]):
        run(mode,args.cases,args.seed)


if __name__ == '__main__':
    main()
