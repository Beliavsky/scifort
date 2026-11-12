#!/usr/bin/env python3
"""Independent SciFort complex sparse SVD comparisons with NumPy and SciPy.

Builds the Fortran driver with GNU Fortran. Tests singular values, both true
residual equations, left/right orthogonality, CSC/CSR and optional complex
initial vectors. Includes duplicated COO entries and explicitly stored zeros.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.linalg import svds as scipy_svds

ROOT = Path(__file__).resolve().parents[2]


def small_case(rng: np.random.Generator, i: int):
    q = int(rng.integers(2, 22))
    p = q + int(rng.integers(0, 12))
    m, n = (p, q) if i % 2 else (q, p)
    a = np.zeros((m, n), dtype=np.complex128)
    if i % 31 == 0:
        pass  # Include the zero matrix.
    elif i % 17 == 0:
        for j in range(q):
            a[j,j] = (q-j)*np.exp(0.31j*j) if j < q-2 else 0
        a[0,0] = a[1,1] = q+2j
    else:
        a = (rng.normal(size=(m,n)) + 1j*rng.normal(size=(m,n))) * (rng.random((m,n)) < .35)
        if i % 9 == 0:
            if m >= n:
                a[:,-1] = a[:,0]
            else:
                a[-1,:] = a[0,:]
        if i % 13 == 0:
            a *= .01
    k = min(q, 1+i%3)
    which = 'SM' if i%2==0 else 'LM'
    v0 = rng.normal(size=q)+1j*rng.normal(size=q) if i%7==0 else None
    return (a,k,q,q,which,2e-8,i%2+1,v0,False)


def partial_case(rng: np.random.Generator, i: int):
    q = (52,68,92,121)[i%4]
    p = q+24
    m,n = (p,q) if i%2 else (q,p)
    a = np.zeros((m,n),dtype=np.complex128)
    diagonal = rng.uniform(.05,.9,q)
    diagonal[:3] = [16,10,6]
    a[np.arange(q),np.arange(q)] = diagonal*np.exp(1j*rng.uniform(-2,2,q))
    for _ in range(q):
        row = int(rng.integers(m))
        col = int(rng.integers(n))
        a[row,col] += .008 * (rng.normal()+1j*rng.normal())
    return (a,2,32,32,'LM',2e-8,i%2+1,None,True)


def build_one(mode: str, build: Path) -> Path:
    flags = (['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3'])
    cc = ['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
    objects = []
    for name,rel in [
        ('kinds','src/kinds.f90'),
        ('sparse','src/sparse.f90'),
        ('sparse_coo','src/sparse_coo.f90'),
        ('sparse_coo_nd','src/sparse_coo_nd.f90'),
        ('sparse_coo_nd_complex','src/sparse_coo_nd_complex.f90'),
        ('sparse_complex','src/sparse_complex.f90'),
        ('sparse_linalg_complex_svds','src/sparse_linalg_complex_svds.f90'),
    ]:
        obj = build / f'{name}.o'
        subprocess.run([*cc,'-c',str(ROOT/rel),'-o',str(obj)],cwd=build,check=True,
                       capture_output=True,text=True)
        objects.append(obj)
    exe = build / 'complex_sparse_svds_driver'
    subprocess.run([*cc,str(ROOT/'tools/sparse_reference/complex_sparse_svds_driver.f90'),
                    *(str(obj) for obj in objects),'-o',str(exe)],cwd=build,check=True,
                   capture_output=True,text=True)
    return exe


def check_one_mode(mode: str, count: int, partial: int, seed: int):
    rng = np.random.default_rng(seed)
    cases = [small_case(rng,i) for i in range(count)]
    cases.extend(partial_case(rng,i) for i in range(partial))
    inputs = [str(len(cases))]
    for a,k,ncv,cap,which,tol,fmt,v0,_ in cases:
        rr,cc = np.nonzero(a)
        coords = [(int(r)+1,int(c)+1,complex(a[r,c])) for r,c in zip(rr,cc,strict=True)]
        if coords:
            r,c,x = coords[0]
            coords[0] = (r,c,x*.5)
            coords.append((r,c,x*.5))
            coords.append((r,c,0j))
        inputs.append(f'{a.shape[0]} {a.shape[1]} {len(coords)} {k} {ncv} '
                      f'{cap} {which} {tol:.17g} {fmt} {int(v0 is not None)}')
        inputs.extend(f'{r} {c} {x.real:.17g} {x.imag:.17g}' for r,c,x in coords)
        if v0 is not None:
            inputs.extend(f'{x.real:.17g} {x.imag:.17g}' for x in v0)
    with tempfile.TemporaryDirectory(prefix='scifort-complex-svds-') as td:
        build = Path(td)
        exe = build_one(mode,build)
        out = subprocess.run([str(exe)],cwd=build,input='\n'.join(inputs)+'\n',
                             capture_output=True,text=True,check=True,timeout=180)
        lines = iter(out.stdout.splitlines())
        scipy_checked = 0
        max_error = 0.0
        for index,(a,k,ncv,cap,which,tol,fmt,v0,partial_search) in enumerate(cases,1):
            assert next(lines).strip() == f'CASE {index}'
            info,steps,num_values = map(int,next(lines).split())
            if info != 0 or num_values != k:
                raise AssertionError(f'{mode} {index}: info={info} steps={steps} '
                                     f'size={num_values} shape={a.shape} '
                                     f'which={which} partial={partial_search}')
            singular = np.empty(k)
            residuals = np.empty(k)
            u = np.empty((a.shape[0],k),dtype=np.complex128)
            vh = np.empty((k,a.shape[1]),dtype=np.complex128)
            for j in range(k):
                singular[j],residuals[j] = np.fromstring(next(lines),sep=' ')
                left = np.fromstring(next(lines),sep=' ')
                right = np.fromstring(next(lines),sep=' ')
                assert left.size == 2*a.shape[0] and right.size == 2*a.shape[1]
                u[:,j] = left.reshape(-1,2)[:,0] + 1j*left.reshape(-1,2)[:,1]
                vh[j,:] = right.reshape(-1,2)[:,0] + 1j*right.reshape(-1,2)[:,1]
            spectrum = np.linalg.svd(a,compute_uv=False)
            expected = np.sort(spectrum[:k] if which=='LM' else spectrum[-k:])
            np.testing.assert_allclose(singular,expected,rtol=1e-6,atol=2e-8,
                                       err_msg=f'{mode} case {index} spectrum')
            if partial_search or (index%7==0 and which=='LM' and min(a.shape)>k+1 and np.count_nonzero(a)>0):
                scipy_s = scipy_svds(csr_array(a),k=k,which=which,return_singular_vectors=False,
                                    tol=1e-9,solver='arpack',maxiter=3000,rng=np.random.default_rng(seed+index))
                np.testing.assert_allclose(singular,np.sort(scipy_s),rtol=2e-6,atol=3e-8)
                scipy_checked += 1
            v = vh.conj().T
            lres = a@v-u*singular
            rres = a.conj().T@u-v*singular
            true_residual = np.hypot(np.linalg.norm(lres,axis=0),np.linalg.norm(rres,axis=0))
            max_error = max(max_error,float(np.max(true_residual)))
            np.testing.assert_allclose(residuals,true_residual,rtol=1e-4,atol=3e-10)
            if np.max(true_residual)>4e-8*max(1,np.linalg.norm(a,ord=2)):
                raise AssertionError(f'{mode} case {index}: residual={true_residual}')
            np.testing.assert_allclose(u.conj().T@u,np.eye(k),rtol=0,atol=3e-7)
            np.testing.assert_allclose(vh@vh.conj().T,np.eye(k),rtol=0,atol=3e-7)
        assert next(lines,None) is None
        print(f'{mode}: {len(cases)} NumPy complex SVD comparisons passed '
              f'({scipy_checked} independently cross-checked against SciPy svds); '
              f'largest true residual {max_error:.3e}',flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    parser.add_argument('--cases',type=int,default=100)
    parser.add_argument('--partial',type=int,default=16)
    parser.add_argument('--seed',type=int,default=20261106)
    args = parser.parse_args()
    if args.cases<0 or args.partial<0:
        parser.error('counts must be nonnegative')
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        check_one_mode(mode,args.cases,args.partial,args.seed)

if __name__=='__main__':
    main()
