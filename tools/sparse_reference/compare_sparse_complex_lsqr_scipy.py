#!/usr/bin/env python3
"""Compare complex SciFort LSQR with SciPy and dense least-squares references.

Run: python tools/sparse_reference/compare_sparse_complex_lsqr_scipy.py --mode both --cases 200
No Fortran-package-manager or external BLAS/LAPACK link dependency.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import csr_array, csc_array
from scipy.sparse.linalg import lsqr

ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    'src/kinds.f90', 'src/sparse.f90', 'src/sparse_coo.f90',
    'src/sparse_coo_nd.f90', 'src/sparse_coo_nd_complex.f90',
    'src/sparse_complex.f90', 'src/sparse_linalg_complex_lsqr.f90',
)


def compile_driver(work: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    opts = ['gfortran', '-std=f2018', *flags, '-J', str(work), '-I', str(work)]
    objects = []
    for path in SOURCE:
        obj = work / (Path(path).stem+'.o')
        subprocess.run([*opts, '-c', str(ROOT/path), '-o', str(obj)], cwd=work,
                       capture_output=True, text=True, check=True)
        objects.append(obj)
    exe = work/'complex_lsqr_driver'
    subprocess.run([*opts, str(ROOT/'tools/sparse_reference/complex_lsqr_driver.f90'),
                    *map(str, objects), '-o', str(exe)], cwd=work, capture_output=True,
                   text=True, check=True)
    return exe


def generate_cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    cases = []
    for i in range(count):
        m = int(rng.integers(0 if i % 31 == 0 else 2, 32))
        n = int(rng.integers(0 if i % 29 == 0 else 2, 28))
        fmt = 1+i%2
        cap = 6*max(m,n,1)
        damp = [0.0, 0.15, 0.9, 0.01][i%4]
        usex0 = i%3 != 0
        a = (rng.normal(size=(m,n))+1j*rng.normal(size=(m,n)))*0.25
        a[rng.random((m,n)) < (0.2 if i%3 else 0.65)] = 0.0
        if i%11 == 0 and n > 1:
            a[:, -1] = a[:, 0]  # Exactly rank deficient.
        if i%13 == 0 and n > 0:
            a[:, -1] = 0.0  # Zero column.
        if i%17 == 0 and m > 1:
            a[-1, :] = a[0,:]  # Repeated row.
        b = rng.normal(size=m)+1j*rng.normal(size=m)
        if i%7 == 0:
            b = a @ (rng.normal(size=n)+1j*rng.normal(size=n))
        x0 = rng.normal(size=n)*0.2+1j*rng.normal(size=n)*0.2
        if i%19 == 0:
            b = np.zeros(m,dtype=complex)
            x0 = np.zeros(n,dtype=complex)
        cases.append((a,b,x0,fmt,cap,float(damp),int(usex0)))
    return cases


def serialize(cases):
    lines = [str(len(cases))]
    for a,b,x0,fmt,cap,damp,usex0 in cases:
        m,n = a.shape
        lines.append(f'{m} {n} {fmt} {cap} {usex0} {damp:.17e} 1e-12 1e-12 1e14')
        for value in np.concatenate((a.ravel(),b,x0)):
            lines.append(f'{value.real:.17e} {value.imag:.17e}')
    return '\n'.join(lines)+'\n'


def check_results(output, cases):
    lines = output.splitlines()
    p = 0
    largest_vs_scipy = 0.0
    largest_vs_dense = 0.0
    largest_stationarity = 0.0
    stop_codes = set()
    for index,(a,b,x0,fmt,cap,damp,usex0) in enumerate(cases):
        m,n = a.shape
        h = lines[p].split()
        p += 1
        info, itn, nx = map(int, h[:3])
        r1, r2, gradnorm = map(float, h[3:])
        if nx != n or info not in (0,1,2,3,7) or itn > cap:
            raise AssertionError(f'case {index}: LSQR status={info}, size={nx}/{n}, iter={itn}')
        got = np.array([complex(*map(float,lines[p+j].split())) for j in range(n)])
        p += n
        start = x0 if usex0 else np.zeros(n,dtype=complex)
        resid = b-a@got
        augmented = np.hypot(np.linalg.norm(resid),damp*np.linalg.norm(got-start))
        grad = a.conj().T@resid - damp**2*(got-start)
        stationarity = np.linalg.norm(grad)
        if not np.allclose([r1,r2,gradnorm],
                           [np.linalg.norm(resid),augmented,stationarity],atol=3e-11,rtol=1e-9):
            raise AssertionError(f'case {index}: bad reported norms')
        if n:
            if damp:
                left = np.vstack((a,damp*np.eye(n)))
                right = np.concatenate((b-a@start,np.zeros(n,dtype=complex)))
            else:
                left, right = a, b-a@start
            reference_dense = start+np.linalg.lstsq(left,right,rcond=None)[0]
        else:
            reference_dense = np.zeros(0,dtype=complex)
        dense_err = np.linalg.norm(got-reference_dense)
        dense_tol = 1e-7*(1+np.linalg.norm(reference_dense))
        if dense_err > dense_tol:
            raise AssertionError(f'case {index}: dense solution diff={dense_err:.3e}, '
                                 f'limit={dense_tol:.3e}, shape={a.shape}, damp={damp}, info={info}')
        largest_vs_dense = max(largest_vs_dense,dense_err)
        largest_stationarity = max(largest_stationarity,stationarity)
        if m and n:
            sparse = csr_array(a) if fmt == 1 else csc_array(a)
            scipy_result = lsqr(sparse,b,damp=damp,atol=1e-12,btol=1e-12,
                                conlim=1e14,iter_lim=cap,x0=start)
            scipy_err = np.linalg.norm(got-scipy_result[0])
            if scipy_err > 2e-7*(1+np.linalg.norm(scipy_result[0])):
                raise AssertionError(f'case {index}: SciPy diff={scipy_err:.3e}')
            largest_vs_scipy = max(largest_vs_scipy,scipy_err)
        stop_codes.add(info)
    assert p == len(lines), (p,len(lines))
    return largest_vs_scipy,largest_vs_dense,largest_stationarity,sorted(stop_codes)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    parser.add_argument('--cases', type=int, default=200)
    parser.add_argument('--seed', type=int, default=20261106)
    args = parser.parse_args()
    if args.cases <= 0:
        parser.error('cases must be positive')
    cases = generate_cases(args.cases,args.seed)
    stream = serialize(cases)
    modes = ('checked','optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-complex-lsqr-') as directory:
            exe = compile_driver(Path(directory),mode)
            process = subprocess.run([str(exe)], input=stream, text=True,
                                     capture_output=True, check=True)
            vs_scipy,vs_dense,stationarity,stop_codes = check_results(process.stdout,cases)
            print(f'{mode}: PASS {len(cases)} complex CSR/CSC LSQR comparisons; '
                  f'vs SciPy={vs_scipy:.3e}, vs dense={vs_dense:.3e}, '
                  f'max normal residual={stationarity:.3e}, stop codes={stop_codes}',flush=True)


if __name__ == '__main__':
    main()
