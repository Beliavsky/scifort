#!/usr/bin/env python3
"""Independent SciPy and dense reference comparisons for complex CSR/CSC LSMR.

python tools/sparse_reference/compare_sparse_complex_lsmr_scipy.py --mode both --cases 250
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import csr_array, csc_array
from scipy.sparse.linalg import lsmr

ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    'src/kinds.f90', 'src/sparse.f90', 'src/sparse_coo.f90',
    'src/sparse_coo_nd.f90', 'src/sparse_coo_nd_complex.f90',
    'src/sparse_complex.f90', 'src/sparse_linalg_complex_lsmr.f90',
)


def compile_driver(work: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    opts = ['gfortran', '-std=f2018', *flags, '-J', str(work), '-I', str(work)]
    objs = []
    for src in SOURCE:
        obj = work/(Path(src).stem+'.o')
        subprocess.run([*opts, '-c', str(ROOT/src), '-o', str(obj)],
                       cwd=work, capture_output=True, text=True, check=True)
        objs.append(obj)
    exe = work/'complex_lsmr_driver'
    subprocess.run([*opts, str(ROOT/'tools/sparse_reference/complex_lsmr_driver.f90'),
                    *map(str,objs), '-o', str(exe)], cwd=work,
                   capture_output=True, text=True, check=True)
    return exe


def generate_cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    cases = []
    for i in range(count):
        m = int(rng.integers(0 if i % 32 == 0 else 2, 33))
        n = int(rng.integers(0 if i % 31 == 0 else 2, 29))
        fmt = 1+i%2
        cap = 9*max(m,n,1)
        damp = [0.0,0.12,0.8,0.02][i%4]
        usex0 = int(i%3 != 0)
        a = (rng.normal(size=(m,n))+1j*rng.normal(size=(m,n)))*0.2
        a[rng.random((m,n))<(0.55 if i%3 == 0 else 0.2)] = 0
        if i%11 == 0 and n > 1:
            a[:,-1] = a[:,0]
        if i%13 == 0:
            a[:,-1] = 0
        if i%17 == 0 and m > 1:
            a[-1,:] = a[0,:]
        if i%23 == 0:
            a *= 0
        target = rng.normal(size=n)+1j*rng.normal(size=n)
        b = rng.normal(size=m)+1j*rng.normal(size=m)
        if i%7 == 0:
            b = a @ target
        x0 = (rng.normal(size=n)+1j*rng.normal(size=n))*0.2
        if i%19 == 0:
            b = np.zeros(m,dtype=complex)
            x0 = np.zeros(n,dtype=complex)
        cases.append((a,b,x0,fmt,cap,float(damp),usex0))
    return cases


def serialize(cases):
    lines = [str(len(cases))]
    for a,b,x0,fmt,cap,damp,usex0 in cases:
        m,n = a.shape
        lines.append(f'{m} {n} {fmt} {cap} {usex0} {damp:.17e} 1e-12 1e-12 1e14')
        for z in np.concatenate((a.ravel(),b,x0)):
            lines.append(f'{z.real:.17e} {z.imag:.17e}')
    return '\n'.join(lines)+'\n'


def check_results(out: str, cases):
    lines = out.splitlines()
    p = 0
    largest_scipy = 0.0
    largest_dense = 0.0
    largest_stationarity = 0.0
    stop_codes = set()
    for i,(a,b,x0,fmt,cap,damp,usex0) in enumerate(cases):
        m,n = a.shape
        head = lines[p].split()
        p += 1
        info,itn,nx = map(int,head[:3])
        rnorm,arnorm,anorm,acond,xnorm = map(float,head[3:])
        if nx != n or info not in range(8) or itn > cap:
            raise AssertionError(f'case {i}: stop={info}, size={nx}/{n}, iterations={itn}')
        got = np.array([complex(*map(float,lines[p+j].split())) for j in range(n)])
        p += n
        center = x0 if usex0 else np.zeros(n,complex)
        resid = b - a @ got
        true_r = np.hypot(np.linalg.norm(resid),damp*np.linalg.norm(got-center))
        grad = a.conj().T @ resid - damp**2*(got-center)
        gradnorm = np.linalg.norm(grad)
        if not np.allclose([rnorm,arnorm,xnorm],
                           [true_r,gradnorm,np.linalg.norm(got)],atol=1e-10,rtol=1e-9):
            raise AssertionError(f'case {i}: diagnostics differ: '
                                 f'{[rnorm,arnorm,xnorm]} vs {[true_r,gradnorm,np.linalg.norm(got)]}')
        if not np.isfinite(anorm) or not np.isfinite(acond):
            raise AssertionError(f'case {i}: nonfinite condition diagnostics')
        if damp:
            left = np.vstack((a,damp*np.eye(n)))
            rhs = np.concatenate((b-a@center,np.zeros(n,complex)))
        else:
            left = a
            rhs = b-a@center
        dense = center + np.linalg.lstsq(left,rhs,rcond=None)[0]
        err = np.linalg.norm(got-dense)
        if err > 8e-7*(1+np.linalg.norm(dense)):
            raise AssertionError(f'case {i}: dense residual={err:.3e} shape={a.shape} '
                                 f'damp={damp}, rank={np.linalg.matrix_rank(a)}, status={info}')
        largest_dense = max(largest_dense,err)
        sparse = csr_array(a) if fmt == 1 else csc_array(a)
        sci = lsmr(sparse,b,atol=1e-12,btol=1e-12,conlim=1e14,
                   maxiter=cap,damp=damp,x0=center)[0]
        err = np.linalg.norm(sci-got)
        if err > 9e-7*(1+np.linalg.norm(sci)):
            raise AssertionError(f'case {i}: vs scipy error {err:.3e}, shape={a.shape}, '
                                 f'status={info}')
        largest_scipy = max(largest_scipy,err)
        largest_stationarity = max(largest_stationarity,gradnorm)
        stop_codes.add(info)
    if p != len(lines):
        raise AssertionError(f'Unexpected output trailing lines: {len(lines)-p}')
    return largest_scipy,largest_dense,largest_stationarity,sorted(stop_codes)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    parser.add_argument('--cases',type=int,default=250)
    parser.add_argument('--seed',type=int,default=20261106)
    args = parser.parse_args()
    cases = generate_cases(args.cases,args.seed)
    stream = serialize(cases)
    for mode in (('checked','optimized') if args.mode == 'both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-complex-lsmr-') as temp:
            exe = compile_driver(Path(temp),mode)
            output = subprocess.run([str(exe)],input=stream,text=True,capture_output=True,check=True)
            sci,dense,stat,codes = check_results(output.stdout,cases)
            print(f'{mode}: PASS {len(cases)} complex CSR/CSC LSMR cases; '
                  f'max vs SciPy={sci:.3e}, dense={dense:.3e}, stationarity={stat:.3e}; '
                  f'stop codes={codes}',flush=True)


if __name__ == '__main__':
    main()
