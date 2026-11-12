#!/usr/bin/env python3
"""Compare SciFort complex triangular/eigsh routines with SciPy/NumPy.

Compiles isolated Fortran modules; does not require FPM or BLAS at link time.
Run: python tools/sparse_reference/compare_sparse_complex_triangular_eigen_scipy.py --mode both --cases 160
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.linalg import solve_triangular
from scipy.sparse import csr_array
from scipy.sparse.linalg import spsolve_triangular, eigsh

ROOT = Path(__file__).resolve().parents[2]
SOURCE = [
    'src/kinds.f90', 'src/sparse.f90', 'src/sparse_coo.f90',
    'src/sparse_coo_nd.f90', 'src/sparse_coo_nd_complex.f90',
    'src/sparse_complex.f90', 'src/sparse_linalg_complex_triangular.f90',
    'src/sparse_linalg_complex_eigen.f90',
]


def compile_driver(work: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    opts = ['gfortran', '-std=f2018', *flags, '-J', str(work), '-I', str(work)]
    objects = []
    for f in SOURCE:
        obj = work / (Path(f).stem + '.o')
        subprocess.run([*opts, '-c', str(ROOT / f), '-o', str(obj)], cwd=work, check=True,
                       capture_output=True, text=True)
        objects.append(obj)
    exe = work / 'complex_triangular_eigen_driver'
    subprocess.run([*opts, str(ROOT / 'tools/sparse_reference/complex_triangular_eigen_driver.f90'),
                    *map(str, objects), '-o', str(exe)], cwd=work, check=True,
                   capture_output=True, text=True)
    return exe


def pair_lines(matrix: np.ndarray) -> list[str]:
    return [f'{z.real:.17e} {z.imag:.17e}' for z in matrix.ravel(order='C')]


def generate_cases(rng: np.random.Generator, count: int):
    cases = []
    modes = ('N', 'T', 'H')
    whiches = ('LA', 'SA', 'LM', 'SM')
    for case_id in range(count):
        n = int(rng.integers(0, 22))
        nrhs = int(rng.integers(0, 5))
        lower = bool(rng.integers(0, 2))
        unit = bool(rng.integers(0, 2))
        trans = modes[case_id % len(modes)]
        fmt = 1 + (case_id % 2)
        d = (rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))*0.25
        d[rng.random((n,n)) < 0.5] = 0
        d = np.tril(d) if lower else np.triu(d)
        if n:
            np.fill_diagonal(d, 0 if unit else 1.0+rng.uniform(0.3,1.7,n) + 1j*rng.uniform(-0.7,0.7,n))
        effective = d.copy()
        if n and unit:
            np.fill_diagonal(effective, 1)
        xtrue = rng.normal(size=(n,nrhs))+1j*rng.normal(size=(n,nrhs))
        op = {'N':effective,'T':effective.T,'H':effective.conj().T}[trans]
        b = op@xtrue
        reference = solve_triangular(effective,b,lower=lower,unit_diagonal=unit,
                                     trans={'N':0,'T':1,'H':2}[trans])
        if trans == 'N' and n and nrhs:
            sx = spsolve_triangular(csr_array(effective),b,lower=lower,unit_diagonal=unit)
            if not np.allclose(sx,reference,atol=1e-11,rtol=1e-11):
                raise AssertionError('SciPy triangular/reference mismatch')
        cases.append(('T',n,nrhs,int(lower),int(unit),trans,fmt,d,b,reference))
        if n < 2:
            continue
        ne = int(rng.integers(2, min(19,n+2)))
        # Hermitian, both positive and negative eigenvalues, with genuinely complex off-diagonals.
        re = rng.normal(size=(ne,ne))
        im = rng.normal(size=(ne,ne))
        re[np.abs(rng.normal(size=(ne,ne))) > 1.0] = 0
        h = (re+1j*im)*0.17
        h[rng.random((ne,ne)) < 0.25] = 0
        h = np.triu(h,1)
        h = h+h.conj().T
        h += np.diag(rng.normal(size=ne)*1.5)
        k = int(rng.integers(1,min(4,ne)))
        which = whiches[case_id%4]
        w, _ = np.linalg.eigh(h)
        order = np.argsort(-w if which == 'LA' else w if which == 'SA'
                           else -np.abs(w) if which == 'LM' else np.abs(w), kind='stable')
        expected = w[order[:k]]
        if ne >= 8 and k <= 2 and which != 'SM' and case_id%17 == 0:
            sw = eigsh(csr_array(h),k=k,which=which,tol=1e-10,return_eigenvectors=False)
            if np.max(np.min(np.abs(sw[:,None]-w[None,:]),axis=1)) > 1e-7:
                raise AssertionError('SciPy eigsh differs from NumPy reference')
        cases.append(('E',ne,k,which,ne,1e-10,fmt,h,expected))
    return cases


def stream_input(cases):
    data = [str(len(cases))]
    for case in cases:
        if case[0] == 'T':
            _, n, nrhs, lower, unit, trans, fmt, matrix, rhs, _ = case
            data.extend(('T',f'{n} {nrhs} {lower} {unit} {trans} {fmt}'))
            data.extend(pair_lines(matrix))
            data.extend(pair_lines(rhs))
        else:
            _, n, k, which, ncv, tol, fmt, matrix, _ = case
            data.extend(('E',f'{n} {k} {which} {ncv} {tol:.17e} {fmt}'))
            data.extend(pair_lines(matrix))
    return '\n'.join(data)+'\n'


def verify(output: str, cases):
    lines = output.splitlines()
    p = 0
    max_err_t = max_err_e = max_res = 0.0
    eig_count = tri_count = 0
    for case in cases:
        h = lines[p].split()
        p += 1
        if h[0] != case[0]:
            raise AssertionError((h, case))
        info = int(h[1])
        if info != 0:
            raise AssertionError(f'{case[0]} info={info}, case={case[:7]}')
        if case[0] == 'T':
            _, n, nrhs, _, _, _, _, _, _, expected = case
            if (int(h[2]),int(h[3])) != (n,nrhs):
                raise AssertionError('triangular output shape')
            got = np.zeros((n,nrhs),dtype=complex)
            for i in range(n):
                for j in range(nrhs):
                    r, im = map(float,lines[p].split())
                    p += 1
                    got[i,j] = r+1j*im
            err = np.max(np.abs(got-expected), initial=0)
            if err > 2e-10:
                raise AssertionError(f'triangular residual {err:.3e} case={case[:7]}')
            max_err_t = max(max_err_t,err)
            tri_count += 1
        else:
            _, n, k, which, _, tol, _, matrix, expected = case
            if (int(h[2]),int(h[3])) != (k,n):
                raise AssertionError('eigen output shape')
            got = np.zeros(k)
            vectors = np.zeros((n,k),dtype=complex)
            resids = []
            for j in range(k):
                got[j], residual = map(float,lines[p].split())
                resids.append(residual)
                p += 1
                for i in range(n):
                    r, im = map(float,lines[p].split())
                    p += 1
                    vectors[i,j] = r+1j*im
            err = np.max(np.abs(got-expected), initial=0)
            rr = np.max(np.linalg.norm(matrix@vectors-vectors*got,axis=0),initial=0)
            oo = np.max(np.abs(vectors.conj().T@vectors-np.eye(k)),initial=0)
            if err > 3e-8 or rr > 2e-8 or oo > 3e-8:
                raise AssertionError(f'{which} eigen errors {err:.3e},{rr:.3e},{oo:.3e} case={case[:7]}')
            if np.max(np.abs(np.array(resids)-np.linalg.norm(matrix@vectors-vectors*got,axis=0)),initial=0)>1e-11:
                raise AssertionError('incorrect reported true residual')
            max_err_e = max(max_err_e,err)
            max_res = max(max_res,rr)
            eig_count += 1
    assert p == len(lines), (p,len(lines))
    return tri_count, eig_count, max_err_t, max_err_e, max_res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    parser.add_argument('--cases', type=int, default=160)
    parser.add_argument('--seed', type=int, default=20261106)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error('cases must be positive')
    rng = np.random.default_rng(args.seed)
    cases = generate_cases(rng,args.cases)
    input_data = stream_input(cases)
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-complex-eigen-') as work:
            exe = compile_driver(Path(work),mode)
            process = subprocess.run([str(exe)],input=input_data,text=True,capture_output=True,check=True)
            tri,eig,mt,me,mr = verify(process.stdout,cases)
            print(f'{mode}: PASS {tri} triangular + {eig} Hermitian eigen comparisons; '
                  f'max triangular error {mt:.3e}, max eigenvalue error {me:.3e}, '
                  f'max true eigen residual {mr:.3e}',flush=True)


if __name__ == '__main__':
    main()
