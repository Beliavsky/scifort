#!/usr/bin/env python3
"""Compare SciFort's randomized interpolative decomposition with SciPy/NumPy.

Usage: python test/compare_linalg_interpolative_random_scipy.py /path/to/driver

Only approximate reconstructions are compared: SciPy and SciFort use different
random streams, selected-column permutations, and tolerance norms.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
import scipy
from scipy.linalg import interpolative as sli


def generate_cases():
    rng = np.random.default_rng(20261110)
    result = []
    for m in (0, 1, 2, 5, 9, 15, 27):
        for n in (0, 1, 3, 8, 13, 25):
            for complex_case in (False, True):
                for variant in range(4):
                    a = rng.standard_normal((m,n))
                    if complex_case:
                        a = a + 1j*rng.standard_normal((m,n))
                    if variant == 1 and m and n:
                        rank = min(2, min(m,n))
                        b = rng.standard_normal((m,rank))
                        c = rng.standard_normal((rank,n))
                        if complex_case:
                            b = b+1j*rng.standard_normal(b.shape)
                            c = c+1j*rng.standard_normal(c.shape)
                        a = b@c
                    elif variant == 2:
                        a = np.zeros_like(a)
                    elif variant == 3:
                        a = a* (1e-80 if n%2 else 1e80)
                    r = min(m,n)
                    if variant in (0,2):
                        mode, rank, eps = 0, r//2, 0.0
                    elif variant == 1:
                        mode, rank, eps = 1, 0, 1.0e-5
                    else:
                        mode, rank, eps = 1, 0, 0.1
                    result.append((a,mode,rank,eps,8,int(rng.integers(1,10000000))))
    for j in range(80):
        m = int(rng.integers(20,55))
        n = int(rng.integers(20,60))
        rank = int(rng.integers(1,7))
        is_complex = j%2 == 0
        b = rng.standard_normal((m,rank))
        c = rng.standard_normal((rank,n))
        if is_complex:
            b = b+1j*rng.standard_normal(b.shape)
            c = c+1j*rng.standard_normal(c.shape)
        a = b@c
        if j%3 == 0:
            noise = rng.standard_normal((m,n))
            if is_complex:
                noise = noise + 1j*rng.standard_normal((m,n))
            a = a+noise*(1e-3*np.linalg.norm(a)/np.sqrt(m*n))
        mode, k, eps = (0,rank,0.0) if j%3 != 2 else (1,0,1e-2)
        result.append((a,mode,k,eps,8,int(rng.integers(1,10000000))))
    return result


def run(executable):
    cases = generate_cases()
    records = [str(len(cases))]
    for a, mode, rank, eps, oversampling, seed in cases:
        m,n = a.shape
        dtype = int(np.iscomplexobj(a))
        records.append(f"{dtype} {m} {n} {mode} {rank} {eps:.17e} {oversampling} {seed}")
        for i in range(m):
            for j in range(n):
                x = a[i,j]
                records.append(f"{x.real:.17e} {x.imag:.17e}")
    process = subprocess.run([executable], input='\n'.join(records)+'\n', text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             check=True, timeout=200)
    output = iter(process.stdout.splitlines())
    max_relative_difference_from_optimal = 0.0
    max_relative_tolerance = 0.0
    max_reconstruction_error = 0.0
    count_scipy = 0
    skipped_scipy = 0
    for case_idx, (a,mode,rank,eps,oversampling,seed) in enumerate(cases,1):
        m,n = a.shape
        parts = next(output).split()
        assert parts[:2] == ['CASE',str(case_idx)],(case_idx,parts)
        status,k = map(int,parts[2:])
        assert status == 0,(case_idx,status,a.shape,mode,rank,eps)
        assert 0 <= k <= min(m,n)
        if mode == 0:
            assert k == rank,(case_idx,k,rank)
        idx = np.fromstring(next(output),dtype=int,sep=' ')
        assert len(idx) == n,(case_idx,'idx',len(idx),n)
        assert sorted(idx.tolist()) == list(range(1,n+1)),(case_idx,idx)
        proj = np.empty((k,n-k),dtype=a.dtype)
        for row in range(k):
            vals = np.fromstring(next(output),sep=' ')
            if np.iscomplexobj(a):
                assert len(vals) == 2*(n-k),(case_idx,'complex proj',len(vals),k,n)
                proj[row,:] = vals[::2]+1j*vals[1::2]
            else:
                assert len(vals) == n-k,(case_idx,'real proj',len(vals),k,n)
                proj[row,:] = vals
        assert np.all(np.isfinite(proj)),case_idx
        indices = idx-1
        b = a[:,indices[:k]]
        reconstructed = np.zeros_like(a)
        reconstructed[:,indices[:k]] = b
        reconstructed[:,indices[k:]] = b @ proj
        assert np.array_equal(reconstructed[:,indices[:k]],b),(case_idx,'skeleton')
        denominator = np.linalg.norm(a)
        err = np.linalg.norm(a-reconstructed)/denominator if denominator else 0.0
        assert np.isfinite(err),(case_idx,err)
        max_reconstruction_error = max(max_reconstruction_error,err)
        if mode == 1:
            assert err <= max(5e-9,eps*1.000001),(case_idx,'tol',err,eps,k,m,n)
            max_relative_tolerance = max(max_relative_tolerance,err)
        if mode == 0:
            sv = np.linalg.svd(a,compute_uv=False)
            optimum = np.linalg.norm(sv[k:])/denominator if denominator else 0.0
            if optimum > 1e-8:
                ratio = err/optimum
                max_relative_difference_from_optimal = max(max_relative_difference_from_optimal,ratio)
                assert ratio < 6.0,(case_idx,'fixed quality',ratio,err,optimum)
            if k == min(m,n):
                assert err < 2e-8,(case_idx,'full rank',err)
        if m == 0 or n == 0:
            continue
        if mode == 0 and k == 0:
            continue
        try:
            if mode == 0:
                si,sp = sli.interp_decomp(a,k,rand=True,rng=np.random.default_rng(seed))
                sk = k
            else:
                sk,si,sp = sli.interp_decomp(a,eps,rand=True,rng=np.random.default_rng(seed))
            sci_rec = sli.reconstruct_matrix_from_id(a[:,si[:sk]],si,sp)
            sci_err = np.linalg.norm(a-sci_rec)/denominator if denominator else 0.0
            assert np.isfinite(sci_err)
            assert len(si) == n
            count_scipy += 1
        except (ValueError, RuntimeError, np.linalg.LinAlgError, TypeError, IndexError, ZeroDivisionError):
            skipped_scipy += 1
    try:
        more = next(output)
        raise AssertionError(f'extra output from Fortran: {more!r}')
    except StopIteration:
        pass
    print(f'SciPy version: {scipy.__version__}; cases passed: {len(cases)}; '
          f'SciPy randomized calls succeeded: {count_scipy}; reference errors: {skipped_scipy}')
    print(f'Max relative reconstruction error (all ranks): {max_reconstruction_error:.3e}; '
          f'max certified adaptive error: {max_relative_tolerance:.3e}')
    print(f'Max fixed-rank Frobenius error divided by optimal SVD tail: '
          f'{max_relative_difference_from_optimal:.6f}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: compare_linalg_interpolative_random_scipy.py DRIVER')
    run(sys.argv[1])
