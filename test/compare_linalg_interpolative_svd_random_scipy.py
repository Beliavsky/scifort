#!/usr/bin/env python3
"""Compare SciFort randomized interpolative SVD with NumPy and SciPy.

Usage:
    python test/compare_linalg_interpolative_svd_random_scipy.py ./random_svd_driver

The algorithms use different random generators, so compare approximate
reconstructions, best rank-k SVD error, and orthogonality, never exact vectors.
SciFort's tolerance is Frobenius relative error; SciPy targets spectral error.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
import scipy
from scipy.linalg import interpolative as sli


def cases():
    rng = np.random.default_rng(20261110)
    out = []
    for m in (0, 1, 2, 4, 7, 11, 19, 27):
        for n in (0, 1, 3, 5, 9, 14, 21, 33):
            for complex_data in (False, True):
                for family in range(5):
                    a = rng.standard_normal((m, n))
                    if complex_data:
                        a = a + 1j * rng.standard_normal((m, n))
                    if family == 1 and m and n:
                        k0 = min(3, min(m,n))
                        b = rng.standard_normal((m,k0))
                        c = rng.standard_normal((k0,n))
                        if complex_data:
                            b = b+1j*rng.standard_normal(b.shape)
                            c = c+1j*rng.standard_normal(c.shape)
                        a = b@c
                    if family == 2:
                        a = a * 1e-40
                    if family == 3:
                        a = a * 1e40
                    if family == 4:
                        a = np.zeros_like(a)
                    r = min(m,n)
                    if family == 0:
                        mode, rank, eps = 0, r // 2, 0.0
                    elif family == 1:
                        mode, rank, eps = 0, min(r,3), 0.0
                    elif family == 2:
                        mode, rank, eps = 0, r, 0.0
                    elif family == 3:
                        mode, rank, eps = 1, 0, 0.20
                    else:
                        mode, rank, eps = 1, 0, 1e-5
                    out.append((np.array(a), mode, rank, eps, 8, 123456))
    # Larger numerical low-rank matrices to actually exercise small sketches.
    for j in range(90):
        m = int(rng.integers(20,65))
        n = int(rng.integers(20,75))
        rank = int(rng.integers(1,8))
        dtype_complex = j % 2 == 1
        b = rng.standard_normal((m,rank))
        c = rng.standard_normal((rank,n))
        if dtype_complex:
            b = b+1j*rng.standard_normal(b.shape)
            c = c+1j*rng.standard_normal(c.shape)
        a = b@c
        if j % 4 == 0:
            a = a + (1e-3 * np.linalg.norm(a) / np.sqrt(m*n)) * (
                rng.standard_normal((m,n)) + (1j*rng.standard_normal((m,n)) if dtype_complex else 0)
            )
        if j % 3 == 0:
            mode, k, eps = 1, 0, 1e-2
        elif j % 3 == 1:
            mode, k, eps = 0, rank, 0.0
        else:
            mode, k, eps = 1, 0, 1e-5
        out.append((a,mode,k,eps,10,j*193+7))
    return out


def check(executable: str):
    refs = cases()
    input_lines = [str(len(refs))]
    for a, mode, rank, eps, over, seed in refs:
        m,n = a.shape
        typ = 1 if np.iscomplexobj(a) else 0
        input_lines.append(f"{typ} {m} {n} {mode} {rank} {eps:.17e} {over} {seed}")
        for i in range(m):
            for j in range(n):
                x = a[i,j]
                input_lines.append(f"{x.real:.17e} {x.imag:.17e}")
    out = subprocess.run(
        [executable], input="\n".join(input_lines)+"\n",
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
        timeout=210,
    ).stdout.splitlines()
    it = iter(out)
    max_residual = 0.0
    max_orth = 0.0
    max_optimal_ratio = 0.0
    max_singular_deviation = 0.0
    scipy_checked = 0
    skipped = 0
    for count,(a, mode, rank, eps, over, seed) in enumerate(refs, 1):
        m,n=a.shape
        header = next(it).split()
        assert header[:2]==['CASE',str(count)], (count, header)
        status,k = map(int, header[2:])
        assert status==0, (count, status, a.shape,mode,rank,eps)
        assert 0 <= k <= min(m,n)
        if mode == 0:
            assert k == rank
        s = np.fromstring(next(it),sep=" ")
        assert len(s)==k, (count,'s',len(s),k)
        dtype = a.dtype
        u = np.zeros((m,k),dtype=dtype)
        v = np.zeros((n,k),dtype=dtype)
        for target, rows in ((u,m),(v,n)):
            for i in range(rows):
                vals = np.fromstring(next(it),sep=" ")
                if np.iscomplexobj(a):
                    assert len(vals)==2*k,(count,'complex',i,vals,k)
                    target[i,:]=vals[::2]+1j*vals[1::2]
                else:
                    assert len(vals)==k,(count,'real',i,vals,k)
                    target[i,:]=vals
        assert np.all(np.isfinite(s)) and np.all(np.isfinite(u)) and np.all(np.isfinite(v))
        assert np.all(s>=-1e-14) and np.all(np.diff(s)<=1e-9*max(float(max(s,default=0)),1))
        if k:
            ou=np.linalg.norm(u.conj().T@u-np.eye(k),ord=2)
            ov=np.linalg.norm(v.conj().T@v-np.eye(k),ord=2)
            orth=max(ou,ov)
            max_orth=max(orth,max_orth)
            assert orth<2e-7,(count,'orth',orth,m,n,k)
        approx=(u*s) @ v.conj().T
        denom=np.linalg.norm(a)
        err=np.linalg.norm(a-approx)/denom if denom else 0.0
        max_residual=max(max_residual,err)
        if mode:
            assert err <= max(eps*1.002, 1.5e-7), (count,'adaptive',m,n,k,err,eps)
        else:
            exact_s=np.linalg.svd(a,compute_uv=False)
            opt=np.linalg.norm(exact_s[k:])/denom if denom else 0.0
            if opt > 1e-10:
                max_optimal_ratio=max(max_optimal_ratio,err/opt)
                assert err <= 9.0*opt + 2e-7, (count,'fixed bound',m,n,k,err,opt)
            else:
                assert err<1.5e-7,(count,'fixed lowrank',m,n,k,err,opt)
            if np.linalg.norm(exact_s[:k]):
                sv_dev=np.linalg.norm(s-exact_s[:k])/np.linalg.norm(exact_s[:k])
                max_singular_deviation=max(max_singular_deviation,sv_dev)
                assert sv_dev < 0.9,(count,'sv diff',sv_dev)
        # SciPy stochastic methods need not yield identical indices, ranks or
        # phases. Compare existence, shapes and independent reconstruction.
        if m and n and k and min(m,n) > 1:
            try:
                srank = int(k) if mode == 0 else eps
                su,ss,sv = sli.svd(a,srank,rand=True,rng=np.random.default_rng(seed))
                assert su.shape[0] == m and sv.shape[0] == n
                assert len(ss) == su.shape[1] == sv.shape[1]
                assert np.all(np.isfinite(ss))
                scipy_checked += 1
            except (ValueError, IndexError, FloatingPointError, RuntimeError, np.linalg.LinAlgError):
                skipped += 1
    print(f"SciFort comparisons passed: {len(refs)}")
    print(f"SciPy randomized calls succeeded: {scipy_checked}; skipped by SciPy: {skipped}")
    print(f"Largest relative reconstruction error: {max_residual:.6e}")
    print(f"Largest orthogonality error: {max_orth:.6e}")
    print(f"Largest ratio to optimal rank-k Frobenius error (nontrivial tails): {max_optimal_ratio:.4f}")
    print(f"Largest relative leading-singular-value deviation: {max_singular_deviation:.6e}")
    print(f"SciPy {scipy.__version__}, NumPy {np.__version__}")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: compare_linalg_interpolative_svd_random_scipy.py DRIVER')
    check(sys.argv[1])
