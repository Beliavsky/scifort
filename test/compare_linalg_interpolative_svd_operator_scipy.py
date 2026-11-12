#!/usr/bin/env python3
"""Validate SciFort matrix-free randomized ID-SVD against SciPy / NumPy.

Usage: python test/compare_linalg_interpolative_svd_operator_scipy.py DRIVER

Fortran and SciPy use independent random generators and numerical algorithms;
tests compare rank, orthogonality, singular values, and approximation quality,
not the phases of individual singular vectors or the random sketch itself.
"""

from __future__ import annotations

import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.linalg import interpolative
from scipy.sparse.linalg import aslinearoperator


def fortran_svd(exe: Path, a: np.ndarray, target: int | float, seed: int = 188137):
    m, n = a.shape
    cmplx = np.iscomplexobj(a)
    mode = int(isinstance(target, float))
    rows = [f"{m} {n} {int(cmplx)} {mode} {target} {seed}\n"]
    for i in range(m):
        if cmplx:
            rows.append(" ".join(f"({z.real:.17e},{z.imag:.17e})" for z in a[i]) + "\n")
        else:
            rows.append(" ".join(f"{x:.17e}" for x in a[i]) + "\n")
    p = subprocess.run([str(exe)], input="".join(rows), capture_output=True,
                       text=True, timeout=30, check=False)
    if p.returncode != 0:
        raise RuntimeError(f"Driver failed: {p.stderr}\n{p.stdout[:800]}")
    lines = p.stdout.splitlines()
    header = lines[0].split()
    assert header[0] == 'RESULT', lines[:2]
    status, k, est = int(header[1]), int(header[2]), float(header[3])
    assert status == 0, (status, m, n, target, lines[:4])
    calls = tuple(map(int, lines[1].split()[1:]))
    dtype = np.complex128 if cmplx else np.float64
    u, v, s = np.zeros((m, k), dtype=dtype), np.zeros((n, k), dtype=dtype), np.zeros(k)
    expected = 2 + k*(1+m+n)
    assert len(lines) == expected, (len(lines), expected)
    for line in lines[2:]:
        data = line.split()
        if data[0] == 'S':
            s[int(data[1])-1] = float(data[2])
        elif data[0] in ('U', 'V'):
            arr = u if data[0] == 'U' else v
            i, j = int(data[1])-1, int(data[2])-1
            arr[i,j] = (float(data[3]) + 1j*float(data[4])) if cmplx else float(data[3])
        else:
            raise AssertionError(data)
    return u, s, v, k, est, calls


def relative_error(a: np.ndarray, approximated: np.ndarray) -> float:
    if not a.size:
        return 0.0
    scale = np.max(np.abs(a))
    if scale == 0:
        return float(np.linalg.norm(approximated))
    return float(np.linalg.norm((a-approximated)/scale)/np.linalg.norm(a/scale))


def check_orthonormal(q):
    k = q.shape[1]
    if k:
        return float(np.max(abs(q.conj().T @ q - np.eye(k))))
    return 0.0


def test(exe: Path):
    rng = np.random.default_rng(20261110)
    count, scipy_ok, scipy_failures, scipy_bad = 0, 0, 0, 0
    scipy_failure_causes: dict[str, int] = {}
    max_reconstruction_ratio = 0.0
    max_ortho = 0.0
    max_singular_relative_error = 0.0
    max_adaptive_error = 0.0
    max_scipy_residual = 0.0
    total_forward, total_adjoint = 0, 0
    for m, n in [(0,0),(0,4),(4,0),(3,5),(5,3),(6,8),(8,6),(9,9),
                 (12,18),(18,12),(16,25),(25,16),(20,29)]:
        for cmplx in (False, True):
            for variant in range(4):
                x = rng.normal(size=(m,n))
                if cmplx:
                    x = x + 1j*rng.normal(size=(m,n))
                if variant == 0:
                    a = x
                elif variant == 1:
                    r = min(3,m,n)
                    left = rng.normal(size=(m,r))
                    right = rng.normal(size=(r,n))
                    if cmplx:
                        left = left + 1j*rng.normal(size=left.shape)
                        right = right + 1j*rng.normal(size=right.shape)
                    a = left @ right
                elif variant == 2:
                    a = x*1.0e-80
                else:
                    a = x*1.0e80
                maxrank = min(m,n)
                targets: list[int | float] = list(dict.fromkeys([0,min(1,maxrank),
                    min(3,maxrank),maxrank])) + [0.1,1e-5]
                for target in targets:
                    u,s,v,k,est,calls = fortran_svd(exe,a,target)
                    approx = (u*s) @ v.conj().T
                    if isinstance(target,int):
                        assert k == target, (target,k)
                        if target:
                            assert calls[0] == target, (calls,target)
                    else:
                        assert np.isfinite(est) and est <= max(target,1e-12), (target,est)
                    assert 0 <= k <= maxrank
                    assert s.shape == (k,)
                    assert np.all(np.isfinite(s)) and np.all(s >= -1e-12)
                    assert not np.any(np.diff(s) > 1e-11*max(1.0,float(np.max(s,initial=0.0))))
                    assert np.all(np.isfinite(u)) and np.all(np.isfinite(v))
                    ortho = max(check_orthonormal(u),check_orthonormal(v))
                    max_ortho = max(max_ortho,ortho)
                    assert ortho < 1e-8, (m,n,variant,target,ortho)
                    err = relative_error(a,approx)
                    if isinstance(target,float):
                        max_adaptive_error = max(max_adaptive_error,err)
                        assert err <= max(1.1*target,2e-8), (m,n,cmplx,variant,err,k)
                    else:
                        if a.size and np.max(abs(a)):
                            spectrum = np.linalg.svd(a/np.max(abs(a)),compute_uv=False)
                            best = float(np.linalg.norm(spectrum[k:])/np.linalg.norm(spectrum))
                        else:
                            best = 0.0
                        ratio = err/max(best,1e-12)
                        if best > 1e-10:
                            max_reconstruction_ratio = max(max_reconstruction_ratio,ratio)
                        assert err < max(3.0*best,1e-8), (m,n,cmplx,variant,k,err,best)
                    if k and a.size:
                        exact_approx_s = np.linalg.svd(approx/ max(np.max(abs(approx)),1e-300),compute_uv=False)
                        if np.max(abs(approx)):
                            singular_error = float(np.max(np.abs(exact_approx_s[:k]-s/np.max(abs(approx)))))
                            max_singular_relative_error = max(max_singular_relative_error,singular_error)
                            assert singular_error < 1e-8
                    count += 1
                    total_forward += calls[0]
                    total_adjoint += calls[1]
                    if m == 0 or n == 0 or target == 0:
                        continue
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore')
                            ur,sr,vr = interpolative.svd(aslinearoperator(a),target,
                                rng=np.random.default_rng(1115))
                    except (ValueError, RuntimeError, IndexError, np.linalg.LinAlgError,
                            TypeError, FloatingPointError) as exc:
                        scipy_failures += 1
                        key = f"{type(exc).__name__}: {str(exc)[:85]}"
                        scipy_failure_causes[key] = scipy_failure_causes.get(key,0) + 1
                    else:
                        sci_approx = (ur*sr) @ vr.conj().T
                        sci_err = relative_error(a,sci_approx)
                        if (not np.isfinite(sci_err) or check_orthonormal(ur) >= 1e-7 or
                                check_orthonormal(vr) >= 1e-7):
                            scipy_bad += 1
                        else:
                            scipy_ok += 1
                            max_scipy_residual = max(max_scipy_residual,sci_err)
    # Exact reproducibility uses the same seed and original operator.
    x = rng.normal(size=(9,12))
    left = fortran_svd(exe,x,3,seed=1234)
    right = fortran_svd(exe,x,3,seed=1234)
    assert all(np.array_equal(a,b) if isinstance(a,np.ndarray) else a == b
               for a,b in zip(left,right))
    print(f"SciFort validated cases: {count}")
    print(f"SciPy LinearOperator SVD reference calls: {scipy_ok} validated, "
          f"{scipy_failures} raised, {scipy_bad} returned invalid factors")
    print(f"SciPy failure categories: {scipy_failure_causes}")
    print(f"Maximum fixed-rank error / optimal Frobenius tail: {max_reconstruction_ratio:.6g}")
    print(f"Maximum factor orthogonality error: {max_ortho:.6g}")
    print(f"Maximum singular-value mismatch of reconstructed approximation: {max_singular_relative_error:.6g}")
    print(f"Maximum adaptive relative residual: {max_adaptive_error:.6g}")
    print(f"Maximum SciPy SVD relative residual: {max_scipy_residual:.6g}")
    print(f"Total forward/adjoint callbacks: {total_forward}/{total_adjoint}")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: compare_linalg_interpolative_svd_operator_scipy.py DRIVER')
    test(Path(sys.argv[1]).resolve())
