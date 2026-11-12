#!/usr/bin/env python3
"""Cross-check deterministic SciFort interpolative SVD against SciPy and NumPy.

Usage: python test/compare_linalg_interpolative_svd_scipy.py ./id_svd_driver

No SciPy code is copied into Fortran: the SciPy 1.17 implementation is an
independent deterministic reference (rand=False). Different valid IDs can
choose different skeletons and phases. Check reconstruction, singular values,
and orthonormality instead of comparing U and V entry by entry.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
import scipy
from scipy.linalg import interpolative as sli


def make_cases() -> list[tuple[np.ndarray, int, int | float]]:
    rng = np.random.default_rng(270919)
    cases: list[tuple[np.ndarray, int, int | float]] = []
    for m in range(0, 14):
        for n in range(0, 15):
            for dtype in (np.float64, np.complex128):
                for family in range(4):
                    a = rng.standard_normal((m, n)).astype(dtype)
                    if dtype == np.complex128:
                        a += 1j * rng.standard_normal((m, n))
                    if family == 1:
                        rank = max(1, min(m, n) // 3)
                        b = rng.standard_normal((m, rank))
                        c = rng.standard_normal((rank, n))
                        if dtype == np.complex128:
                            b = b + 1j * rng.standard_normal(b.shape)
                            c = c + 1j * rng.standard_normal(c.shape)
                        a = b @ c
                    elif family == 2:
                        a *= 1e-35
                    elif family == 3:
                        a *= 1e35
                    k = min(m, n)
                    if family == 0:
                        cases.append((a, 0, 0))
                    elif family == 1:
                        cases.append((a, 0, max(1,k//2) if k else 0))
                    elif family == 2:
                        cases.append((a, 0, k))
                    else:
                        cases.append((a, 1, 0.08))
    # Multiple spectral gaps: adaptive ranks do not always match SciPy's.
    for i in range(240):
        m, n = int(rng.integers(2, 17)), int(rng.integers(2, 19))
        dtype = (np.float64, np.complex128)[i % 2]
        left = rng.standard_normal((m,m))
        right = rng.standard_normal((n,n))
        if dtype == np.complex128:
            left = left + 1j*rng.standard_normal((m,m))
            right = right + 1j*rng.standard_normal((n,n))
        q1, _ = np.linalg.qr(left)
        q2, _ = np.linalg.qr(right)
        singular = np.geomspace(1.0, 1e-9, min(m,n))
        a = (q1[:,:len(singular)] * singular) @ q2[:,:len(singular)].conj().T
        if i % 3 == 0:
            cases.append((a,0,int(rng.integers(1,min(m,n)+1))))
        else:
            cases.append((a,1,[1e-1,1e-3,1e-6][i % 3]))
    # Rank-deficient and all-zero matrices including multiple precision choices.
    for m, n in [(0,8),(8,0),(0,0),(4,7),(7,4),(6,6)]:
        for dtype in (np.float64,np.complex128):
            a = np.zeros((m,n),dtype=dtype)
            cases.append((a,0,min(m,n)))
            cases.append((a,1,1e-4))
    return cases


def check(executable: str) -> None:
    cases = make_cases()
    lines = [str(len(cases))]
    for a, mode, spec in cases:
        m,n = a.shape
        dtype = int(np.iscomplexobj(a))
        rank = int(spec) if mode == 0 else 0
        eps = float(spec) if mode == 1 else 0.0
        lines.append(f'{dtype} {m} {n} {mode} {rank} {eps:.17g}')
        for x in a.flat:
            lines.append(f'{x.real:.17e} {x.imag:.17e}')
    result = subprocess.run([executable], input='\n'.join(lines)+'\n',
                            text=True, capture_output=True, check=True, timeout=180)
    output = iter(result.stdout.splitlines())
    max_relative_residual = 0.0
    max_abs_orthogonality = 0.0
    max_reference_s_diff = 0.0
    fixed_compared = 0
    adaptive_compared = 0
    degenerate_compared = 0
    for case_no,(a,mode,spec) in enumerate(cases,1):
        header = next(output).split()
        if len(header) != 4 or header[:2] != ['CASE',str(case_no)]:
            raise AssertionError((case_no,'unexpected header',header))
        status, k = map(int,header[2:])
        if status != 0:
            raise AssertionError((case_no,'unexpected Fortran status',status))
        m,n = a.shape
        if k < 0 or k > min(m,n):
            raise AssertionError((case_no,'invalid selected rank',k))
        if mode == 0 and k != spec:
            raise AssertionError((case_no,'wrong fixed rank',k,spec))
        s = np.fromstring(next(output),sep=' ')
        if s.shape != (k,):
            raise AssertionError((case_no,'singular value count',s.shape,k))
        def read_factor(rows: int) -> np.ndarray:
            factor = np.empty((rows,k),dtype=a.dtype)
            for i in range(rows):
                entry = np.fromstring(next(output),sep=' ')
                if np.iscomplexobj(a):
                    if entry.size != 2*k:
                        raise AssertionError((case_no,'complex factor row',entry.size,k))
                    factor[i,:] = entry[::2]+1j*entry[1::2]
                else:
                    if entry.size != k:
                        raise AssertionError((case_no,'real factor row',entry.size,k))
                    factor[i,:] = entry
            return factor
        u, v = read_factor(m), read_factor(n)
        if not (np.isfinite(s).all() and np.isfinite(u).all() and np.isfinite(v).all()):
            raise AssertionError((case_no,'nonfinite result'))
        scale = max(float(np.linalg.norm(a)),1e-90)
        if k:
            if np.any(s < -1e-13*scale) or np.any(np.diff(s) > 1e-9*scale):
                raise AssertionError((case_no,'singular values unordered'))
            orth = max(np.max(np.abs(u.conj().T@u-np.eye(k))),
                       np.max(np.abs(v.conj().T@v-np.eye(k))))
            max_abs_orthogonality = max(max_abs_orthogonality,orth)
            if orth > 2e-8:
                raise AssertionError((case_no,'factor not orthonormal',orth))
        approx = (u*s) @ v.conj().T
        rel = float(np.linalg.norm(a-approx)/scale)
        max_relative_residual = max(max_relative_residual,rel)
        if mode == 1:
            adaptive_compared += 1
            if rel > max(1e-9, float(spec)*1.03):
                raise AssertionError((case_no,'tolerance not met',rel,spec,k))
        elif k == 0 or min(m,n) == 0:
            degenerate_compared += 1
            if np.any(approx):
                raise AssertionError((case_no,'zero rank returned nonzero product'))
        else:
            fixed_compared += 1
            # Direct deterministic SciPy ID-based SVD comparison when eligible.
            try:
                us, ss, vs = sli.svd(a,int(k),rand=False)
            except (ValueError,ZeroDivisionError,RuntimeError) as exc:
                # Report, do not assert an invalid reference on degenerate inputs.
                if np.linalg.norm(a) > 0:
                    raise AssertionError((case_no,'SciPy reference failed',repr(exc))) from exc
                continue
            sci_approx = (us*ss)@vs.conj().T
            sci_res = float(np.linalg.norm(a-sci_approx)/scale)
            # Relative approximation quality, permitting different pivot choices.
            if rel > max(1e-9,12*sci_res):
                raise AssertionError((case_no,'worse than SciPy ID',rel,sci_res))
            # Weyl/Hoffman-Wielandt: values for two approximations can differ
            # by at most the sum of their matrix reconstruction errors.
            s_err = float(np.linalg.norm(s-ss))/scale
            max_reference_s_diff = max(max_reference_s_diff,s_err)
            if s_err > rel+sci_res+1e-9:
                raise AssertionError((case_no,'SVD singular values',s_err,rel,sci_res))
    if next(output,None) is not None:
        raise AssertionError('unexpected trailing output')
    print(f'SciPy {scipy.__version__}: {len(cases)} independent interpolative SVD comparisons passed')
    print(f'Fixed rank: {fixed_compared}; adaptive: {adaptive_compared}; zero rank/dimension: {degenerate_compared}')
    print(f'Max relative reconstruction error (truncated cases included): {max_relative_residual:.3e}')
    print(f'Max absolute factor orthogonality error: {max_abs_orthogonality:.3e}')
    print(f'Max scaled difference from SciPy ID singular values: {max_reference_s_diff:.3e}')


if __name__ == '__main__':
    if len(sys.argv)!=2:
        raise SystemExit('Usage: compare_linalg_interpolative_svd_scipy.py DRIVER')
    check(sys.argv[1])
