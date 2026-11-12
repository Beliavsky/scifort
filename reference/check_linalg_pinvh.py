#!/usr/bin/env python3
"""SciPy reference comparisons for SciFort's Hermitian pseudoinverse.

Usage: python reference/check_linalg_pinvh.py /path/to/driver_linalg_pinvh
The driver is compiled from reference/driver_linalg_pinvh.f90 and SciFort.
The script does not use the SciFort eigenvalue algorithm to make inputs.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.linalg import pinvh as scipy_pinvh


def cases():
    rng = np.random.default_rng(20261109)
    output = []
    for n in range(0, 14):
        for is_complex in (False, True):
            dtype = np.complex128 if is_complex else np.float64
            for variant in range(10):
                for repetition in range(8):
                    if n:
                        x = rng.normal(size=(n, n))
                        if is_complex:
                            x = x + 1j * rng.normal(size=(n,n))
                        q, _ = np.linalg.qr(x)
                        signs = np.where(np.arange(n) % 2, -1, 1)
                        if variant == 0:  # well conditioned indefinite
                            eigenvalues = signs * (1 + np.arange(n) / max(1, n))
                        elif variant == 1:  # positive semidefinite, low rank
                            eigenvalues = np.zeros(n)
                            eigenvalues[:min(n, max(1, n // 2))] = np.linspace(0.5, 2, min(n, max(1, n // 2)))
                        elif variant == 2:  # negative semidefinite, low rank
                            eigenvalues = np.zeros(n)
                            eigenvalues[:min(n, max(1, n // 2))] = -np.linspace(0.5, 2, min(n, max(1, n // 2)))
                        elif variant == 3:  # repeated eigenvalues
                            eigenvalues = signs * 2.0
                        elif variant == 4:  # all zeros
                            eigenvalues = np.zeros(n)
                        elif variant == 5:  # modest condition number
                            eigenvalues = signs * 10 ** (-np.arange(n) / 3.0)
                        elif variant == 6:  # designed tolerance cutoff
                            eigenvalues = signs * 10 ** (-np.arange(n) / 2.0)
                        elif variant == 7:  # very small / very large scales
                            eigenvalues = signs * (1 + np.arange(n) / max(1,n))
                        elif variant == 8:  # diagonal input
                            eigenvalues = signs * (1 + np.arange(n) / max(1,n))
                            q = np.eye(n, dtype=dtype)
                        else:  # fully dense unstructured Hermitian matrix
                            a = 0.5 * (x + x.conj().T)
                            if repetition % 3 == 0:
                                a = a + 0.5 * np.eye(n)
                            eigenvalues = None
                        if eigenvalues is not None:
                            a = (q * eigenvalues) @ q.conj().T
                            if not is_complex:
                                a = a.real.astype(dtype)
                        if variant == 7:
                            a = a * (1.e30 if repetition % 2 == 0 else 1.e-30)
                        if variant in (0,3,5,6,8) and repetition == 0:
                            # Poison only the *ignored* triangle, using finite numbers.
                            a = a.copy()
                            a[np.triu_indices(n,1)] = (100 + 1j * 30 if is_complex else 100)
                        if variant in (0,3,5,6,8) and repetition == 1:
                            a = a.copy()
                            a[np.tril_indices(n,-1)] = (100 + 1j * 30 if is_complex else 100)
                    else:
                        a = np.zeros((0,0),dtype=dtype)
                    if variant == 6:
                        atol, rtol = 0.0, 2.e-2
                    elif variant == 7 and repetition % 2 == 1:
                        atol, rtol = 1.e-35, 1.e-7
                    elif variant == 1 or variant == 2:
                        atol, rtol = 0.0, 1.e-10
                    elif variant == 4:
                        atol, rtol = 0.0, 0.0
                    else:
                        atol, rtol = 0.0, max(n,1)*np.finfo(float).eps
                    lower = repetition != 1
                    if repetition == 0 and variant in (0,3,5,6,8):
                        lower = True
                    if repetition == 1 and variant in (0,3,5,6,8):
                        lower = False
                    # With random near-rank-deficient matrices, avoid testing
                    # eigenvalues lying on a floating-point cutoff boundary.
                    output.append((a, is_complex, lower, atol, rtol))
    return output


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('driver', type=Path, help='path to the compiled Fortran reference driver')
    args = ap.parse_args()
    all_cases = cases()
    with tempfile.TemporaryDirectory(prefix='scifort-pinvh-') as tmp:
        input_path = Path(tmp)/'input.txt'
        output_path = Path(tmp)/'output.txt'
        references = []
        with input_path.open('w',encoding='utf8') as f:
            for case_id, (a, is_complex, lower, atol, rtol) in enumerate(all_cases):
                ref, rank = scipy_pinvh(a, atol=atol, rtol=rtol, lower=lower, return_rank=True)
                references.append((ref, rank))
                print(case_id, a.shape[0], int(is_complex), int(lower),
                      format(atol,'.17g'),format(rtol,'.17g'),file=f)
                for v in a.flatten(order='F'):
                    print(format(v.real,'.17g'),format(v.imag,'.17g'), file=f)
        subprocess.run([str(args.driver.resolve()), str(input_path), str(output_path)],
                       check=True,timeout=120)
        with output_path.open('r',encoding='utf8') as f:
            worst = 0.0
            worst_case = None
            worst_penrose = 0.0
            for case_id, (a, is_complex, lower, atol, rtol) in enumerate(all_cases):
                rec = f.readline().split()
                if len(rec) != 3:
                    raise AssertionError(f'case {case_id}: malformed status {rec}')
                got_id, status, got_rank = map(int,rec)
                if got_id != case_id or status != 0:
                    raise AssertionError(f'case {case_id}: status={status}, id={got_id}')
                ref, rank = references[case_id]
                if got_rank != rank:
                    raise AssertionError(f'case {case_id}: rank={got_rank}, SciPy rank={rank}, n={a.shape[0]}, '
                                         f'complex={is_complex}, lower={lower}, atol={atol}, rtol={rtol}')
                v = np.empty((a.shape[0], a.shape[0]),dtype=np.complex128 if is_complex else np.float64)
                for j in range(a.shape[0]):
                    for i in range(a.shape[0]):
                        parts = f.readline().split()
                        if len(parts) != 2:
                            raise AssertionError(f'case {case_id}: malformed entry')
                        re, im = map(float,parts)
                        v[i,j] = re+1j*im if is_complex else re
                if v.size:
                    rel = np.linalg.norm(v-ref)/max(np.linalg.norm(ref), np.finfo(float).tiny)
                    physical = (np.tril(a) + np.tril(a,-1).conj().T if lower else
                                np.triu(a) + np.triu(a,1).conj().T)
                    worst_penrose = max(worst_penrose,
                        np.linalg.norm(physical @ v @ physical-physical)/
                        max(np.linalg.norm(physical),np.finfo(float).tiny))
                    if rel > worst:
                        worst = rel
                        worst_case = (case_id,a.shape[0],is_complex,lower,atol,rtol)
                    if rel > 3.e-7:
                        raise AssertionError(f'case {case_id}: relative error={rel:.3g} '
                                             f'n={a.shape[0]} complex={is_complex} lower={lower} rank={rank}')
            if f.read().strip():
                raise AssertionError('driver produced unexpected extra output')
    print(f'PASS: {len(all_cases)} SciPy pinvh reference cases; '
          f'max relative difference {worst:.3e} at {worst_case}; '
          f'max relative Penrose residual {worst_penrose:.3e}')


if __name__ == '__main__':
    main()
