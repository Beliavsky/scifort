#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check nonadjacent repeated Schur eigenvalues against SciPy matrix functions.

The SciPy *specialized* functions provide the references: scipy.linalg.funm
itself is not reliable for defective matrices with repeated eigenvalues.
"""
from __future__ import annotations

import argparse
import subprocess
import numpy as np
from scipy import linalg


def cases(seed: int, count: int):
    rng = np.random.default_rng(seed)
    records = []
    for i in range(count):
        n = (3, 4, 6, 9, 12, 16, 24, 32, 40)[i % 9]
        is_complex = (i % 2 == 0)
        func = (1, 2, 3, 4)[(i // 2) % 4]  # exp, log, sin, square
        # Several eigenvalues recur far apart along the Schur diagonal.
        # Random triangular couplings generate nontrivial Jordan chains via
        # intermediate distinct eigenvalues as well as direct couplings.
        distinct = [0.35, 1.2, 2.4, 0.7]
        eigen = np.array([distinct[(3*j + i//3) % 4] for j in range(n)])
        if i % 5 == 0 and n >= 3:
            eigen[-1] = eigen[0]  # endpoint of an explicitly separated pair
        if i % 7 == 0:
            eigen += np.linspace(-1.0e-13, 1.0e-13, n)
        if is_complex:
            eigen = eigen * (1.0 + 0.2j)
        mat = np.diag(eigen.astype(complex))
        upper = np.triu(rng.normal(size=(n, n)), 1).astype(complex)
        if is_complex:
            upper += 1j * np.triu(rng.normal(size=(n, n)), 1)
        mat += 0.04 * upper
        if not is_complex:
            mat = mat.real
        similarity = (i % 11 == 4 and n <= 6)
        if similarity:
            # Orthogonal/unitary similarity also probes the upstream Schur QR.
            # Defective spectra become much more sensitive to roundoff.
            q = rng.normal(size=(n, n))
            if is_complex:
                q = q+1j*rng.normal(size=(n, n))
            q, _ = np.linalg.qr(q)
            mat = q @ mat @ q.conj().T
            if not is_complex:
                mat = mat.real
        if func == 1:
            reference = linalg.expm(mat)
        elif func == 2:
            reference = linalg.logm(mat)
        elif func == 3:
            reference = linalg.sinm(mat)
        else:
            reference = mat @ mat
        records.append((mat, func, reference, similarity))
    return records


def compare(executable: str, ncases=120, seed=20261107) -> None:
    records = cases(seed, ncases)
    data = [str(len(records))]
    for mat, fid, _, _ in records:
        n = mat.shape[0]
        data.append(f'{2 if np.iscomplexobj(mat) else 1} {n} {fid} 1')
        for j in range(n):
            for i in range(n):
                value = complex(mat[i, j])
                data.append(f'{value.real:.17e} {value.imag:.17e}')
    proc = subprocess.run([executable], input='\n'.join(data)+'\n',
                          text=True, capture_output=True, check=True)
    out = proc.stdout.splitlines()
    pos = 0
    max_error = 0.0
    max_similarity_error = 0.0
    failures = []
    for ix, (mat, fid, expected, similarity) in enumerate(records):
        n = len(mat)
        code, residual = out[pos].split()
        pos += 1
        result = np.zeros((n,n),dtype=complex)
        for j in range(n):
            for i in range(n):
                re,im = map(float, out[pos].split())
                pos += 1
                result[i,j] = re+1j*im
        if int(code) != 0:
            failures.append(f'case {ix}: order {n} function {fid} status {code}')
            continue
        error = np.linalg.norm(result-expected, ord='fro')/max(1.0, np.linalg.norm(expected,ord='fro'))
        if similarity:
            max_similarity_error = max(max_similarity_error, error)
        else:
            max_error = max(max_error, error)
        # The Schur QR can perturb a nonnormal, defective similarity even
        # though the triangular test cases are accurate to much tighter bounds.
        limit = 3.0e-7 if similarity else 2.0e-8
        if not np.isfinite(error) or error > limit:
            failures.append(f'case {ix}: order {n} function {fid} '
                            f'relative err {error:.3e} (limit {limit:.1e})')
    print(f'Reordered-cluster SciPy comparisons: {len(records)-len(failures)}/{len(records)} passed; '
          f'max triangular relative error {max_error:.3e}; '
          f'max dense similarity relative error {max_similarity_error:.3e}',flush=True)
    if failures:
        print('\n'.join(failures[:20]),flush=True)
        raise AssertionError(f'{len(failures)} reordered-cluster failures')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases',type=int,default=120)
    parser.add_argument('--seed',type=int,default=20261107)
    args = parser.parse_args()
    compare(args.executable,args.cases,args.seed)
