#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Stress-test repeated spectra, nonnormal coupling, and extreme scaling."""
from __future__ import annotations

import subprocess
import numpy as np
from scipy import linalg


def compare(executable: str) -> None:
    cases = []
    # Taylor blocks of various sizes, with arbitrary triangular couplings.
    for n in (2, 3, 5, 8, 12, 20, 40):
        for complex_input in (False, True):
            lam = 0.3+0.4j if complex_input else 0.3
            a = np.eye(n,dtype=complex)*lam
            for j in range(n-1):
                a[j,j+1] = 0.45+0.15j if complex_input else 0.45
            for j in range(n-2):
                a[j,j+2] = -0.13+0.07j if complex_input else -0.13
            if not complex_input:
                a = a.real
            cases.append((a, 1, 1, linalg.expm(a), 0, 'Jordan exp'))
            cases.append((a, 4, 1, a@a, 0, 'Jordan square'))
            if n < 9:
                cases.append((a, 2, 1, linalg.logm(a), 0, 'Jordan log'))
            # Without derivatives a nontrivial repeated Jordan block must fail.
            cases.append((a, 1, 0, None, 4, 'Jordan no derivatives'))
    # Two distinct repeated eigenvalue blocks, with cross-coupling.
    a = np.array([[0.6, 0.9, 0.1, 0.4], [0.0, 0.6, -0.3, 0.6],
                  [0.0, 0.0, 1.8, -0.7], [0.0, 0.0, 0.0, 1.8]])
    cases.append((a, 1, 1, linalg.expm(a), 0, 'two Jordan blocks'))
    cases.append((a.astype(complex)*(1.0+0.2j), 1, 1,
                  linalg.expm(a*(1.0+0.2j)), 0, 'complex Jordan blocks'))
    # Noncontiguous repeated eigenvalues now supported via Schur reordering.
    a = np.array([[1.0,1.0,0.1], [0.0,2.0,1.0], [0.0,0.0,1.0]])
    cases.append((a, 1, 1, linalg.expm(a), 0, 'nonadjacent repeated eigenvalues'))
    cases.append((a, 1, 0, None, 4, 'nonadjacent without derivatives'))
    # Very small/very large inputs: scaling avoids overflowing the Schur QR.
    for n in (2,3,6):
        for scale in (1e-200,1e-100,1e100,1e200):
            eig = np.arange(1,n+1,dtype=float)*0.7+0.3
            a = np.diag(eig)
            for j in range(n-1):
                a[j,j+1] = 0.3
            a = scale*a
            # log(scale*B) = log(scale)*I + log(B), using an independent
            # moderate-norm SciPy logarithm rather than a huge/small input.
            expected = linalg.logm(a/scale)+np.log(scale)*np.eye(n)
            cases.append((a, 2, 0, expected, 0, 'scaled logarithm'))
            cases.append((a, 1, 0, linalg.expm(a) if scale<1 else None,
                          0 if scale<1 else 4, 'scaled exponential'))
    # Larger nonnormal triangular matrices with separated repeated diagonals.
    rng = np.random.default_rng(20261107)
    for n in (8,16,32,40):
        eigen = np.linspace(0.2,2.5,n)
        a = np.diag(eigen)+0.06*np.triu(rng.normal(size=(n,n)),1)
        cases.append((a, 1, 0, linalg.expm(a), 0, 'large triangular exp'))

    lines = [str(len(cases))]
    for a, fid, deriv, _, _, _ in cases:
        n = len(a)
        lines.append(f'{2 if np.iscomplexobj(a) else 1} {n} {fid} {deriv}')
        for j in range(n):
            for i in range(n):
                z = complex(a[i,j])
                lines.append(f'{z.real:.17e} {z.imag:.17e}')
    proc = subprocess.run([executable], input='\n'.join(lines)+'\n', text=True,
                          capture_output=True, check=True)
    outputs = proc.stdout.splitlines()
    pos = 0
    failure = []
    max_error = 0.0
    for index,(a, fid, deriv, expected, status_expected, name) in enumerate(cases):
        n = len(a)
        code, residual = outputs[pos].split()
        pos += 1
        r = np.zeros((n,n),complex)
        for j in range(n):
            for i in range(n):
                re,im = map(float, outputs[pos].split())
                pos += 1
                r[i,j] = re+1j*im
        if int(code) != status_expected:
            failure.append(f'{index}: {name}, n={n}, expected status {status_expected}, got {code}')
            continue
        if expected is not None and status_expected == 0:
            err = np.linalg.norm(r-expected)/max(1.0,np.linalg.norm(expected))
            max_error = max(max_error,err)
            if not np.isfinite(err) or err > 2e-8:
                failure.append(f'{index}: {name}, n={n}, relative err {err:.3e}')
    print(f'Stress cases: {len(cases)-len(failure)}/{len(cases)} passed; '
          f'max relative error {max_error:.3e}',flush=True)
    if failure:
        print('\n'.join(failure[:20]),flush=True)
        raise AssertionError(f'{len(failure)} stress failures')


if __name__=='__main__':
    import sys
    compare(sys.argv[1])
