#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent numerical checks for SciFort dense analytic matrix functions.

Uses SciPy's public funm and, for Jordan blocks, independently implemented
matrix functions (SciPy funm is unreliable on defective repeated spectra).
"""
from __future__ import annotations

import argparse
import subprocess
from typing import Callable
import numpy as np
from scipy import linalg

FUNCTIONS: dict[int, Callable] = {
    1: np.exp, 2: np.log, 3: np.sin, 4: lambda z: z*z,
    5: lambda z: 1.0/(1.0+z), 6: np.cos, 7: lambda z: np.exp(-z*z),
    8: np.sqrt, 9: np.sinh,
}


def make_case(rng, n, fid, is_complex, scenario):
    eig = rng.uniform(0.35, 1.75, size=n)
    if n > 1:
        eig += 0.24*np.arange(n)
        rng.shuffle(eig)
    if is_complex:
        eig = eig+1j*rng.uniform(-0.5,0.5,size=n)
    if scenario in (3,4):
        eig = np.full(n, 0.5+0.2j if is_complex else 0.5)
        if scenario == 4 and n > 1:
            eig = eig+np.linspace(-2e-13,2e-13,n)
    elif scenario == 5 and n > 1:
        eig[1] = eig[0]
    base = np.diag(eig)
    if scenario in (1,2,3,4):
        upper = np.triu(rng.normal(size=(n,n)), k=1)*0.25
        if is_complex:
            upper = upper+1j*np.triu(rng.normal(size=(n,n)), k=1)*0.2
        base += upper
    if scenario == 2 or scenario == 5:
        q = rng.normal(size=(n,n))
        if is_complex:
            q = q+1j*rng.normal(size=(n,n))
        q, _ = np.linalg.qr(q)
        base = q @ base @ q.conj().T
    return base


def true_reference(a, fid, scenario):
    if fid == 4:
        return a@a
    if fid == 5:
        return np.linalg.solve(np.eye(a.shape[0])+a, np.eye(a.shape[0]))
    if scenario in (3,4):
        if fid == 1:
            return linalg.expm(a)
        if fid == 2:
            return linalg.logm(a)
        if fid == 3:
            return linalg.sinm(a)
        if fid == 6:
            return linalg.cosm(a)
        if fid == 8:
            return linalg.sqrtm(a)
        if fid == 9:
            return linalg.sinhm(a)
    return linalg.funm(a, FUNCTIONS[fid], disp=False)[0]


def compare(executable, cases=300, seed=20261107, max_order=7, include_jordan=True):
    rng = np.random.default_rng(seed)
    records = []
    successes = 0
    skips = 0
    for k in range(cases):
        n = int(rng.integers(1, max_order+1))
        fid = int(rng.integers(1, 10))
        is_complex = bool(rng.integers(0,2))
        scenario = int(rng.integers(0, 6 if include_jordan else 3))
        if scenario in (3,4) and fid not in (1,2,3,4):
            fid = [1,2,3,4][k % 4]
        # SciPy does not reliably handle repeated non-diagonalizable matrices.
        # The callback supplies analytic derivatives for the supported functions.
        with_derivative = int(scenario in (3,4))
        a = make_case(rng,n,fid,is_complex,scenario)
        if not is_complex:
            a = a.real
        expected = true_reference(a, fid, scenario)
        records.append((a, fid, is_complex, scenario, with_derivative, expected))
    items = [str(len(records))]
    for a, fid, complex_input, _, deriv, _ in records:
        n = a.shape[0]
        items.append(f'{2 if complex_input else 1} {n} {fid} {deriv}')
        for j in range(n):
            for i in range(n):
                z = complex(a[i,j])
                items.append(f'{z.real:.17e} {z.imag:.17e}')
    proc = subprocess.run([str(executable)], input='\n'.join(items)+'\n',
                          text=True, capture_output=True, check=True)
    output = proc.stdout.splitlines()
    pos = 0
    max_error = 0.0
    errors = []
    for index, (a, fid, _, scenario, _, expected) in enumerate(records):
        n = len(a)
        status, residual = output[pos].split()
        pos += 1
        result = np.empty((n,n),complex)
        for j in range(n):
            for i in range(n):
                re, im = map(float, output[pos].split())
                result[i,j] = re+1j*im
                pos += 1
        if int(status) != 0:
            errors.append(f'case={index}, n={n}, fid={fid}, scenario={scenario}: status={status}')
            continue
        error = np.linalg.norm(result-expected,ord='fro')/max(1.0,np.linalg.norm(expected,ord='fro'))
        max_error = max(max_error, error)
        if not np.isfinite(error) or error > 5e-8:
            errors.append(f'case={index}, n={n}, fid={fid}, scenario={scenario}: err={error:.3e}')
        elif float(residual)>1e-7:
            errors.append(f'case={index}: commutator residual={residual}')
        else:
            successes += 1
    print(f'SciPy comparisons: {successes}/{len(records)} passed, '
          f'max relative Frobenius error={max_error:.3e}', flush=True)
    if errors:
        print('\n'.join(errors[:25]),flush=True)
        raise AssertionError(f'{len(errors)} comparisons failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases', type=int, default=300)
    parser.add_argument('--seed', type=int, default=20261107)
    parser.add_argument('--max-order', type=int, default=7)
    parser.add_argument('--no-jordan', action='store_true')
    args = parser.parse_args()
    compare(args.executable, args.cases, args.seed, args.max_order,
            not args.no_jordan)
