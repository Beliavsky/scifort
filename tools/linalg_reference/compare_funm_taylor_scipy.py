#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check opt-in adaptive clustered-block Taylor evaluation against SciPy.

Only specialized SciPy matrix functions, rather than scipy.linalg.funm, are
used as references because SciPy's generic funm is unreliable for tightly
clustered, nonnormal eigenvalues.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy import linalg


def reference(a: np.ndarray, function_id: int) -> np.ndarray:
    return {
        1: linalg.expm,
        2: linalg.logm,
        3: linalg.sinm,
        4: lambda z: z @ z,
        6: linalg.cosm,
        9: linalg.sinhm,
    }[function_id](a)


def compare(executable: str, cases: int = 240, seed: int = 20261108) -> None:
    rng = np.random.default_rng(seed)
    inputs = []
    for k in range(cases):
        n = int(rng.integers(2, 13))
        fid = (1, 2, 3, 4, 6, 9)[k % 6]
        is_complex = bool(rng.integers(0, 2))
        similarity = bool(rng.integers(0, 2))
        # All eigenvalues are within the requested 2% clustering threshold,
        # but they are not equal: the old finite Taylor polynomial could
        # reject these blocks or return the wrong answer.
        width = float(rng.uniform(2e-4, 3e-3))
        center = 1.3 + (0.25j if is_complex else 0.0)
        spectrum = center + rng.uniform(-width, width, size=n)
        if is_complex:
            spectrum = spectrum + 1j * rng.uniform(-width, width, size=n)
        a = np.diag(spectrum).astype(complex)
        amplitude = 0.18 if n <= 6 else 0.055
        a += np.triu(rng.normal(size=(n, n)), 1) * amplitude
        if is_complex:
            a += 1j * np.triu(rng.normal(size=(n, n)), 1) * amplitude
        if similarity:
            q = rng.normal(size=(n, n))
            if is_complex:
                q = q + 1j * rng.normal(size=(n, n))
            q, _ = np.linalg.qr(q)
            a = q @ a @ q.conj().T
        if not is_complex:
            a = a.real
        inputs.append((a, fid, reference(a, fid)))

    # Deterministic larger, weakly nonnormal spectra.  The second family
    # includes two different close clusters with inter-block couplings,
    # exercising both adaptive diagonal blocks and block Sylvester solves.
    for n in (16, 24, 32, 40):
        for is_complex in (False, True):
            for grouped in (False, True):
                middle = n // 2
                spectrum = np.full(n, 1.3+0.2j if is_complex else 1.3,
                                   dtype=complex)
                spectrum += rng.uniform(-4e-4, 4e-4, size=n)
                if is_complex:
                    spectrum += 1j*rng.uniform(-4e-4, 4e-4, size=n)
                if grouped:
                    spectrum[middle:] += 0.8
                    rng.shuffle(spectrum)
                a = np.diag(spectrum)
                a += 0.012*np.triu(rng.normal(size=(n, n)), 1)
                if is_complex:
                    a += 0.012j*np.triu(rng.normal(size=(n, n)), 1)
                else:
                    a = a.real
                for fid in (1, 2, 3, 6):
                    inputs.append((a, fid, reference(a, fid)))

    total = len(inputs)
    lines = [str(total)]
    for a, fid, _ in inputs:
        n = len(a)
        lines.append(f'{2 if np.iscomplexobj(a) else 1} {n} {fid} 2')
        for j in range(n):
            for i in range(n):
                z = complex(a[i, j])
                lines.append(f'{z.real:.17e} {z.imag:.17e}')
    run = subprocess.run([str(executable)], input='\n'.join(lines) + '\n',
                         text=True, capture_output=True, check=True)
    output = run.stdout.splitlines()
    index = 0
    max_error = 0.0
    failures = []
    for case, (a, fid, expected) in enumerate(inputs):
        n = len(a)
        status, residual = output[index].split()
        index += 1
        result = np.empty((n, n), dtype=complex)
        for j in range(n):
            for i in range(n):
                re, im = map(float, output[index].split())
                result[i, j] = re + 1j * im
                index += 1
        if int(status) != 0:
            failures.append(f'{case}: order={n}, function={fid}, status={status}')
            continue
        error = np.linalg.norm(result-expected, 'fro') / max(1.0, np.linalg.norm(expected, 'fro'))
        max_error = max(max_error, error)
        if not np.isfinite(error) or error > 4e-9 or float(residual) > 2e-8:
            failures.append(f'{case}: order={n}, function={fid}, error={error:.3e}, '
                            f'commutator={residual}')
    print(f'Adaptive Taylor SciPy comparisons: {total-len(failures)}/{total} passed '
          f'({cases} random, {total-cases} larger/mixed blocks), '
          f'max relative Frobenius error={max_error:.3e}', flush=True)
    if failures:
        print('\n'.join(failures[:25]), flush=True)
        raise AssertionError(f'{len(failures)} comparisons failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases', type=int, default=240)
    parser.add_argument('--seed', type=int, default=20261108)
    args = parser.parse_args()
    compare(args.executable, args.cases, args.seed)
