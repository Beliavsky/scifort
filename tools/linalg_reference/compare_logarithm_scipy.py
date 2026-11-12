#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent SciPy comparisons for small dense matrix logs and powers.

Calls the Fortran reference executable with a batch of matrices over stdin.
This script is intentionally a *test*; no SciPy code is compiled into SciFort.
"""
import argparse
import os
import subprocess

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from scipy.linalg import expm, fractional_matrix_power, logm


def cases_generator(count: int, seed: int, max_order: int = 5):
    rng = np.random.default_rng(seed)
    for k in range(count):
        n = int(rng.integers(1, max_order + 1))
        group = k % 8
        real = group in (0, 1, 2, 4, 6)
        power = (k // 8) % 2 == 1
        t = float(rng.choice([-2.0, -1.0, -0.5, 0.2, 0.5, 1.5, 2.0, 3.0]))
        if group == 0:
            z = rng.normal(size=(n, n))
            q, _ = np.linalg.qr(z)
            eigen = np.exp(rng.uniform(-2.0, 2.0, size=n))
            a = (q * eigen) @ q.T
        elif group == 1:
            z = rng.normal(size=(n, n))
            a = expm(z * 0.5)
        elif group == 2:
            # Nonnormal upper triangular with clustered (occasionally repeated) eigenvalues.
            a = np.eye(n) + np.triu(rng.uniform(-1.0, 1.0, (n, n)), 1)
        elif group == 3:
            z = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
            a = expm(z * 0.4)
        elif group == 4:
            # A real matrix with a negative eigenvalue, an otherwise benign spectrum.
            a = np.diag(rng.uniform(0.5, 2.0, size=n))
            a[0, 0] *= -1.0
        elif group == 5:
            # Complex diagonalizable matrix with negative-real-side eigenvalues.
            z = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
            q, _ = np.linalg.qr(z)
            eigen = np.exp(rng.uniform(-1.5, 1.5, size=n) +
                           1j*rng.uniform(-3.0, 3.0, size=n))
            a = (q * eigen) @ q.conj().T
        elif group == 6:
            scales = [-90, -60, 60, 90] if power else [-150, -60, 60, 150]
            a = np.eye(n) * 10.0**float(rng.choice(scales))
        else:
            z = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
            a = np.eye(n) + z*1e-3
        a = np.asarray(a, dtype=np.complex128 if not real else np.float64)
        mode = (3 if real else 4) if power else (1 if real else 2)
        yield mode, a, t, group


def compare(executable: str, cases: int, seed: int, verbose: bool = False,
            max_order: int = 5):
    entries = list(cases_generator(cases, seed, max_order))
    lines = [str(cases)]
    for mode, a, t, group in entries:
        n = len(a)
        lines.append(f'{mode} {n} {t:.17g}')
        for j in range(n):
            for i in range(n):
                z = a[i, j]
                lines.append(f'{np.real(z):.17e} {np.imag(z):.17e}')
    proc = subprocess.run([executable], input='\n'.join(lines)+'\n', text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    output = iter(proc.stdout.splitlines())
    errors = []
    bad = []
    failures = []
    for index, (mode, a, t, group) in enumerate(entries):
        n = len(a)
        status_line = next(output).split()
        status = int(status_line[0])
        residual = float(status_line[1])
        x = np.zeros((n, n), dtype=complex)
        for j in range(n):
            for i in range(n):
                re, im = next(output).split()
                x[i, j] = float(re) + 1j * float(im)
        if status != 0:
            failures.append((index, group, mode, n, t, status))
            continue
        ref = np.asarray(logm(a) if mode in (1, 2) else fractional_matrix_power(a, t),
                         dtype=complex)
        diff = np.linalg.norm(x-ref, ord=1)/max(1.0, np.linalg.norm(ref, ord=1))
        errors.append(diff)
        if diff > 2e-7:
            bad.append((index, group, mode, n, t, float(diff), residual))
    print(f'count={cases}, success={len(errors)}, failures={len(failures)}, '
          f'>2e-7 errors={len(bad)}, worst={max(errors, default=0):.3g}, '
          f'median={np.median(errors) if errors else 0:.3g}')
    if verbose or bad or failures:
        print('largest mismatches:', sorted(bad, key=lambda b: -b[-2])[:12])
        print('first failures:', failures[:18])
    if failures or bad:
        raise AssertionError('some dense matrix logarithm/power comparisons did not pass')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('--cases', type=int, default=240)
    parser.add_argument('--seed', type=int, default=20261107)
    parser.add_argument('--max-order', type=int, default=5)
    args = parser.parse_args()
    compare(args.executable, args.cases, args.seed, True, args.max_order)


if __name__ == '__main__':
    main()
