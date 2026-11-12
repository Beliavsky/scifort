#!/usr/bin/env python3
"""Independent SciPy comparisons for SciFort dense symmetric/Hermitian eigh.

Compile native checked and optimized Fortran drivers, then compare ascending
standard eigenvalues and eigenpair residuals without assuming matching bases
when eigenvalues repeat. The unused input triangle is deliberately poisoned in
some cases to exercise lower/upper semantics.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    'src/kinds.f90',
    'src/linalg.f90',
    'src/linalg_banded.f90',
    'src/linalg_hessenberg.f90',
    'src/linalg_banded_tridiagonal.f90',
    'src/linalg_banded_eigen.f90',
    'src/linalg_dense.f90',
    'src/linalg_dense_eigen.f90',
]


def compile_fortran(work: Path, name: str, mode: str) -> Path:
    flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace',
              '-ffpe-trap=invalid,zero,overflow'] if mode == 'checked' else ['-O3'])
    path = work / (name + '_' + mode)
    source = ('test/test_linalg_dense_eigen.f90' if name == 'native'
              else 'tools/linalg_reference/driver_dense_eigen.f90')
    command = ['gfortran', '-std=f2018', *flags,
               *(str(ROOT / s) for s in SOURCES), str(ROOT / source), '-o', str(path)]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError('Fortran build failed:\n' + result.stderr)
    return path


def test_matrices(cases: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = np.array(sorted(set([1, 2, 3, 4, 5, 7, 10, 16, 24, 32, max_order])))
    for k in range(cases):
        n = int(rng.choice(sizes))
        is_complex = bool(k % 2)
        lower = bool((k // 2) % 2)
        family = (k // 4) % 13
        a = rng.standard_normal((n, n))
        if is_complex:
            a = a + 1j * rng.standard_normal((n, n))
        a = (a + a.conj().T) / 2
        if family == 1:
            a = a @ a.conj().T  # positive semidefinite
        elif family == 2:
            a = np.eye(n, dtype=a.dtype)
        elif family == 3:
            a = np.zeros_like(a)
        elif family == 4:
            a = np.diag(rng.choice([-2.0, 0.0, 4.0], size=n)).astype(a.dtype)
        elif family == 5:
            a = np.diag(np.linspace(-3, 3, n)).astype(a.dtype)
            if n > 1:
                ix = np.arange(n - 1)
                a[ix, ix + 1] = 0.00000001j if is_complex else 0.00000001
                a[ix + 1, ix] = a[ix, ix + 1].conj()
        elif family == 6:
            a *= 1e150
        elif family == 7:
            a *= 1e-150
        elif family == 8:
            a *= 1e-8
            a += np.eye(n) * 8
        elif family == 9:
            # Force a well-defined rank-deficient matrix with zero eigenvalues.
            if n > 1:
                a[-1, :] = 0
                a[:, -1] = 0
        elif family == 10:
            a = np.diag(rng.standard_normal(n)).astype(a.dtype)
        elif family == 11:
            # Nearly repeated eigenvalues after a unitary similarity.
            q, _ = linalg.qr(rng.standard_normal((n, n)) +
                             (1j * rng.standard_normal((n, n)) if is_complex else 0))
            spectrum = np.ones(n)
            if n > 1:
                spectrum[n // 2:] = 1.0 + 1e-7
            a = (q * spectrum) @ q.conj().T
        elif family == 12:
            # Wide range of eigenvalues, moderate enough to resolve in binary64.
            q, _ = linalg.qr(rng.standard_normal((n, n)) +
                             (1j * rng.standard_normal((n, n)) if is_complex else 0))
            spectrum = np.geomspace(1e-8, 1.0, n)
            a = (q * spectrum) @ q.conj().T
        a = np.asarray(a, dtype=np.complex128 if is_complex else np.float64)
        # Copy only the triangle SciFort is meant to use, the other is undefined.
        raw = a.copy()
        if (k % 5) == 0:
            tri = np.triu_indices(n, 1) if lower else np.tril_indices(n, -1)
            raw[tri] = np.nan
        if (is_complex and k % 7 == 1):
            raw[np.diag_indices(n)] = [complex(v.real, float("nan"))
                                       for v in raw.diagonal()]
        want_z = (k % 3 != 0)
        first, last = 0, 0
        if k % 4 in (1, 3):
            first = int(rng.integers(1, n+1))
            last = int(rng.integers(first, n+1))
        yield is_complex, lower, a, raw, want_z, first, last


def compare(executable: Path, examples, tolerance: float):
    lines = [str(len(examples))]
    for is_complex, lower, a, raw, want_z, first, last in examples:
        n = a.shape[0]
        lines.append(f'{2 if is_complex else 1} {n} {int(lower)} {first} {last} {int(want_z)}')
        for value in raw.ravel(order='F'):
            lines.append(f'{value.real:.17e} {value.imag:.17e}')
    env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1',
           'MKL_NUM_THREADS': '1'}
    result = subprocess.run([str(executable)], input='\n'.join(lines)+'\n',
                            capture_output=True, text=True, timeout=240, env=env)
    if result.returncode:
        raise RuntimeError('Fortran numerical run failed: ' + result.stderr[-2000:])
    output = iter(result.stdout.splitlines())
    worst = {'eigenvalues': 0.0, 'residual': 0.0, 'orthogonality': 0.0}
    for k, (is_complex, lower, a, raw, want_z, first, last) in enumerate(examples):
        n = a.shape[0]
        lo, hi = (first - 1, last) if first else (0, n)
        m = hi - lo
        status = int(next(output))
        if status != 0:
            raise AssertionError(f'case {k}: n={n} family={(k//4)%13}, status={status}')
        actual = np.array([float(next(output)) for _ in range(m)])
        expected = linalg.eigvalsh(a, lower=lower, subset_by_index=[lo, hi - 1],
                                   check_finite=True)
        norm_a = max(np.max(abs(a)), np.finfo(float).tiny)
        error = float(np.max(abs(expected - actual)) / norm_a)
        worst['eigenvalues'] = max(worst['eigenvalues'], error)
        if error > tolerance:
            raise AssertionError(f'case {k}: eigenvalues {error} n={n}')
        if want_z:
            entries = []
            for _ in range(n*m):
                re, im = next(output).split()
                entries.append(complex(float(re), float(im)))
            vectors = np.array(entries).reshape((n, m), order='F')
            relative_residual = np.linalg.norm(a @ vectors - vectors * actual) / (
                max(np.linalg.norm(a), np.finfo(float).tiny) * max(1.0, np.sqrt(m)))
            orth = np.linalg.norm(vectors.conj().T @ vectors - np.eye(m)) / max(1, m)
            worst['residual'] = max(worst['residual'], float(relative_residual))
            worst['orthogonality'] = max(worst['orthogonality'], float(orth))
            if relative_residual > tolerance * 10 or orth > tolerance * 10:
                raise AssertionError(f'case {k}: residual={relative_residual}, '
                                     f'orthogonality={orth}, n={n}')
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['both', 'checked', 'optimized'], default='both')
    parser.add_argument('--cases', type=int, default=300)
    parser.add_argument('--max-order', type=int, default=64)
    parser.add_argument('--seed', type=int, default=20261108)
    parser.add_argument('--tol', type=float, default=1e-9)
    args = parser.parse_args()
    if args.cases < 1 or args.max_order < 1:
        parser.error('positive --cases and --max-order required')
    cases = list(test_matrices(args.cases, args.max_order, args.seed))
    modes = ['checked', 'optimized'] if args.mode == 'both' else [args.mode]
    with tempfile.TemporaryDirectory(prefix='scifort_dense_eigen_') as temp:
        work = Path(temp)
        for mode in modes:
            native = compile_fortran(work, 'native', mode)
            result = subprocess.run([str(native)], cwd=work, capture_output=True,
                                    text=True, timeout=120)
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
            driver = compile_fortran(work, 'driver', mode)
            errors = compare(driver, cases, args.tol)
            print(f'{mode}: {len(cases)} SciPy comparisons passed; ' +
                  ', '.join(f'{k}={v:.3e}' for k, v in errors.items()))
            print(result.stdout.strip())


if __name__ == '__main__':
    main()
