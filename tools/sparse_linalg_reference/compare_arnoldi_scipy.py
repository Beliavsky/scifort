#!/usr/bin/env python3
"""Compare SciFort thick-restarted Arnoldi eigenpairs with SciPy eigenvalues.

The matrices have separated target eigenvalues, avoiding selection ambiguity
from eigenvalue ties and conjugate-pair tie-breaking in ARPACK.  The tests
check the *true* eigenpair residual in the native solver and the requested
spectral ordering, not the arbitrary eigenvector phases or order within ties.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs

ROOT = Path(__file__).resolve().parents[2]


def make_matrix(rng: np.random.Generator, n: int, mode: str) -> np.ndarray:
    a = np.diag(rng.uniform(-0.6, 0.6, n))
    rows = np.arange(n - 1)
    a[rows, rows + 1] = rng.uniform(-0.04, 0.04, n - 1)
    a[rows + 1, rows] = rng.uniform(-0.04, 0.04, n - 1)
    if mode in ('LM', 'LR'):
        a[0, 0] = rng.uniform(5, 9)
        a[1, 1] = rng.uniform(4, 4.5)
    elif mode == 'SR':
        a[0, 0] = -rng.uniform(5, 9)
        a[1, 1] = -rng.uniform(4, 4.5)
    else:
        a[:4, :4] = 0
        u, v = rng.uniform(5, 9), rng.uniform(3, 4.5)
        a[0, 1], a[1, 0] = -u, u
        a[2, 3], a[3, 2] = -v, v
    return a


def priority(values: np.ndarray, which: str) -> np.ndarray:
    return {'LM': lambda v: abs(v), 'SM': lambda v: -abs(v),
            'LR': lambda v: v.real, 'SR': lambda v: -v.real,
            'LI': lambda v: v.imag, 'SI': lambda v: -v.imag}[which](values)


def run_one(exe: Path, a: np.ndarray, k: int, which: str, ncv: int,
            cycles: int = 120) -> tuple[np.ndarray, np.ndarray, int, int]:
    n = a.shape[0]
    text = f'{n} {k} {ncv} {cycles} {which} 1e-8\n'
    text += ' '.join(repr(float(v)) for v in a.flatten(order='F')) + '\n'
    result = subprocess.run([str(exe)], input=text, capture_output=True,
                            text=True, check=True, timeout=35)
    lines = result.stdout.splitlines()
    if len(lines) != k + 1:
        raise AssertionError(f'bad output: {result.stdout[:250]}')
    status, steps, locked, matvecs = map(int, lines[0].split())
    vals = np.array([complex(float(fields[0]), float(fields[1]))
                     for line in lines[1:] for fields in [line.split()]])
    errs = np.array([float(line.split()[2]) for line in lines[1:]])
    if status != 0:
        raise AssertionError(f'nonconvergence: {status=} {which=} {steps=}')
    if locked != k or max(errs) > 3e-6:
        raise AssertionError(f'incorrect true residual: {which=} {errs=} {locked=}')
    if matvecs < ncv:
        raise AssertionError(f'matvec counter too small: {matvecs}')
    return vals, errs, steps, matvecs


def test_case(exe: Path, a: np.ndarray, k: int, mode: str, p: int) -> int:
    vals, errors, steps, matvecs = run_one(exe, a, k, mode, p)
    dense = eigvals(a)
    targets = dense[np.argsort(-priority(dense, mode))[:k]]
    rows, cols = linear_sum_assignment(abs(vals[:, None] - targets[None, :]))
    delta = max(abs(vals[rows] - targets[cols]))
    if delta > 1e-4:
        raise AssertionError(f'eigenvalue mismatch {mode=} {delta=} {targets=} {vals=}')
    # The same well-separated eigenvalues should be identified by ARPACK.
    if mode in ('LM', 'LR', 'SR') and p < a.shape[0]:
        scipy_vals = eigs(csr_matrix(a), k=k, which=mode, tol=1e-11,
                          ncv=min(a.shape[0], max(2 * k + 2, p)),
                          v0=np.sin(np.arange(1, a.shape[0] + 1, dtype=float)),
                          maxiter=1500, return_eigenvectors=False)
        rows, cols = linear_sum_assignment(abs(vals[:, None] - scipy_vals[None, :]))
        if max(abs(vals[rows] - scipy_vals[cols])) > 1e-4:
            raise AssertionError(f'ARPACK eigenvalue mismatch {mode=}')
    return steps


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    parser.add_argument('--cases', type=int, default=120)
    parser.add_argument('--seed', type=int, default=20261104)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error('--cases must be positive')
    modes = ('LM', 'LR', 'SR', 'LI', 'SI')
    for build_mode in (('checked', 'optimized') if args.mode == 'both' else (args.mode,)):
        rng = np.random.default_rng(args.seed)
        with tempfile.TemporaryDirectory(prefix='scifort-arnoldi-') as temp:
            build = Path(temp)
            flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if build_mode == 'checked' else ['-O3']
            cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
            objects = []
            for name, rel in [('kinds', 'src/kinds.f90'),
                              ('sparse', 'src/sparse.f90'),
                              ('lsmr', 'src/optimize/lsmr.f90'),
                              ('sparse_linalg', 'src/sparse_linalg.f90')]:
                obj = build / (name + '.o')
                subprocess.run([*cc, '-c', str(ROOT / rel), '-o', str(obj)],
                               cwd=build, check=True)
                objects.append(obj)
            exe = build / 'arnoldi_reference'
            subprocess.run([*cc, str(ROOT / 'tools/sparse_linalg_reference'
                                      '/sparse_arnoldi_reference_driver.f90'),
                            *(str(o) for o in objects), '-o', str(exe)],
                           cwd=build, check=True)
            cycles = []
            for i in range(args.cases):
                mode = modes[i % len(modes)]
                n = (24, 48, 64)[(i // len(modes)) % 3]
                p = 10 if n < 64 else 12
                a = make_matrix(rng, n, mode)
                cycles.append(test_case(exe, a, 2, mode, p))
            # Clustered real eigenvalues with symmetric/nonsymmetric similarity
            # transformations; a small ncv requires several restart cycles.
            for symmetric in (True, False):
                for gap in (0.01, 1e-5):
                    n = 60
                    diagonal = rng.uniform(-4, 4, n)
                    diagonal[:3] = (8, 8 - gap, 8 - 2 * gap)
                    if symmetric:
                        transform, _ = np.linalg.qr(rng.standard_normal((n, n)))
                        a = (transform * diagonal) @ transform.T
                    else:
                        transform = np.eye(n) + 0.02 * rng.standard_normal((n, n))
                        a = (transform * diagonal) @ np.linalg.inv(transform)
                    cycles.append(test_case(exe, a, 3, 'LR', 8))
            print(f'{build_mode}: {args.cases} randomized comparisons, '
                  f'4 clustered-spectrum cases PASSED; '
                  f'{sum(c > 1 for c in cycles)} restarted cases; '
                  f'maximum {max(cycles)} cycles', flush=True)


if __name__ == '__main__':
    main()
