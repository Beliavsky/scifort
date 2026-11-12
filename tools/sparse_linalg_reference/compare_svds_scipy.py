#!/usr/bin/env python3
"""Check SciFort's partial sparse SVD against NumPy and SciPy.

Includes full-basis dense-reference problems, rank deficiency, repeated
singular values, and genuinely partial Golub--Kahan searches on tall/wide
sparse matrices. Every returned singular triplet gets a true two-sided
residual and orthogonality check, regardless of agreement in singular values.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import svds as scipy_svds

ROOT = Path(__file__).resolve().parents[2]


def small_case(rng: np.random.Generator, i: int):
    q = int(rng.integers(3, 32))
    p = q + int(rng.integers(0, 20))
    m, n = (p, q) if i % 2 else (q, p)
    if i % 19 == 0:
        a = np.zeros((m, n))
    elif i % 17 == 0:
        a = np.zeros((m, n))
        for j in range(q):
            a[j, j] = float(q - j) if j < q - 2 else 0.0
        a[0, 0] = a[1, 1] = float(q + 2)
    else:
        a = rng.standard_normal((m, n)) * (rng.random((m, n)) < .35)
        if i % 9 == 0:
            if n <= m:
                a[:, -1] = a[:, 0]
            else:
                a[-1, :] = a[0, :]
        if i % 13 == 0:
            a *= 0.01
    k = min(q, 1 + (i % 3))
    mode = 'SM' if i % 2 == 0 else 'LM'
    return a, k, q, q, mode, 2e-8, False


def partial_case(rng: np.random.Generator, i: int):
    q = (55, 85, 125, 210)[i % 4]
    p = q + 30
    m, n = (p, q) if i % 2 else (q, p)
    a = np.zeros((m, n))
    diagonal = rng.uniform(.05, 1.0, q)
    a[np.arange(q), np.arange(q)] = diagonal
    a[0, 0] = 14.0
    a[1, 1] = 9.0
    a[2, 2] = 6.0
    for _ in range(q):
        row = int(rng.integers(m))
        col = int(rng.integers(n))
        a[row, col] += rng.uniform(-.04, .04)
    return a, 2, 20, 20, 'LM', 2e-8, True


def test_one_mode(mode: str, count: int, partial: int, seed: int):
    with tempfile.TemporaryDirectory(prefix='scifort-svds-') as temp:
        build = Path(temp)
        flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace']
                 if mode == 'checked' else ['-O3'])
        cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
        objs = []
        for name, rel in [('kinds', 'src/kinds.f90'),
                          ('sparse', 'src/sparse.f90'),
                          ('lsmr', 'src/optimize/lsmr.f90'),
                          ('sparse_linalg', 'src/sparse_linalg.f90')]:
            obj = build / (name + '.o')
            subprocess.run([*cc, '-c', str(ROOT / rel), '-o', str(obj)],
                           cwd=build, check=True)
            objs.append(obj)
        exe = build / 'sparse_svds_reference'
        subprocess.run([*cc, str(ROOT / 'tools/sparse_linalg_reference'
                                      '/sparse_svds_reference_driver.f90'),
                        *(str(o) for o in objs), '-o', str(exe)],
                       cwd=build, check=True)
        rng = np.random.default_rng(seed)
        cases = [small_case(rng, i) for i in range(count)]
        cases += [partial_case(rng, i) for i in range(partial)]
        data = [str(len(cases))]
        for a, k, width, cap, which, tol, _ in cases:
            rows, cols = np.nonzero(a)
            values = a[rows, cols]
            # Add a duplicate pair and a stored zero without changing A.
            # COO construction must combine those entries canonically.
            coordinates = list(zip(rows + 1, cols + 1, values, strict=True))
            if values.size:
                row, col, value = coordinates[0]
                coordinates[0] = row, col, value / 2
                coordinates.append((row, col, value / 2))
                coordinates.append((row, col, 0.0))
            data.append(f'{a.shape[0]} {a.shape[1]} {k} {width} {cap} '
                        f'{which} {tol:.17g} {len(coordinates)}')
            data.extend(f'{r} {c} {float(v):.17g}' for r, c, v in coordinates)
        result = subprocess.run([str(exe)], cwd=build,
                                input='\n'.join(data) + '\n',
                                text=True, capture_output=True, check=True)
        output = iter(result.stdout.splitlines())
        total_steps, partial_success = 0, 0
        for case_id, (a, k, width, cap, which, tol, partial_search) in enumerate(cases, 1):
            assert next(output).strip() == f'CASE {case_id}'
            status, steps = map(int, next(output).split())
            singular, residuals, left, right = [], [], [], []
            for _ in range(k):
                record = np.fromstring(next(output), sep=' ')
                assert record.size == a.shape[0] + a.shape[1] + 2
                singular.append(record[0])
                residuals.append(record[1])
                left.append(record[2:2+a.shape[0]])
                right.append(record[2+a.shape[0]:])
            singular = np.asarray(singular)
            residuals = np.asarray(residuals)
            left = np.column_stack(left)
            right_t = np.vstack(right)
            if status != 0:
                raise AssertionError(f'{mode}: svds case {case_id}, shape={a.shape}, '
                                     f'{which}, ncv={width}: info={status}, '
                                     f'iters={steps}, residuals={residuals}')
            if partial_search:
                expected = np.sort(scipy_svds(csr_matrix(a), k=k, which='LM',
                                               tol=1e-10,
                                               return_singular_vectors=False))
            else:
                spectrum = np.linalg.svd(a, compute_uv=False)
                expected = np.sort(spectrum[:k] if which == 'LM' else spectrum[-k:])
            np.testing.assert_allclose(singular, expected, rtol=3e-6, atol=5e-8,
                                       err_msg=f'{mode}: singular values case {case_id}')
            assert np.all(singular[:-1] <= singular[1:] + 1e-12)
            true_left = a @ right_t.T - left * singular[None, :]
            true_right = a.T @ left - right_t.T * singular[None, :]
            true_residuals = np.hypot(np.linalg.norm(true_left, axis=0),
                                      np.linalg.norm(true_right, axis=0))
            np.testing.assert_allclose(true_residuals, residuals, rtol=2e-5,
                                       atol=2e-10)
            assert max(true_residuals) <= 4e-8*max(1, np.linalg.norm(a, 2))
            np.testing.assert_allclose(left.T @ left, np.eye(k), atol=2e-7)
            np.testing.assert_allclose(right_t @ right_t.T, np.eye(k), atol=2e-7)
            total_steps += steps
            partial_success += int(partial_search and steps < min(a.shape))
        assert next(output, None) is None
        print(f'{mode}: {count} full-basis and {partial} partial sparse SVD '
              f'comparisons passed; {partial_success} partial searches; '
              f'average {total_steps/len(cases):.1f} bidiagonalization steps',
              flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'),
                        default='both')
    parser.add_argument('--cases', type=int, default=160)
    parser.add_argument('--partial', type=int, default=60)
    parser.add_argument('--seed', type=int, default=20261104)
    args = parser.parse_args()
    if args.cases < 0 or args.partial < 0:
        parser.error('case counts must be nonnegative')
    modes = ('checked', 'optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        test_one_mode(mode, args.cases, args.partial, args.seed)


if __name__ == '__main__':
    main()
