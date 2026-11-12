#!/usr/bin/env python3
"""Randomized dense/SciPy verification of SciFort thick-restarted real Lanczos."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]


def make_case(rng: np.random.Generator, i: int):
    n = (28, 40, 56, 70)[(i // 4) % 4]
    k = 2
    mode = ('LA', 'SA', 'LM', 'SM')[i % 4]
    # Moderate off-diagonal coupling, separated wanted eigenvalues, and a
    # varied background spectrum; selected matrices use a dense rotation.
    diag = rng.uniform(-1.0, 1.0, n)
    if mode == 'LA':
        diag[:2] = [8.0, 6.5]
    elif mode == 'SA':
        diag[:2] = [-8.0, -6.5]
    elif mode == 'LM':
        diag[:2] = [-9.0, 7.5]
    else:
        diag[0:2] = [-0.01, 0.03]
        diag[2:] = rng.choice([-1.0, 1.0], size=n-2) * rng.uniform(1.5, 5, n-2)
    a = np.diag(diag)
    for r in range(n - 1):
        a[r, r + 1] = a[r + 1, r] = rng.uniform(-.07, .07)
    if i % 6 == 2:
        u, _ = np.linalg.qr(rng.standard_normal((n, n)))
        a = u @ a @ u.T
    p = n if mode == 'SM' else min(n, (12, 16, 20, 24)[(i//4) % 4])
    return a, k, p, 75, mode, 2.0e-8


def difficult_case(mode):
    n = 100
    diagonal = 2.0 + np.arange(1,n+1)/30.0
    a = np.diag(diagonal) + np.diag(np.full(n-1,-.1),1) + np.diag(np.full(n-1,-.1),-1)
    return a, 2, 20, 80, mode, 2.0e-8


def expected(a, k, mode):
    spectrum = np.linalg.eigvalsh(a)
    if mode == 'LA':
        return spectrum[-k:]
    if mode == 'SA':
        return spectrum[:k]
    if mode == 'LM':
        return np.sort(spectrum[np.argsort(abs(spectrum))[-k:]])
    return np.sort(spectrum[np.argsort(abs(spectrum))[:k]])


def test_mode(build_mode, count, seed):
    with tempfile.TemporaryDirectory(prefix='scifort-lanczos-') as tmp:
        build = Path(tmp)
        flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace']
                 if build_mode == 'checked' else ['-O3'])
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
        exe = build / 'sparse_lanczos_reference'
        subprocess.run([*cc, str(ROOT / 'tools/sparse_linalg_reference'
                                      '/sparse_lanczos_reference_driver.f90'),
                        *(str(o) for o in objs), '-o', str(exe)],
                       cwd=build, check=True)
        rng = np.random.default_rng(seed)
        cases = [make_case(rng, i) for i in range(count)]
        cases.extend([difficult_case('LA'), difficult_case('SA')])
        data = [str(len(cases))]
        for a, k, p, cap, mode, tol in cases:
            data.append(f'{a.shape[0]} {k} {p} {cap} {mode} {tol:.17g}')
            data.append(' '.join(f'{x:.17g}' for x in a.flatten(order='F')))
        result = subprocess.run([str(exe)], cwd=build, text=True,
                                input='\n'.join(data)+'\n', capture_output=True, check=True)
        output = iter(result.stdout.splitlines())
        total_cycles, restarted = 0, 0
        for i, (a, k, p, cap, mode, tol) in enumerate(cases, start=1):
            assert next(output).strip() == f'CASE {i}'
            status, cycles, locks, products = map(int, next(output).split())
            vals, errors, vectors = [], [], []
            for _ in range(k):
                numbers = np.fromstring(next(output), sep=' ')
                assert numbers.size == a.shape[0] + 2
                vals.append(numbers[0])
                errors.append(numbers[1])
                vectors.append(numbers[2:])
            vals = np.asarray(vals)
            errors = np.asarray(errors)
            vectors = np.column_stack(vectors)
            if status != 0:
                raise AssertionError(f'{build_mode} case={i} mode={mode} n={len(a)} '
                                     f'ncv={p}: nonconvergence status={status}, '
                                     f'cycles={cycles}, residuals={errors}')
            np.testing.assert_allclose(vals, expected(a,k,mode), rtol=1e-6, atol=1e-6,
                                       err_msg=f'{build_mode} eigenvalues case {i} {mode}')
            true_resid = np.linalg.norm(a@vectors-vectors*vals[None,:],axis=0)
            np.testing.assert_allclose(true_resid,errors,rtol=1e-5,atol=1e-11)
            assert max(true_resid) <= tol*max(1.0,np.linalg.norm(a,np.inf))*1.01
            np.testing.assert_allclose(vectors.T@vectors,np.eye(k),atol=5e-8)
            assert locks == k and products >= p and 1 <= cycles <= cap
            # The same spectral targets should match scipy ARPACK, with near-ties
            # decided by the independent dense eigensolver instead.
            spectrum = np.linalg.eigvalsh(a)
            priorities = (np.abs(spectrum) if mode in ('LM','SM') else spectrum
                          if mode == 'LA' else -spectrum)
            priorities = np.sort(priorities)
            gap = (priorities[-k] - priorities[-k-1] if mode != 'SM'
                   else priorities[k] - priorities[k-1])
            if gap > .05 and mode != 'SM' and p < len(a):
                scipy_vals = np.sort(eigsh(csr_matrix(a),k,which=mode,
                                           tol=1e-10,return_eigenvectors=False))
                np.testing.assert_allclose(vals,scipy_vals,rtol=1e-5,atol=1e-5)
            total_cycles += cycles
            restarted += cycles > 1
        assert next(output,None) is None
        print(f'{build_mode}: {count} randomized + 2 difficult Lanczos comparisons passed; '
              f'{restarted} used multiple cycles; average cycles {total_cycles/len(cases):.2f}',
              flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['checked','optimized','both'], default='both')
    parser.add_argument('--cases', type=int, default=80)
    parser.add_argument('--seed', type=int, default=20261104)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error('--cases must be positive')
    for mode in (('checked','optimized') if args.mode == 'both' else (args.mode,)):
        test_mode(mode,args.cases,args.seed)


if __name__ == '__main__':
    main()
