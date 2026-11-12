#!/usr/bin/env python3
"""Check SciFort matrix_balance using independent SciPy/LAPACK reference.

The balancing diagonal and permutation are generally nonunique, so compare
similarity identities, power-of-two factors, eigenvalue invariants, and the
amount of off-diagonal balancing rather than requiring identical factors.
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
SOURCES = ('src/kinds.f90', 'src/linalg.f90', 'src/linalg_balance.f90')


def random_cases(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = sorted(set([0, 1, 2, 3, 4, 5, 8, 12, 16, 24, max_order]))
    cases = []
    for idx in range(count):
        n = int(rng.choice(sizes))
        complex_case = (idx // 13) % 2 == 1
        a = rng.normal(size=(n, n))
        if complex_case:
            a = a + 1j * rng.normal(size=(n, n))
        family = idx % 20
        if family == 0:
            a.fill(0)
        elif family == 1:
            a = np.diag(np.diag(a))
        elif family == 2:
            a = np.triu(a)
        elif family == 3:
            a = np.tril(a)
        elif family == 4:
            a = np.triu(a, 1)
        elif family == 5 and n >= 2:
            a[0, :] = 0
            a[0, 0] = 2
        elif family == 6 and n >= 2:
            a[:, -1] = 0
            a[-1, -1] = -3
        elif family == 7 and n >= 4:
            a[0, :] = 0
            a[:, -1] = 0
            a[0, -1] = 1
        elif family == 8:
            a *= 1e250
        elif family == 9:
            a *= 1e-250
        elif family == 10:
            a = (a + a.conj().T) / 2
        elif family == 11:
            a = a * 2.0**30
        elif family == 12:
            a = a * 2.0**-30
        elif family == 13 and n:
            a[0, :] *= 2.0**-40
        elif family == 14 and n:
            a[:, 0] *= 2.0**-30
        elif family == 15 and n > 1:
            # Explicitly ill-scaled similarity A = D C D^{-1}.
            exponents = rng.integers(-35, 36, size=n)
            d = np.exp2(exponents.astype(float))
            a = (a * d[:, None]) / d[None, :]
        elif family == 16 and n > 1:
            a = np.triu(a) + 1e-12 * np.tril(a, -1)
        elif family == 17 and n > 2:
            a[1:-1, :] = 0
            a[1:-1, 1:-1] = np.eye(n-2)
        elif family == 18 and n > 1:
            a = a * np.exp2(np.arange(n, dtype=float)[None, :] * 2)
        elif family == 19 and n > 1:
            a = a * np.exp2(np.arange(n, dtype=float)[:, None] * 2)
        permute = idx % 7 != 0
        scale = idx % 9 != 0
        cases.append((np.asarray(a), permute, scale))
    return cases


def compile_program(tmp: Path, mode: str, source: str, name: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace',
             '-ffpe-trap=invalid,zero,overflow'] if mode == 'checked' else ['-O3']
    exe = tmp / (name + '_' + mode)
    cmd = ['gfortran', '-std=f2018', *flags,
           *(str(ROOT / s) for s in SOURCES), str(ROOT / source), '-o', str(exe)]
    proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, check=False)
    if proc.returncode:
        raise RuntimeError(f'{name} compilation failed:\n' + proc.stderr)
    return exe


def run_native(tmp: Path, mode: str):
    exe = compile_program(tmp, mode, 'test/test_linalg_balance.f90', 'balance_native')
    proc = subprocess.run([str(exe)], cwd=tmp, capture_output=True, text=True, check=False)
    if proc.returncode:
        raise RuntimeError(f'native tests failed ({mode}):\n{proc.stdout}\n{proc.stderr}')
    print(mode + ': ' + proc.stdout.strip())


def run_cases(exe: Path, cases, tol: float):
    text = [str(len(cases))]
    for a, permute, scale in cases:
        n = len(a)
        text.append(f'{2 if np.iscomplexobj(a) else 1} {n} {int(permute)} {int(scale)}')
        for z in a.ravel(order='F'):
            text.append(f'{z.real:.17e} {z.imag:.17e}')
    proc = subprocess.run([str(exe)], input='\n'.join(text)+'\n',
                          capture_output=True, text=True, timeout=160, check=False,
                          env={**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
    if proc.returncode:
        raise RuntimeError('Fortran driver error:\n' + proc.stdout[-500:] + proc.stderr[-2500:])
    results = iter(proc.stdout.splitlines())
    worst = dict(similarity=0., sciPy_similarity=0., trace=0., offdiagonal_ratio=0.)
    improved = 0
    for k, (a, permute, use_scale) in enumerate(cases):
        n = len(a)
        status = int(next(results))
        if status:
            raise AssertionError(f'case {k} n={n}: SciFort status {status}')
        b = np.array([complex(*map(float, next(results).split())) for _ in range(n*n)],
                     dtype=np.complex128).reshape((n, n), order='F')
        pairs = [next(results).split() for _ in range(n)]
        factors = np.array([float(row[0]) for row in pairs])
        perm = np.array([int(row[1]) - 1 for row in pairs], dtype=int)
        if not np.array_equal(np.sort(perm), np.arange(n)):
            raise AssertionError(f'case {k}: invalid permutation {perm}')
        if np.any(factors <= 0) or np.any(~np.isfinite(factors)):
            raise AssertionError(f'case {k}: invalid scaling factors {factors}')
        if n and not np.all(np.frexp(factors)[0] == 0.5):
            raise AssertionError(f'case {k}: scale factors are not powers of two')
        if not use_scale and np.any(factors != 1):
            raise AssertionError(f'case {k}: scale=False changed factors')
        if not permute and not np.array_equal(perm, np.arange(n)):
            raise AssertionError(f'case {k}: permute=False changed ordering')
        max_abs = np.max(np.abs(a)) if n else 0.
        max_abs = max_abs if max_abs else 1.
        a_unit = a / max_abs
        b_unit = b / max_abs
        expected = a_unit[np.ix_(perm, perm)] / factors[:, None] * factors[None, :]
        denom = max(1., linalg.norm(expected))
        similarity = linalg.norm(b_unit - expected) / denom
        trace = abs(np.trace(b_unit) - np.trace(a_unit)) / max(1., linalg.norm(a_unit))
        if not np.isfinite(similarity) or similarity > tol:
            raise AssertionError(f'case {k} n={n}: similarity={similarity:.3e}')
        if not np.isfinite(trace) or trace > tol:
            raise AssertionError(f'case {k} n={n}: trace={trace:.3e}')
        if not np.iscomplexobj(a) and np.any(b.imag):
            raise AssertionError(f'case {k}: real input returned imaginary matrix')
        if not permute and not use_scale and linalg.norm(b_unit-a_unit) > tol:
            raise AssertionError(f'case {k}: no-op flags altered matrix')
        # Independent SciPy/LAPACK reference, with a possibly different valid T.
        scipy_b, (scipy_factors, scipy_perm) = linalg.matrix_balance(
            a, permute=permute, scale=use_scale, separate=True)
        s_perm = np.asarray(scipy_perm, dtype=int)
        ref_unit = scipy_b / max_abs
        reference = a_unit[np.ix_(s_perm, s_perm)] / scipy_factors[:,None] * scipy_factors[None,:]
        reference_error = linalg.norm(ref_unit - reference)/max(1., linalg.norm(reference))
        if reference_error > 10*tol:
            raise AssertionError(f'case {k}: SciPy reference discrepancy {reference_error:.3e}')
        # A balanced matrix may be nonunique. Compare the 1-norm of the
        # off-diagonal entries as a useful quality diagnostic, not a promise
        # that different permutation strategies will give identical results.
        if n > 1 and max_abs > 0:
            base = a_unit.copy()
            np.fill_diagonal(base, 0)
            test = b_unit.copy()
            np.fill_diagonal(test, 0)
            reference_b = ref_unit.copy()
            np.fill_diagonal(reference_b, 0)
            quality = np.sum(np.abs(test.real) + np.abs(test.imag))
            reference_quality = np.sum(np.abs(reference_b.real) + np.abs(reference_b.imag))
            base_quality = np.sum(np.abs(base.real) + np.abs(base.imag))
            if quality < base_quality * (1. - 1e-12):
                improved += 1
            if reference_quality > 0:
                worst['offdiagonal_ratio'] = max(worst['offdiagonal_ratio'],
                                                  quality/reference_quality)
        worst['similarity'] = max(worst['similarity'], float(similarity))
        worst['sciPy_similarity'] = max(worst['sciPy_similarity'], float(reference_error))
        worst['trace'] = max(worst['trace'], float(trace))
    return worst, improved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    parser.add_argument('--cases', type=int, default=300)
    parser.add_argument('--max-order', type=int, default=48)
    parser.add_argument('--seed', type=int, default=20261107)
    parser.add_argument('--tol', type=float, default=2e-12)
    args = parser.parse_args()
    if args.cases < 1 or args.max_order < 1 or args.tol <= 0:
        parser.error('positive --cases, --max-order and --tol required')
    cases = random_cases(args.cases,args.max_order,args.seed)
    modes = ('checked','optimized') if args.mode == 'both' else (args.mode,)
    with tempfile.TemporaryDirectory(prefix='scifort_balance_') as path:
        tmp = Path(path)
        for mode in modes:
            run_native(tmp,mode)
            exe = compile_program(tmp,mode,'tools/linalg_reference/driver_balance.f90', 'balance')
            worst, improved = run_cases(exe,cases,args.tol)
            print(f'{mode}: {len(cases)} SciPy comparisons passed, {improved} improved '
                  'off-diagonal norms; ' + ', '.join(f'{k}={v:.3e}' for k,v in worst.items()))


if __name__ == '__main__':
    main()
