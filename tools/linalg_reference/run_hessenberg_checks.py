#!/usr/bin/env python3
"""Validate real/complex Hessenberg reductions against SciPy 1.17+.

Householder reductions are not unique. Comparisons therefore check the
similarity identity, unitary factors, structural zeros and invariant norms,
not an elementwise match to SciPy's particular LAPACK reflector convention.
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
SOURCES = ('src/kinds.f90', 'src/linalg.f90', 'src/linalg_hessenberg.f90')


def random_cases(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = sorted(set([0, 1, 2, 3, 4, 5, 8, 12, 16, max_order]))
    cases = []
    for idx in range(count):
        n = int(rng.choice(sizes))
        complex_case = (idx // 16) % 2 == 1
        a = rng.standard_normal((n, n))
        if complex_case:
            a = a + 1j * rng.standard_normal((n, n))
        family = idx % 16
        if family == 1:
            a[:] = 0
        elif family == 2:
            a = np.triu(a)
        elif family == 3:
            a = np.diag(np.diag(a))
        elif family == 4:
            a = (a + a.conj().T) / 2
        elif family == 5:
            a = np.outer(a[:, 0], a[-1, :]) if n else a
        elif family == 6:
            a = np.triu(a, -1)
        elif family == 7:
            a = np.triu(a) + 1e-12 * np.tril(a, -1)
        elif family == 8:
            a = a - a.conj().T
        elif family == 9:
            a = np.triu(a, 1)
        elif family == 10:
            a = np.eye(n, dtype=a.dtype) * 2.5
            if n > 2:
                a[0, -1] = 0.75
        elif family == 11:
            a = a * 1e300
        elif family == 12:
            a = a * 1e-300
        elif family == 13 and n > 1:
            # Large scale disparities that are not an algorithmic failure.
            a[0, :] *= 1e-12
        elif family == 14 and n > 1:
            # Ill-scaled but finite array with both large and small entries.
            a[0, :] *= 1e-100
        elif family == 15 and n > 1:
            # Repeated eigenvalues and a nontrivial nilpotent component.
            a = np.eye(n, dtype=a.dtype)
            a += np.diag(np.ones(n - 1, dtype=a.dtype), 1)
        cases.append(np.asarray(a))
    return cases


def compile_driver(tmp: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    exe = tmp / ('hessenberg_' + mode)
    args = ['gfortran', '-std=f2018', *flags,
            *(str(ROOT / s) for s in SOURCES),
            str(ROOT / 'tools/linalg_reference/driver_hessenberg.f90'),
            '-o', str(exe)]
    proc = subprocess.run(args, cwd=tmp, capture_output=True, text=True, check=False)
    if proc.returncode:
        raise RuntimeError('Fortran compilation failed:\n' + proc.stderr)
    return exe


def run_native(tmp: Path, mode: str) -> None:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    exe = tmp / ('hessenberg_native_' + mode)
    cmd = ['gfortran', '-std=f2018', *flags,
           *(str(ROOT / s) for s in SOURCES),
           str(ROOT / 'test/test_linalg_hessenberg.f90'), '-o', str(exe)]
    proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, check=False)
    if proc.returncode:
        raise RuntimeError('native compilation failed:\n' + proc.stderr)
    proc = subprocess.run([str(exe)], cwd=tmp, capture_output=True, text=True, check=False)
    if proc.returncode:
        raise RuntimeError('native tests failed:\n' + proc.stdout + proc.stderr)
    print(f'{mode}: {proc.stdout.strip()}')


def run_cases(exe: Path, cases, tol: float):
    lines = [str(len(cases))]
    for a in cases:
        n = a.shape[0]
        lines.append(f'{2 if np.iscomplexobj(a) else 1} {n}')
        for value in a.ravel(order='F'):
            lines.append(f'{value.real:.17e} {value.imag:.17e}')
    proc = subprocess.run([str(exe)], input='\n'.join(lines) + '\n',
                          capture_output=True, text=True, timeout=180, check=False,
                          env={**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
    if proc.returncode:
        raise RuntimeError('Fortran driver failed:\n' + proc.stderr[-3000:])
    values = iter(proc.stdout.splitlines())
    worst = {'reconstruction': 0., 'orthogonality': 0., 'frobenius': 0., 'trace': 0.}
    for idx, a in enumerate(cases):
        n = a.shape[0]
        status = int(next(values).strip())
        if status:
            raise AssertionError(f'case {idx}: status={status} n={n}')
        h = np.array([complex(*map(float, next(values).split())) for _ in range(n*n)],
                     dtype=complex).reshape((n,n), order='F')
        q = np.array([complex(*map(float, next(values).split())) for _ in range(n*n)],
                     dtype=complex).reshape((n,n), order='F')
        # Reference SciPy both checks the intended public contract and is used
        # to compare the numerical quality of the similarity reconstruction.
        ref_h, ref_q = linalg.hessenberg(a, calc_q=True)
        amax = float(np.max(abs(a))) if n else 0.
        scale = amax if amax else 1.
        aa, hh, rh = a/scale, h/scale, ref_h/scale
        denom = max(1., linalg.norm(aa, 'fro'))
        reconstruction = linalg.norm(q@hh@q.conj().T - aa, 'fro')/denom
        ref_residual = linalg.norm(ref_q@rh@ref_q.conj().T - aa,'fro')/denom
        orthogonality = linalg.norm(q.conj().T@q - np.eye(n), 'fro')/max(1,n)
        trace_error = abs(np.trace(hh) - np.trace(aa))/denom
        fro_error = abs(linalg.norm(hh,'fro') - linalg.norm(aa,'fro'))/denom
        values_to_check = {'reconstruction': reconstruction,
                           'orthogonality': orthogonality,
                           'frobenius': fro_error, 'trace': trace_error}
        for key, value in values_to_check.items():
            worst[key] = max(worst[key],float(value))
            if not np.isfinite(value) or value > tol:
                raise AssertionError(f'case {idx} n={n} {key}={value:g} reference={ref_residual:g}')
        if np.count_nonzero(np.tril(h, -2)):
            raise AssertionError(f'case {idx} has nonzero entries below subdiagonal')
        if not np.iscomplexobj(a) and (np.any(h.imag) or np.any(q.imag)):
            raise AssertionError(f'case {idx} returned complex output for real input')
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['checked','optimized','both'], default='both')
    parser.add_argument('--cases',type=int,default=250)
    parser.add_argument('--seed',type=int,default=20261107)
    parser.add_argument('--max-order',type=int,default=48)
    parser.add_argument('--tol',type=float,default=3e-12)
    args = parser.parse_args()
    if args.cases < 1 or args.max_order < 1 or args.tol <= 0:
        parser.error('positive --cases, --max-order and --tol required')
    cases = random_cases(args.cases,args.max_order,args.seed)
    modes = ['checked','optimized'] if args.mode == 'both' else [args.mode]
    with tempfile.TemporaryDirectory(prefix='scifort_hessenberg_') as path:
        tmp = Path(path)
        for mode in modes:
            run_native(tmp,mode)
            exe = compile_driver(tmp,mode)
            worst = run_cases(exe,cases,args.tol)
            print(f'{mode}: {len(cases)} SciPy comparisons passed; '
                  + ', '.join(f'{key}={value:.3e}' for key,value in worst.items()))


if __name__ == '__main__':
    main()
