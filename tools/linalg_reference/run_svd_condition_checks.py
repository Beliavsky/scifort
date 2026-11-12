#!/usr/bin/env python3
"""Compare SciFort singular-value-only Jacobi, cond2 and matrix_rank to SciPy/NumPy.

The numerical condition number is compared when the smallest singular value is
well resolved; exact or numerically singular matrices may return infinity.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('src/kinds.f90', 'src/linalg.f90', 'src/linalg_svd.f90')


def make_cases(count: int, seed: int, max_order: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    sizes = sorted({1, 2, 3, 5, 7, 11, min(20, max_order), max_order})
    result = []
    for case in range(count):
        m, n = int(rng.choice(sizes)), int(rng.choice(sizes))
        k = min(m, n)
        is_complex = bool(case % 2)

        def noise(shape: tuple[int, ...]) -> np.ndarray:
            v = rng.standard_normal(shape)
            if is_complex:
                v = v + 1j * rng.standard_normal(shape)
            return v

        kind = case % 9
        if kind in (1, 2, 3, 4, 5, 6):
            u, _ = linalg.qr(noise((m, k)), mode='economic')
            v, _ = linalg.qr(noise((n, k)), mode='economic')
            if kind == 1:
                singular = np.geomspace(1.0, 1.0e-8, k)
            elif kind == 2:
                singular = np.ones(k)
            elif kind == 3:
                singular = np.zeros(k)
            elif kind == 4:
                singular = np.r_[np.ones(max(1, k // 2)), np.zeros(k - max(1, k // 2))]
            elif kind == 5:
                singular = np.linspace(0.4, 1.6, k)
            else:
                singular = np.geomspace(1, 1e-11, k)
            a = (u * singular) @ v.conj().T
            if kind == 3:
                a[:] = 0
        else:
            a = noise((m, n)) / np.sqrt(k)
            if kind == 7:
                if m >= n:
                    a[:, -1] = 0
                else:
                    a[-1, :] = 0
        if case % 43 == 0:
            a *= 1e-105
        elif case % 47 == 0:
            a *= 1e105
        result.append(np.asarray(a))
    return result


def compile_driver(tmp: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    exe = tmp / ('svd_condition_' + mode)
    cmd = ['gfortran', '-std=f2018', *flags,
           *(str(ROOT / src) for src in SOURCES),
           str(ROOT / 'tools/linalg_reference/driver_svd_condition.f90'),
           '-o', str(exe)]
    proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError('Fortran compilation failed:\n' + proc.stderr)
    return exe


def run(exe: Path, cases: list[np.ndarray], tolerances: tuple[float, float]) -> None:
    lines = [str(len(cases))]
    for a in cases:
        s = linalg.svdvals(a)
        tolerance = 2e-4 * float(s[0])
        m, n = a.shape
        lines.append(f'{2 if np.iscomplexobj(a) else 1} {m} {n} {tolerance:.17e}')
        lines.extend(f'{z.real:.17e} {z.imag:.17e}' for z in a.ravel(order='F'))
    proc = subprocess.run([str(exe)], input='\n'.join(lines) + '\n',
                          text=True, capture_output=True, timeout=240)
    if proc.returncode:
        raise RuntimeError(f'Fortran driver failed: {proc.stderr[-2000:]}')
    output = iter(proc.stdout.splitlines())
    worst_spectrum = 0.0
    worst_condition = 0.0
    for index, a in enumerate(cases):
        m, n = a.shape
        k = min(m, n)
        try:
            sts, stc, strank, sttol, rank, rank_tol = map(int, next(output).split())
            condition = float(next(output))
            sigma = np.array([float(next(output)) for _ in range(k)])
        except (StopIteration, ValueError) as exc:
            raise AssertionError(f'case {index}: malformed driver output') from exc
        if (sts, stc, strank, sttol) != (0, 0, 0, 0):
            raise AssertionError(f'case {index}: statuses {sts, stc, strank, sttol}')
        ref = linalg.svdvals(a)
        scale = max(float(np.linalg.norm(ref)), np.finfo(float).tiny)
        error_s = float(np.linalg.norm(sigma - ref) / scale)
        worst_spectrum = max(worst_spectrum, error_s)
        if error_s > tolerances[0]:
            raise AssertionError(f'case {index}: spectrum error {error_s}, shape {a.shape}')
        ref_rank = int(np.linalg.matrix_rank(a))
        tol = float(ref[0]) * 2e-4
        ref_tol_rank = int(np.linalg.matrix_rank(a, tol=tol))
        if (rank, rank_tol) != (ref_rank, ref_tol_rank):
            raise AssertionError(f'case {index}: rank {(rank, rank_tol)} vs '
                                 f'{(ref_rank, ref_tol_rank)} shape={a.shape}')
        # A near-zero singular value is not reliably resolved by double
        # precision Jacobi or by an independent LAPACK SVD.
        if ref[-1] / max(ref[0], np.finfo(float).tiny) >= 1e-10:
            ref_condition = float(ref[0] / ref[-1])
            err_c = abs(condition - ref_condition) / ref_condition
            worst_condition = max(worst_condition, err_c)
            if not np.isfinite(condition) or err_c > tolerances[1]:
                raise AssertionError(f'case {index}: cond {condition} vs {ref_condition}')
        elif not (np.isinf(condition) or (condition >= 1e9)):
            raise AssertionError(f'case {index}: unexpectedly small condition {condition}')
    if next(output, None) is not None:
        raise AssertionError('Unexpected extra output')
    print(f'{len(cases)} SciPy/NumPy cases passed, '
          f'max rel spectrum error {worst_spectrum:.3g}, '
          f'max resolved cond error {worst_condition:.3g}')


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=['both', 'checked', 'optimized'], default='both')
    p.add_argument('--cases', type=int, default=250)
    p.add_argument('--seed', type=int, default=903)
    p.add_argument('--max-order', type=int, default=35)
    p.add_argument('--spectrum-tol', type=float, default=5e-10)
    p.add_argument('--cond-tol', type=float, default=2e-5)
    args = p.parse_args()
    cases = make_cases(args.cases, args.seed, args.max_order)
    with tempfile.TemporaryDirectory(prefix='scifort_svd_condition_') as d:
        tmp = Path(d)
        for mode in ('checked', 'optimized') if args.mode == 'both' else (args.mode,):
            print(mode + ':', end=' ', flush=True)
            run(compile_driver(tmp, mode), cases, (args.spectrum_tol, args.cond_tol))


if __name__ == '__main__':
    main()
