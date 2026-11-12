#!/usr/bin/env python3
"""Randomized SciPy comparisons for SciFort dense polar and exponential Frechet algorithms.

Requires Python NumPy, SciPy, and GNU Fortran; independent reference via
scipy.linalg.polar / scipy.linalg.expm_frechet. No plotting is performed.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.linalg import expm, expm_frechet, polar

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [ROOT / "src" / name for name in
           ("kinds.f90", "linalg.f90", "linalg_polar.f90", "linalg_expm_frechet.f90")]
DRIVER = ROOT / "tools/linalg_reference/driver_polar_frechet.f90"


def make_cases(ncases: int, seed: int):
    rng = np.random.default_rng(seed)
    cases = []
    for k in range(ncases):
        dtype = 1 if (k // 2) % 2 == 0 else 2
        if k % 2 == 0:
            m = int(rng.integers(1, 9))
            n = int(rng.integers(1, 9))
            side = "left" if k % 4 == 0 else "right"
            a = rng.standard_normal((m, n))
            if dtype == 2:
                a = a + 1j*rng.standard_normal((m, n))
            if k % 10 == 0:
                a[:] = 0
            elif k % 10 == 2 and n > 1:
                a[:, -1] = a[:, 0] * (0.25 + (0.2j if dtype == 2 else 0))
            elif k % 10 == 4:
                a *= 1e-8
            elif k % 10 == 6:
                a *= 200
            cases.append((1, dtype, a, None, side))
        else:
            n = int(rng.integers(1, 8))
            a = rng.standard_normal((n, n))
            e = rng.standard_normal((n, n))
            if dtype == 2:
                a = a + 1j*rng.standard_normal((n, n))
                e = e + 1j*rng.standard_normal((n, n))
            a *= (0.0, 0.03, 0.4, 1, 3, 7, 12)[k % 7]
            if k % 13 == 0:
                e[:] = 0
            elif k % 11 == 0:
                e = np.eye(n, dtype=a.dtype)
            cases.append((2, dtype, a, e, "right"))
    return cases


def make_input(cases):
    lines = [str(len(cases))]
    for op, dtype, a, e, side in cases:
        lines.append(f"{op} {dtype} {a.shape[0]} {a.shape[1]} {int(side == 'left')}")
        for mat in (a,) if e is None else (a, e):
            for val in mat.flat:
                lines.append(f"{val.real:.17e} {val.imag:.17e}")
    return "\n".join(lines) + "\n"


def parse_output(out, cases):
    lines = iter(out.splitlines())
    results = []
    for op, dtype, a, e, side in cases:
        header = next(lines).split()
        status = int(header[0])
        if op == 1:
            pside = a.shape[1] if side == 'right' else a.shape[0]
            shapes = (a.shape, (pside, pside))
        else:
            shapes = (a.shape, a.shape)
        result = []
        for m, n in shapes:
            vals = [complex(*map(float, next(lines).split())) for _ in range(m*n)]
            result.append(np.asarray(vals).reshape(m, n))
        results.append((status, result, tuple(map(float, header[1:]))))
    return results


def check(cases, results):
    max_polar = max_frechet = 0.0
    counts = [0, 0]
    for idx, ((op, dtype, a, e, side), (status, outputs, metrics)) in enumerate(zip(cases, results)):
        if status != 0:
            raise AssertionError(f"case {idx}, {op=}, {a.shape=}, Fortran status {status}")
        if op == 1:
            u, p = outputs
            uref, pref = polar(a, side=side)
            reconstruction = u @ p if side == 'right' else p @ u
            error = np.linalg.norm(reconstruction - a) / max(np.linalg.norm(a), 1e-300)
            perr = np.linalg.norm(p-pref) / max(1, np.linalg.norm(pref))
            gram = u.conj().T @ u if a.shape[0] >= a.shape[1] else u @ u.conj().T
            unit_err = np.linalg.norm(gram-np.eye(min(a.shape)))
            herm_err = np.linalg.norm(p-p.conj().T)
            psd_min = np.linalg.eigvalsh(p).min()
            if max(error, perr, unit_err, herm_err) > 5e-9 or psd_min < -1e-8 * max(1, np.linalg.norm(p)):
                raise AssertionError(f"polar {idx}, shape {a.shape}, side={side}: "
                                     f"reconstruct={error:g}, p={perr:g}, unit={unit_err:g}, "
                                     f"hermitian={herm_err:g}, eigmin={psd_min:g}")
            max_polar = max(max_polar, error, perr, unit_err)
            counts[0] += 1
        else:
            x, d = outputs
            refx, refd = expm_frechet(a, e)
            if idx % 31 == 1:
                n = len(a)
                block = np.zeros((2*n, 2*n), dtype=np.result_type(a, e))
                block[:n, :n] = a
                block[n:, n:] = a
                block[:n, n:] = e
                blockexp = expm(block)
                if np.linalg.norm(d-blockexp[:n,n:]) / max(1, np.linalg.norm(refd)) > 3e-10:
                    raise AssertionError(f'block-exponential mismatch in case {idx}')
            error_x = np.linalg.norm(x-refx) / max(1, np.linalg.norm(refx))
            error_d = np.linalg.norm(d-refd) / max(1, np.linalg.norm(refd))
            if max(error_x, error_d) > 3e-10:
                raise AssertionError(f"frechet {idx}, n={len(a)}: exp={error_x:g}, frechet={error_d:g}")
            max_frechet = max(max_frechet, error_x, error_d)
            counts[1] += 1
    return counts, max_polar, max_frechet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    parser.add_argument('--cases', type=int, default=300, help='total combined polar and Frechet cases per build')
    parser.add_argument('--seed', type=int, default=20261106)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error('--cases must be positive')
    modes = ('checked', 'optimized') if args.mode == 'both' else (args.mode,)
    cases = make_cases(args.cases, args.seed)
    with tempfile.TemporaryDirectory(prefix='scifort-polar-frechet-') as tmp:
        for mode in modes:
            flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
            exe = Path(tmp) / f'pf_{mode}'
            cmd = ['gfortran', *flags, '-ffree-line-length-none', *map(str, SOURCES), str(DRIVER), '-o', str(exe)]
            subprocess.run(cmd, cwd=tmp, check=True, capture_output=True, text=True)
            output = subprocess.run([str(exe)], input=make_input(cases), text=True,
                                    capture_output=True, check=True, cwd=tmp)
            counts, p_error, f_error = check(cases, parse_output(output.stdout, cases))
            print(f'{mode}: {counts[0]} polar / {counts[1]} expm_frechet cases PASS; '
                  f'max relative polar error={p_error:.3e}, Frechet error={f_error:.3e}')


if __name__ == '__main__':
    main()
