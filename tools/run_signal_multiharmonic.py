#!/usr/bin/env python3
"""Validate SciFort's weighted multi-harmonic Lomb-Scargle periodogram.

The one-harmonic specialization is checked against SciPy's public
lombscargle. Multiple harmonics are checked independently with weighted
NumPy least squares, not by calling any SciFort helper.

Run from any directory: python tools/run_signal_multiharmonic.py
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.signal import lombscargle

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261028)


def make_cases():
    cases = []
    for j in range(384):
        terms = 1 + j % 4
        n = int(RNG.integers(32, 112))
        m = int(RNG.integers(3, 11))
        x = np.sort(RNG.uniform(-5.0, 17.0, n))
        freqs = RNG.uniform(0.23, 7.5, m)
        floating = j % 2 == 0
        center = j % 3 == 0
        weighted = j % 5 != 0
        normalized = j % 2 == 1
        weights = RNG.uniform(0.1, 2.5, n)
        if j % 7 == 0:
            weights[:n//8] = 0.0
        if j % 11 == 0:
            weights *= 1e4
        omega = float(freqs[0])
        y = RNG.normal(0.0, 0.15, n)
        for h in range(1, terms + 1):
            y += RNG.normal() * np.sin(h*omega*x)
            y += RNG.normal() * np.cos(h*omega*x)
        y += RNG.uniform(-2.0, 3.0)
        cases.append((x, y, freqs, weights, terms, floating,
                      center, weighted, normalized))
    return cases


def expected(case):
    x, y, freqs, weights, terms, floating, center, weighted, normalized = case
    yy = y.copy()
    if center:
        yy -= np.mean(y)
    w = weights.copy() if weighted else np.ones_like(yy)
    w /= np.sum(w)
    if floating:
        yy -= np.sum(w*yy)
    target = yy*np.sqrt(w)
    denominator = np.sum(target**2)
    results = []
    for frequency in freqs:
        columns = []
        for h in range(1, terms+1):
            for col in [np.cos(h*frequency*x), np.sin(h*frequency*x)]:
                if floating:
                    col = col - np.sum(w*col)
                columns.append(col*np.sqrt(w))
        coefficients, _, _, _ = np.linalg.lstsq(np.column_stack(columns), target, rcond=None)
        improvement = np.sum((np.column_stack(columns)@coefficients)**2)
        results.append(improvement/denominator if normalized else 0.5*len(y)*improvement)
    numpy_expected = np.asarray(results)
    if terms == 1:
        kw = {'precenter': center, 'floating_mean': floating, 'normalize':
              'normalize' if normalized else 'power'}
        if weighted:
            kw['weights'] = weights
        scipy_expected = lombscargle(x, y, freqs, **kw)
        np.testing.assert_allclose(numpy_expected, scipy_expected,
                                   atol=7e-10, rtol=7e-10)
    return numpy_expected


def input_cases(cases):
    chunks = [str(len(cases))]
    for x, y, f, w, terms, floating, center, weighted, normalized in cases:
        chunks.append(f'{len(x)} {len(f)} {terms} {int(center)} '
                      f'{int(floating)} {int(weighted)} {int(normalized)}')
        for vector in (x, y, f, w):
            chunks.append(' '.join(f'{val:.17e}' for val in vector))
    return '\n'.join(chunks)+'\n'


def check_driver(exe, cases, refs):
    proc = subprocess.run([str(exe)], input=input_cases(cases), text=True,
                          capture_output=True, check=True)
    lines = iter(proc.stdout.splitlines())
    worst = 0.0
    count = 0
    for idx, ref in enumerate(refs):
        status, length = map(int, next(lines).split())
        if status != 0 or length != len(ref):
            raise AssertionError((idx, status, length))
        got = np.array([float(next(lines)) for _ in range(length)])
        worst = max(worst, float(np.max(np.abs(got-ref))))
        count += length
        np.testing.assert_allclose(got, ref, atol=4e-10, rtol=4e-10,
                                   err_msg=f'case {idx}, terms {cases[idx][4]}')
    assert next(lines, None) is None
    return count, worst


def build_and_test(path, optimize, cases, refs):
    path.mkdir()
    flags = ['-std=f2018', '-O3', '-march=native'] if optimize else [
        '-std=f2018', '-O0', '-g', '-fcheck=all', '-fbacktrace']
    base = ['gfortran', *flags, '-Wall', '-Wextra', '-Werror=implicit-interface',
            '-J', str(path), '-I', str(path)]
    src = [str(ROOT/'src/kinds.f90'), str(ROOT/'src/signal_irregular.f90')]
    for idx, test in enumerate([ROOT/'test/test_signal_multiharmonic.f90',
                                ROOT/'tools/reference_signal_multiharmonic.f90']):
        exe = path/f'driver_{idx}'
        subprocess.run([*base, *src, str(test), '-o', str(exe)],
                       cwd=path, check=True, capture_output=True, text=True)
        if idx == 0:
            result = subprocess.run([str(exe)], check=True, capture_output=True, text=True)
            print(' ',result.stdout.strip())
        else:
            count, worst = check_driver(exe, cases, refs)
            print(f'  PASS {len(cases)} reference cases / {count} bins; max |error| = {worst:.3e}')
    # The legacy irregular-sampling tests use only these two Fortran modules
    # in their direct-module counterpart. Its full facade test is part of CI.
    return count, worst


def main():
    cases = make_cases()
    refs = [expected(case) for case in cases]
    with tempfile.TemporaryDirectory(prefix='scifort_multiharmonic_') as tmp:
        for optimize in [False, True]:
            print('optimized (-O3):' if optimize else 'checked (-fcheck=all):')
            build_and_test(Path(tmp)/('opt' if optimize else 'checked'), optimize, cases, refs)
    print('Validated NumPy weighted least squares and SciPy single harmonic.')


if __name__ == '__main__':
    main()
