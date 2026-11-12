#!/usr/bin/env python3
"""Verify near-rank-deficient FIRLS designs against SciPy 1.17+.

The coefficients of nearly singular designs are not unique at double
precision, so compare weighted in-band responses rather than individual taps.
Requires NumPy, SciPy, and gfortran; no external Fortran library or FPM.
"""
from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]


def make_cases() -> list[tuple[int, np.ndarray, np.ndarray, np.ndarray, float]]:
    rng = np.random.default_rng(20261018)
    cases: list[tuple[int, np.ndarray, np.ndarray, np.ndarray, float]] = []
    for n in (9, 15, 31, 51):
        for width in (.03, .005, .001, .0001, .00001):
            cases.append((n, np.array([.4, .4 + width]),
                          np.ones(2), np.ones(1), 2.))

    for i in range(350):
        n = int(rng.choice([9, 13, 19, 29, 41, 57]))
        nb = int(rng.integers(1, 4))
        starts = np.sort(rng.uniform(.07, .88, nb))
        if nb > 1 and np.min(np.diff(starts)) < .025:
            continue
        widths = 10**rng.uniform(-4.5, -1.3, nb)
        b = np.empty(2*nb)
        for j in range(nb):
            b[2*j:2*j+2] = (starts[j], min(starts[j] + widths[j], .995))
        if np.any(b[2:-1:2] < b[1:-2:2]):
            continue
        desired = rng.uniform(0, 1.5, 2*nb)
        if i % 3 == 0:
            desired[1::2] = desired[::2]
        weight = 10**rng.uniform(-4, 4, nb)
        fs = float(rng.choice([2., 40., 900.]))
        cases.append((n, b, desired, weight, fs))
    return cases


def input_for(cases: list[tuple[int, np.ndarray, np.ndarray, np.ndarray, float]]) -> str:
    lines = []
    for n, bands, desired, weight, rate in cases:
        lines += ['LS', f'{n} {len(bands)} {rate:.17g}']
        for array in (bands*rate/2, desired, weight):
            lines.append(' '.join(f'{v:.17e}' for v in array))
    return '\n'.join(lines) + '\n'


def response_error(h: np.ndarray, n: int, bands: np.ndarray,
                   desired: np.ndarray, weights: np.ndarray) -> float:
    errors = []
    for j, (low, high) in enumerate(bands.reshape(-1, 2)):
        omega = np.linspace(low, high, 129)*np.pi
        actual = signal.freqz(h, worN=omega)[1] * np.exp(1j*omega*(n-1)/2)
        wanted = np.linspace(desired[2*j], desired[2*j+1], len(omega))
        errors.append(weights[j] * (actual.real - wanted))
    return float(np.sqrt(np.mean(np.concatenate(errors)**2)))


def main() -> None:
    cases = make_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-firls-conditioned-') as tmp:
        work = Path(tmp)
        flags = ['-std=f2018', '-O2', '-fcheck=all', '-Werror=implicit-interface']
        objects = []
        for name in ('kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral', 'signal_fir_advanced'):
            obj = work / f'{name}.o'
            subprocess.run(['gfortran', *flags, '-I', str(work), '-J', str(work),
                            '-c', str(ROOT/'src'/f'{name}.f90'), '-o', str(obj)],
                           check=True, cwd=work)
            objects.append(str(obj))
        executable = work / 'firls_driver'
        subprocess.run(['gfortran', *flags, '-I', str(work), '-J', str(work),
                        str(ROOT/'tools/signal_reference/fir_advanced_driver.f90'),
                        *objects, '-o', str(executable)], check=True, cwd=work)
        process = subprocess.run([str(executable)], input=input_for(cases),
                                 check=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, text=True, cwd=work)

    tokens = iter(process.stdout.split())
    max_excess = 0.
    largest_error = 0.
    for case_index, (n, bands, desired, weight, fs) in enumerate(cases):
        try:
            status, count = int(next(tokens)), int(next(tokens))
            taps = np.array([float(next(tokens)) for _ in range(count)])
        except (StopIteration, ValueError) as exc:
            raise AssertionError(f'incomplete response, case {case_index}') from exc
        if status != 0 or count != n or not np.all(np.isfinite(taps)):
            raise AssertionError(f'case {case_index}: invalid solution, status {status}')
        if not np.allclose(taps, taps[::-1], atol=1e-12, rtol=0):
            raise AssertionError(f'case {case_index}: coefficients are not symmetric')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            reference = signal.firls(n, bands*fs/2, desired, weight=weight, fs=fs)
        if not np.all(np.isfinite(reference)):
            raise AssertionError(f'nonfinite SciPy reference for case {case_index}')
        normalized_weights = weight/max(weight)
        got_error = response_error(taps, n, bands, desired, normalized_weights)
        ref_error = response_error(reference, n, bands, desired, normalized_weights)
        # Least-squares solutions can be nonunique when the Gram matrix
        # is numerically rank deficient; small absolute deviations are allowed.
        tolerance = 3*ref_error + 1.e-5*(1. + max(desired))
        if got_error > tolerance:
            raise AssertionError(f'case {case_index}: weighted RMS {got_error:.3e} '
                                 f'exceeds SciPy {ref_error:.3e} plus slack {tolerance:.3e}')
        max_excess = max(max_excess, got_error - ref_error)
        largest_error = max(largest_error, got_error)
    try:
        next(tokens)
        raise AssertionError('trailing output')
    except StopIteration:
        pass
    print(f'PASS: {len(cases)} rank-sensitive FIRLS cases against SciPy '
          f'{__import__("scipy").__version__}')
    print('  narrow constant-gain cases: 20')
    print(f'  maximum absolute excess weighted RMS: {max_excess:.3e}')
    print(f'  maximum weighted RMS (includes steep impossible targets): {largest_error:.3e}')


if __name__ == '__main__':
    main()
