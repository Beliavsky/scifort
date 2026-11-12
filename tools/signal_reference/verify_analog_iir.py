#!/usr/bin/env python3
"""Randomized comparison of SciFort analog IIR designs with SciPy 1.17.

Requires gfortran, NumPy and SciPy; adds no dependency to SciFort itself.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
import numpy as np
import scipy
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
MODULES = [
    'kinds', 'signal_analog_design', 'signal_dpss', 'signal_windows_extra', 'signal_spectral',
    'signal_design', 'signal_elliptic', 'signal_iir_extra',
]
SEED = 20261016
COUNT = 1200


def main():
    rng = np.random.default_rng(SEED)
    cases = []
    for i in range(COUNT):
        family = ('butter', 'cheby1', 'cheby2')[i % 3]
        mode = ('lowpass', 'highpass', 'bandpass', 'bandstop')[(i // 3) % 4]
        order = int(rng.integers(1, 9))
        db = float(rng.uniform(0.15, 3.5) if family == 'cheby1'
                   else rng.uniform(15.0, 75.0))
        w1 = 10.0 ** float(rng.uniform(-0.5, 1.0))
        if mode.startswith('band'):
            w2 = w1 * float(rng.uniform(1.4, 4.0))
            critical = [w1, w2]
        else:
            w2 = 0.0
            critical = w1
        ref = getattr(signal, family)
        args = [order]
        if family != 'butter':
            args.append(db)
        args.append(critical)
        b, a = ref(*args, btype=mode, analog=True)
        cases.append((family, mode, order, db, w1, w2, b, a))

    data = ''.join(f'{fam} {n} {mode} {db:.17e} {w1:.17e} {w2:.17e}\n'
                   for fam, mode, n, db, w1, w2, _, _ in cases)
    max_coefficient_error = 0.0
    max_response_error = 0.0
    max_meaningful_relative_error = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-analog-iir-') as temp:
        work = Path(temp)
        objs = []
        flags = ['-std=f2018', '-Wimplicit-interface', '-Werror=implicit-interface',
                 '-O0', '-fcheck=all', '-fbacktrace']
        for mod in MODULES:
            path = ROOT / 'src' / f'{mod}.f90'
            obj = work / f'{mod}.o'
            subprocess.run(['gfortran', *flags, '-I', str(work), '-J', str(work),
                            '-c', str(path), '-o', str(obj)], cwd=work, check=True)
            objs.append(str(obj))
        driver = ROOT / 'tools/signal_reference/analog_iir_driver.f90'
        exe = work / 'analog_iir_driver'
        subprocess.run(['gfortran', *flags, '-I', str(work), str(driver),
                        *objs, '-o', str(exe)], cwd=work, check=True)
        results = subprocess.run([str(exe)], input=data, text=True,
                                 capture_output=True, cwd=work, check=True)
        lines = iter(results.stdout.splitlines())
        for idx, (family, mode, order, db, w1, w2, bref, aref) in enumerate(cases):
            status, nb, na = map(int, next(lines).split())
            if status != 0:
                raise AssertionError(f'case {idx}: invalid filter for {family} {mode}')
            b = np.fromstring(next(lines), sep=' ')
            a = np.fromstring(next(lines), sep=' ')
            if len(b) != nb or len(a) != na:
                raise AssertionError(f'case {idx}: malformed filter lengths')
            # Coefficients are derived by complex root multiplication;
            # near-zero terms can have relative cancellation error.
            for observed, reference in ((b, bref), (a, aref)):
                if observed.shape != reference.shape:
                    raise AssertionError(f'case {idx}: shape mismatch')
                relative = np.max(np.abs(observed-reference) /
                                  max(1.0, np.max(np.abs(reference))))
                max_coefficient_error = max(max_coefficient_error, relative)
                if relative > 5e-10:
                    raise AssertionError(f'case {idx} ({family} {mode}) coefficient error {relative}')

            # Check stable, non-singular parts of the analog transfer function.
            freqs = np.geomspace(max(w1 / 20.0, 1e-3),
                                 max(w1, w2) * 20.0, 15)
            _, href = signal.freqs(bref, aref, worN=freqs)
            _, h = signal.freqs(b, a, worN=freqs)
            abs_err = np.abs(h-href)
            max_response_error = max(max_response_error, float(np.max(abs_err)))
            # Absolute differences near transmission zeros can dominate the
            # relative error of transfer-function polynomials.
            eligible = np.abs(href) > 1e-3
            if np.any(eligible):
                rel_err = float(np.max(abs_err[eligible] / np.abs(href[eligible])))
                max_meaningful_relative_error = max(max_meaningful_relative_error, rel_err)
            else:
                rel_err = 0.0
            if not np.allclose(h, href, rtol=3e-8, atol=5e-9):
                raise AssertionError(f'case {idx} ({family} {mode}) response error {rel_err}')
        if next(lines, None) is not None:
            raise AssertionError('Unexpected extra driver output')
    print(f'PASS: {COUNT}/{COUNT} randomized analog IIR designs vs SciPy {scipy.__version__}')
    print(f'Max scaled coefficient discrepancy: {max_coefficient_error:.3e}')
    print(f'Max complex response absolute discrepancy: {max_response_error:.3e}')
    print(f'Max relative response error where |H| > 1e-3: {max_meaningful_relative_error:.3e}')


if __name__ == '__main__':
    main()
