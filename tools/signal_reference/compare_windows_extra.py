#!/usr/bin/env python3
"""Validate SciFort's additional get_window families against SciPy 1.17.

Usage: python tools/signal_reference/compare_windows_extra.py
Requires gfortran, NumPy and SciPy; compiles a standalone reference driver.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import subprocess

import numpy as np
from scipy.signal import get_window

ROOT = Path(__file__).resolve().parents[2]


def make_cases():
    names = [
        'boxcar', 'hann', 'hamming', 'blackman', 'bartlett', 'flattop',
        'triang', 'parzen', 'bohman', 'blackmanharris', 'nuttall',
        'barthann', 'cosine', 'lanczos', 'tukey', 'exponential',
        'kaiser', 'gaussian', 'general_gaussian', 'general_hamming',
    ]
    shapes = {
        'tukey': [None, 0., 0.1, 0.5, 0.8, 1., 2.],
        'exponential': [None, 0.25, 1., 3.],
        'kaiser': [0., 2., 5.5, 12., 40., 70., 100.],
        'gaussian': [0.2, 0.5, 2., 12.],
        'general_gaussian': [(0.5, 1.5), (1., 3.), (2.5, 7.)],
        'general_hamming': [0.25, 0.54, 0.7],
    }
    cases = []
    for name in names:
        for param in shapes.get(name, [None]):
            for n in (1, 2, 3, 4, 5, 6, 7, 8, 9, 16, 33, 64, 127):
                for periodic in (True, False):
                    if name == 'exponential':
                        scipy_arg = name if param is None else (name, None, param)
                    elif isinstance(param, tuple):
                        scipy_arg = (name, *param)
                    else:
                        scipy_arg = name if param is None else (name, param)
                    cases.append((name, n, periodic, param, scipy_arg))
    return cases


def main():
    cases = make_cases()
    with TemporaryDirectory(prefix='scifort_windows_') as td:
        td = Path(td)
        exe = td / ('driver.exe' if os.name == 'nt' else 'driver')
        cmd = [
            'gfortran', '-O2', '-ffree-line-length-none', '-J' + str(td),
            '-I' + str(td), str(ROOT / 'src/kinds.f90'),
            str(ROOT / 'src/signal_windows_extra.f90'),
            str(ROOT / 'src/signal_spectral.f90'),
            str(ROOT / 'tools/signal_reference/window_reference_driver.f90'),
            '-o', str(exe),
        ]
        subprocess.run(cmd, check=True, cwd=td, capture_output=True, text=True)
        lines = [str(len(cases))]
        for name, n, periodic, param, _ in cases:
            v = param if isinstance(param, tuple) else (() if param is None else (param,))
            p1 = v[0] if len(v) >= 1 else 0.
            p2 = v[1] if len(v) >= 2 else 0.
            lines.append(f'{name} {n} {int(periodic)} {int(len(v)>=1)} {p1} {int(len(v)>=2)} {p2}')
        output = subprocess.run(
            [str(exe)], input='\n'.join(lines)+'\n', check=True,
            cwd=td, capture_output=True, text=True,
        ).stdout.splitlines()

    index = 0
    largest = 0.
    failed = []
    for case in cases:
        name, n, periodic, param, scipy_arg = case
        fields = output[index].split()
        index += 1
        status, got_n = map(int, fields)
        actual = np.array([float(s) for s in output[index:index+got_n]])
        index += got_n
        expected = get_window(scipy_arg, n, fftbins=periodic)
        difference = np.max(abs(actual-expected), initial=0.)
        largest = max(largest, difference)
        if status or got_n != n or not np.allclose(actual, expected, atol=3.e-11, rtol=2.e-12):
            failed.append((case, status, difference, actual, expected))
    if index != len(output):
        raise AssertionError('Unconsumed driver output')
    print(f'{len(cases)-len(failed)}/{len(cases)} SciPy window comparisons pass')
    print(f'Maximum absolute window error: {largest:.3g}')
    for case in failed[:8]:
        print('FAIL', case)
    if failed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
