#!/usr/bin/env python3
"""Cross-check SciFort median/Wiener filters against SciPy 1.17.

Compile the driver from src/kinds.f90, src/signal_denoise.f90, and
 tools/signal_reference/denoise_reference_driver.f90 before running.
Usage: python compare_signal_denoise_scipy.py /path/to/denoise_reference_driver
"""
from __future__ import annotations

import subprocess
import sys
import warnings

import numpy as np
import scipy
from scipy import signal


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit('usage: compare_signal_denoise_scipy.py FORTRAN_DRIVER')
    rng = np.random.default_rng(20261019)
    cases = []
    for method in (1, 2, 3):  # medfilt, wiener, medfilt2d
        for dtype in (0, 1):  # double, integer
            for ndim in ((2,) if method == 3 else (1, 2, 3)):
                for rep in range(75):
                    shape = tuple(int(z) for z in rng.integers(1, 9, size=ndim))
                    if method in (1, 3):
                        ks = tuple(int(z) for z in rng.choice([1, 3, 5, 7], size=ndim))
                        noise_flag, noise = 0, 0.0
                    else:
                        ks = tuple(int(z) for z in rng.integers(1, 7, size=ndim))
                        noise_flag = int(rep % 3 != 0)
                        noise = float(rng.uniform(0, 5)) if noise_flag else 0.0
                        if rep % 7 == 0:
                            noise_flag, noise = 1, 0.0
                    if dtype:
                        x = rng.integers(-7, 11, size=shape, dtype=np.int32)
                    else:
                        x = rng.normal(size=shape)
                        if rep % 13 == 0:
                            x *= 0.0
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', RuntimeWarning)
                        warnings.simplefilter('ignore', UserWarning)
                        if method == 1:
                            expected = signal.medfilt(x, kernel_size=ks)
                        elif method == 3:
                            expected = signal.medfilt2d(x, kernel_size=ks)
                        else:
                            expected = signal.wiener(x, mysize=ks,
                                noise=noise if noise_flag else None)
                    cases.append((method, dtype, ndim, shape, ks, noise_flag,
                                  noise, x, expected))

    content = [str(len(cases))]
    for method, dtype, ndim, shape, ks, noise_flag, noise, x, _ in cases:
        content.extend([
            f'{method} {dtype} {ndim} {noise_flag}',
            ' '.join(map(str, shape)),
            ' '.join(map(str, ks)),
            repr(noise),
            ' '.join(format(float(v), '.17g') for v in x.ravel(order='F'))
        ])
    run = subprocess.run([sys.argv[1]], input='\n'.join(content) + '\n',
                         text=True, capture_output=True, check=True)
    lines = run.stdout.splitlines()
    if len(lines) != len(cases):
        raise AssertionError(f'got {len(lines)} output lines for {len(cases)} cases: '
                             f'{run.stderr[-500:]}')
    max_abs_error = 0.0
    for idx, (line, case) in enumerate(zip(lines, cases), 1):
        method, dtype, ndim, shape, ks, noise_flag, noise, x, expected = case
        words = line.split()
        if int(words[0]) != 0:
            raise AssertionError(f'Fortran status {words[0]} in case {idx}')
        actual = np.array([float(v.replace('D', 'E')) for v in words[1:]])
        actual = actual.reshape(shape, order='F')
        if np.shape(expected) != np.shape(actual):
            raise AssertionError(f'shape mismatch in case {idx}')
        try:
            np.testing.assert_allclose(actual, expected, atol=2e-11,
                                       rtol=1e-11, equal_nan=True)
        except AssertionError as exc:
            raise AssertionError(
                f'case {idx}: method={method}, dtype={dtype}, shape={shape}, '
                f'kernel={ks}, noise={noise if noise_flag else None}\n{exc}'
            ) from exc
        finite = np.isfinite(actual) & np.isfinite(expected)
        if np.any(finite):
            max_abs_error = max(max_abs_error,
                                float(np.max(np.abs(actual[finite]-expected[finite]))))
    print(f'SciPy {scipy.__version__}: {len(cases)} / {len(cases)} denoise cases passed; '
          f'max absolute error {max_abs_error:.3e}')


if __name__ == '__main__':
    main()
