#!/usr/bin/env python3
"""Compare SciFort convolve2d/correlate2d to scipy.signal (seeded cases).

Compile the companion driver using scifort_kinds, scifort_signal_convolve,
and scifort_signal_convolve2d, then pass --exe /path/to/driver.
"""
from __future__ import annotations

import argparse
import subprocess
import numpy as np
from scipy import signal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', required=True)
    parser.add_argument('--cases', type=int, default=840)
    args = parser.parse_args()
    rng = np.random.default_rng(20261018)
    modes = ('full', 'same', 'valid')
    boundaries = ('fill', 'wrap', 'symm')
    reference = []
    lines = []
    for case in range(args.cases):
        method = case % 2
        mode = modes[(case // 2) % 3]
        boundary = boundaries[(case // 6) % 3]
        # Include singleton dimensions, even kernels, and both valid orderings.
        n1, n2, m1, m2 = rng.integers(1, 7, size=4).tolist()
        if mode == 'valid':
            if (case // 2) % 2 == 0:
                n1, m1 = min(n1, m1), max(n1, m1)
                n2, m2 = min(n2, m2), max(n2, m2)
            else:
                n1, m1 = max(n1, m1), min(n1, m1)
                n2, m2 = max(n2, m2), min(n2, m2)
        a = rng.standard_normal((n1, n2)) + 1j * rng.standard_normal((n1, n2))
        b = rng.standard_normal((m1, m2)) + 1j * rng.standard_normal((m1, m2))
        fill = complex(rng.standard_normal(), rng.standard_normal())
        f = (signal.correlate2d if method else signal.convolve2d)
        expected = f(a, b, mode=mode, boundary=boundary, fillvalue=fill)
        reference.append((expected,case,method,mode,boundary))
        lines.append(f'{method + 1} {n1} {n2} {m1} {m2} {modes.index(mode)+1} {boundaries.index(boundary)+1}\n')
        lines.append(f'{fill.real:.17g} {fill.imag:.17g}\n')
        for arr in (a,b):
            for value in arr.ravel(order='F'):
                lines.append(f'{value.real:.17g} {value.imag:.17g}\n')
    result = subprocess.run([args.exe], input=''.join(lines), text=True,
                            check=True, capture_output=True)
    raw = iter(result.stdout.splitlines())
    error = 0.0
    for expected, case, method, mode, boundary in reference:
        status, rows, cols = map(int, next(raw).split())
        assert status == 0, (case, status)
        assert (rows, cols) == expected.shape, (case, (rows, cols), expected.shape)
        v = np.empty(rows*cols, dtype=complex)
        for i in range(len(v)):
            pair = next(raw).split()
            v[i] = complex(float(pair[0]), float(pair[1]))
        got = v.reshape((rows, cols), order='F')
        diff = np.max(np.abs(got-expected))
        scaled = diff / max(1.,float(np.max(np.abs(expected))))
        if not np.isfinite(scaled) or scaled > 2e-12:
            raise AssertionError((case, method,mode,boundary, got,expected, scaled))
        error = max(error,float(scaled))
    print(f'PASS: {len(reference)} SciPy 1.17 convolve2d/correlate2d cases; '
          f'max scaled error {error:.3g}')


if __name__ == '__main__':
    main()
