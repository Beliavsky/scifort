#!/usr/bin/env python3
"""Compare arbitrary-length SciFort FFT with SciPy's FFT implementation.

Usage: python compare_signal_fft_bluestein.py /path/to/driver
"""
import sys
import subprocess
from pathlib import Path

import numpy as np
from scipy.fft import fft


def run_driver(exe, x):
    data = f"{len(x)}\n" + "".join(f"{v.real:.17e} {v.imag:.17e}\n" for v in x)
    result = subprocess.run([str(exe)], input=data, text=True, capture_output=True,
                            check=True, timeout=20)
    a = np.fromstring(result.stdout, sep=" ")
    assert a.size == 2 * len(x), (len(x), a.size, result.stderr)
    return a.reshape(-1, 2)[:, 0] + 1j * a.reshape(-1, 2)[:, 1]


def main():
    exe = Path(sys.argv[1]).resolve()
    rng = np.random.default_rng(20261026)
    sizes = [1, 2, 3, 5, 7, 11, 15, 16, 17, 19, 23, 31, 33, 35, 45,
             63, 64, 65, 97, 101, 127, 129, 257, 509, 1009, 2053, 4093, 8191, 16381]
    cases, values, largest_abs, largest_relative = 0, 0, 0., 0.
    for n in sizes:
        k = np.arange(n, dtype=np.float64)
        cases_x = [np.ones(n, dtype=np.complex128),
                   np.exp(2j*np.pi*7*k/n) + 0.13j,
                   rng.standard_normal(n) + 1j*rng.standard_normal(n),
                   rng.standard_normal(n) + 1j*rng.standard_normal(n)]
        for x in cases_x:
            actual, expected = run_driver(exe, x), fft(x)
            err = float(np.max(np.abs(actual - expected)))
            scaled = err/max(1., float(np.max(np.abs(expected))))
            largest_abs = max(largest_abs, err)
            largest_relative = max(largest_relative, scaled)
            if err > 1.e-10 * max(1., float(np.max(np.abs(expected)))):
                raise AssertionError(f"n={n} absolute={err}, scaled={scaled}")
            cases += 1
            values += n
    print(f"PASS {cases} SciPy FFT comparisons; {values} outputs; "
          f"max absolute {largest_abs:.3e}; max scaled {largest_relative:.3e}")


if __name__ == '__main__':
    main()
