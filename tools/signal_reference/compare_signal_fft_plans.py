#!/usr/bin/env python3
"""Independent reference for repeated planned transforms (SciPy FFT)."""
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.fft import fft


def compare(exe, x, repeats):
    data = f"{len(x)} {repeats}\n" + "".join(f"{v.real:.17e} {v.imag:.17e}\n" for v in x)
    proc = subprocess.run([str(exe)], input=data, capture_output=True,
                          text=True, check=True, timeout=30)
    values = np.fromstring(proc.stdout, sep=" ")
    assert len(values) == 2 * len(x), (len(x), proc.stderr, len(values))
    actual = values.reshape(-1, 2)[:, 0] + 1j * values.reshape(-1, 2)[:, 1]
    expected = fft(x * complex(repeats, 0.3 * repeats))
    error = np.max(np.abs(actual - expected))
    return error, error / max(1., float(np.max(np.abs(expected))))


def main():
    exe = Path(sys.argv[1]).resolve()
    rng = np.random.default_rng(20261026)
    sizes = [1, 2, 3, 5, 11, 16, 17, 19, 31, 33, 63, 64, 65, 97,
             101, 257, 509, 1009, 2053, 4093, 8191, 16381]
    count, values, max_absolute, max_scaled = 0, 0, 0., 0.
    for n in sizes:
        k = np.arange(n, dtype=float)
        signals = [rng.standard_normal(n) + 1j*rng.standard_normal(n),
                   np.exp(2j*np.pi*7*k/n) + 0.1j,
                   np.ones(n, dtype=np.complex128)]
        for x in signals:
            for repeats in (1, 2, 5):
                a, r = compare(exe, x, repeats)
                max_absolute, max_scaled = max(max_absolute, a), max(max_scaled, r)
                if r > 3e-12:
                    raise AssertionError(f"n={n}, repeats={repeats}, scaled error={r}")
                values += n
                count += 1
    print(f"PASS {count} planned FFT SciPy comparisons; {values} output coefficients; "
          f"max absolute {max_absolute:.3e}; max scaled {max_scaled:.3e}")


if __name__ == '__main__':
    main()
