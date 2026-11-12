#!/usr/bin/env python3
"""Compare pair-sum radix-3/5/7 FFT kernels with independent SciPy FFTs."""
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.fft import fft


LENGTHS = (21, 25, 27, 35, 49, 63, 75, 105, 125, 147, 175,
           243, 343, 375, 625, 729, 875, 1225, 1701, 2401, 3125,
           5040, 16807)


def call_fortran(executable, x, repeat):
    payload = f"{len(x)} {repeat}\n" + "".join(
        f"{sample.real:.17e} {sample.imag:.17e}\n" for sample in x
    )
    proc = subprocess.run([str(executable)], input=payload, text=True,
                          capture_output=True, check=True, timeout=60)
    data = np.fromstring(proc.stdout, sep=" ")
    if data.size != 2*len(x):
        raise AssertionError(f"n={len(x)}: {data.size} output floats")
    return data[::2] + 1j*data[1::2]


def main():
    executable = Path(sys.argv[1]).resolve()
    rng = np.random.default_rng(20261026)
    max_abs = max_scaled = 0.0
    cases = coefficients = 0
    for n in LENGTHS:
        k = np.arange(n)
        impulse = np.zeros(n, dtype=complex)
        impulse[n//3] = 1.25 - 0.75j
        inputs = (
            rng.normal(size=n) + 1j*rng.normal(size=n),
            impulse,
            np.exp(2j*np.pi*(11*k/n + 0.001*k*k/n)),
        )
        for x in inputs:
            for repeat in (1, 3):
                found = call_fortran(executable, x, repeat)
                reference = fft(x*complex(repeat, 0.3*repeat))
                error = float(np.max(np.abs(found-reference)))
                scaled = error/max(1.0, float(np.max(np.abs(reference))))
                max_abs = max(max_abs, error)
                max_scaled = max(max_scaled, scaled)
                if scaled > 2.0e-12:
                    raise AssertionError(f"n={n}, repeat={repeat}, scaled error={scaled:g}")
                cases += 1
                coefficients += n
    print(f"PASS {cases} independent SciPy butterfly comparisons, "
          f"{coefficients} coefficients; max abs={max_abs:.3e}, "
          f"max scaled={max_scaled:.3e}")


if __name__ == '__main__':
    main()
