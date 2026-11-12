#!/usr/bin/env python3
"""Independent numerical comparison of mixed-radix FFT plans and SciPy FFT."""
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.fft import fft


def run_fortran(exe, values, repeat):
    payload = f"{len(values)} {repeat}\n" + "".join(
        f"{x.real:.17e} {x.imag:.17e}\n" for x in values
    )
    proc = subprocess.run([str(exe)], input=payload, text=True,
                          capture_output=True, check=True, timeout=40)
    result = np.fromstring(proc.stdout, sep=" ")
    if len(result) != 2*len(values):
        raise AssertionError((len(values), proc.stderr, len(result)))
    return result.reshape(-1, 2)[:, 0] + 1j*result.reshape(-1, 2)[:, 1]


def main():
    exe = Path(sys.argv[1]).resolve()
    rng = np.random.default_rng(20261026)
    sizes = [18,21,25,27,28,35,42,49,63,75,81,98,105,120,125,147,175,
             210,225,256,343,375,420,480,625,729,960,1000,1024,1225,1728,
             2048,2401,3072,4096,5040,8192,16807, 17,101,1009]
    count = coefficients = 0
    max_abs = max_rel = 0.0
    for n in sizes:
        t = np.arange(n)
        samples = (rng.standard_normal(n)+1j*rng.standard_normal(n),
                   0.2 + np.exp(2j*np.pi*7*t/n),
                   np.ones(n, dtype=np.complex128))
        for x in samples:
            for repeat in (1,3):
                actual = run_fortran(exe, x, repeat)
                reference = fft(x * complex(repeat,0.3*repeat))
                diff = np.max(np.abs(actual-reference))
                rel = diff/max(1.0, float(np.max(np.abs(reference))))
                max_abs = max(max_abs, float(diff))
                max_rel = max(max_rel,float(rel))
                if rel > 5e-12:
                    raise AssertionError(f"FFT error at n={n}, repeat={repeat}: {rel:g}")
                coefficients += n
                count += 1
    print(f"PASS {count} SciPy FFT comparisons, {coefficients} coefficients; "
          f"max absolute={max_abs:.3e}, max scaled={max_rel:.3e}")


if __name__ == "__main__":
    main()
