"""Validate SciFort check_cola/check_nola against independent SciPy calls.

Usage: python test/compare_signal_overlap_checks_scipy.py [gfortran]
"""
from __future__ import annotations

import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import check_COLA, check_NOLA, windows

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / "src" / name for name in
           ("kinds.f90", "signal_dpss.f90", "signal_windows_extra.f90", "signal_spectral.f90", "signal_istft.f90")]
DRIVER = ROOT / "tools/signal_reference/signal_overlap_checks_reference_driver.f90"


def main() -> None:
    compiler = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FC", "gfortran")
    rng = np.random.default_rng(20261019)
    cases = []
    with tempfile.TemporaryDirectory(prefix="scifort-cola-nola-") as directory:
        tmp = Path(directory)
        exe = tmp / ("overlap.exe" if os.name == "nt" else "overlap")
        subprocess.run([compiler, "-std=f2018", "-O2", "-fcheck=all",
                        "-Werror=implicit-interface", *map(str, SOURCES),
                        str(DRIVER), "-o", str(exe)], cwd=tmp, check=True)
        input_path = tmp / "input.bin"
        output_path = tmp / "output.bin"
        with input_path.open("wb") as output:
            for i in range(560):
                n = int(rng.integers(2, 70))
                no = int(rng.integers(0, n))
                family = i % 8
                if family == 0:
                    w = windows.hann(n, sym=False)
                elif family == 1:
                    w = windows.hamming(n, sym=False)
                elif family == 2:
                    w = windows.hann(n, sym=True)
                elif family == 3:
                    w = windows.boxcar(n)
                elif family == 4:
                    w = windows.blackman(n, sym=False)
                elif family == 5:
                    w = rng.normal(size=n)
                elif family == 6:
                    w = np.ones(n)
                    w[::2] = 0
                else:
                    w = np.r_[np.ones(n // 2), np.zeros(n - n//2)]
                tolerance = (1e-10, 1e-8, 1e-6)[i % 3]
                c = bool(check_COLA(w, n, no, tol=tolerance))
                d = bool(check_NOLA(w, n, no, tol=tolerance))
                output.write(struct.pack("ii", n, no))
                output.write(struct.pack("d", tolerance))
                output.write(np.asarray(w, dtype=np.float64).tobytes())
                cases.append((c, d))
        subprocess.run([str(exe), str(input_path), str(output_path)], check=True)
        result = np.fromfile(output_path, dtype=np.int32).reshape(-1, 2)
        if result.shape[0] != len(cases):
            raise AssertionError("Missing overlap-check output")
        for i, (observed, reference) in enumerate(zip(result, cases)):
            if tuple(map(bool, observed)) != reference:
                raise AssertionError(f"case {i}: SciFort {observed}, SciPy {reference}")
    print(f"SciPy COLA/NOLA comparisons: {2*len(cases)} passed ({len(cases)} windows)")


if __name__ == "__main__":
    main()
