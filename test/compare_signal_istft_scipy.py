"""Compare SciFort ISTFT reconstruction against SciPy 1.17.

Usage: python test/compare_signal_istft_scipy.py [gfortran]
Builds an isolated binary-stream reference driver with GNU Fortran.
"""
from __future__ import annotations

import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import warnings

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "tools/signal_reference/signal_istft_reference_driver.f90"
SOURCES = [ROOT / "src" / name for name in
           ("kinds.f90", "signal_dpss.f90", "signal_windows_extra.f90", "signal_spectral.f90", "signal_istft.f90")]
WINDOWS = ["hann", "hamming", "boxcar", "blackman", "symmetric", "alternating"]


def window_spec(code: int, nperseg: int):
    if code < 4:
        return WINDOWS[code]
    if code == 4:
        return signal.windows.hann(nperseg, sym=True)
    w = np.ones(nperseg)
    w[::2] = 0.0
    return w


def main() -> None:
    rng = np.random.default_rng(20261019)
    compiler = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FC", "gfortran")
    maximum = 0.0
    count = 0
    warning_cases = 0
    with tempfile.TemporaryDirectory(prefix="scifort-istft-") as directory:
        tmp = Path(directory)
        executable = tmp / ("istft.exe" if os.name == "nt" else "istft")
        subprocess.run([compiler, "-std=f2018", "-O2", "-fcheck=all",
                        "-Werror=implicit-interface", *map(str, SOURCES),
                        str(DRIVER), "-o", str(executable)], cwd=tmp, check=True)
        input_path = tmp / "input.bin"
        output_path = tmp / "output.bin"
        cases = []
        with input_path.open("wb") as writer:
            for case in range(480):
                nperseg = int(rng.integers(3, 22))
                noverlap = int(rng.integers(0, nperseg))
                nfft = nperseg + int(rng.integers(0, 10))
                one_sided = (case % 2 == 0)
                complex_output = not one_sided
                boundary = (case % 3 != 0)
                padded = (case % 5 != 0)
                scaling = ("spectrum", "psd")[case % 4 == 0]
                window_id = case % len(WINDOWS)
                window = window_spec(window_id, nperseg)
                fs = float(rng.uniform(0.5, 120))
                length = int(rng.integers(max(25, nperseg), 105))
                x = rng.normal(size=length)
                if complex_output:
                    x = x + 1j*rng.normal(size=length)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    _, _, zxx = signal.stft(
                        x, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap,
                        nfft=nfft, detrend=False, boundary="zeros" if boundary else None,
                        padded=padded, scaling=scaling,
                        return_onesided=one_sided)
                    _, reference = signal.istft(
                        zxx, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap,
                        nfft=nfft, input_onesided=one_sided, boundary=boundary,
                        scaling=scaling)
                header = [zxx.shape[0], zxx.shape[1], nperseg, noverlap, nfft,
                          int(one_sided), int(complex_output), int(boundary),
                          int(scaling == "psd"), window_id]
                writer.write(np.asarray(header, dtype=np.int32).tobytes())
                writer.write(struct.pack("d", fs))
                writer.write(np.asfortranarray(zxx, dtype=np.complex128).tobytes(order="F"))
                cases.append((reference, header))
        subprocess.run([str(executable), str(input_path), str(output_path)], check=True)
        raw = output_path.read_bytes()
        offset = 0
        for i, (reference, header) in enumerate(cases):
            status, nout = struct.unpack_from("ii", raw, offset)
            offset += 8
            got = np.frombuffer(raw, dtype=np.complex128, count=nout, offset=offset)
            offset += nout * 16
            if status not in (0, 2):
                raise AssertionError(f"case {i}: invalid output status={status}, header={header}")
            if status == 2:
                warning_cases += 1
            if got.shape != reference.shape:
                raise AssertionError(f"case {i}: shape {got.shape} != {reference.shape}; {header}")
            difference = float(np.max(np.abs(got - reference))) if nout else 0.0
            maximum = max(maximum, difference)
            if not np.allclose(got, reference, atol=4e-11, rtol=4e-11):
                raise AssertionError(f"case {i}: max error={difference:.3e}; {header}")
            count += 1
        if offset != len(raw):
            raise AssertionError("unexpected bytes remaining in Fortran output")
        print(f"SciPy ISTFT comparisons: {count} passed; max abs error {maximum:.3e}; "
              f"{warning_cases} incomplete-overlap statuses")


if __name__ == "__main__":
    main()
