#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent SciPy reference comparisons for streaming STFT/spectrogram.

Usage: python test/compare_signal_stft_stream_scipy.py PATH_TO_FORTRAN_EXECUTABLE
"""
from __future__ import annotations

import itertools
from pathlib import Path
import subprocess
import sys
import warnings

import numpy as np
from scipy import signal


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("pass compiled reference_signal_stft_stream executable")
    executable = Path(sys.argv[1]).resolve()
    if not executable.is_file():
        raise SystemExit(f"executable not found: {executable}")
    rng = np.random.default_rng(31971)
    combinations = itertools.product(
        (0, 1),  # stft vs spectrogram
        (0, 1),  # real vs complex
        (3, 5, 8, 11),  # window length
        (0, 1),  # nonzero overlap or none
        (0, 1),  # equal or longer nfft
        ("boxcar", "hann", "default", "tukey"),
        ("none", "constant", "linear"),
    )
    total_cases = total_values = 0
    maximum_error = 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        for mode, is_complex, nw, overlap_flag, pad_flag, window, trend in combinations:
            no = (nw - 1) if overlap_flag else nw // 3
            nf = nw + 2 if pad_flag else nw
            n = 4 * nw + 7 + (nw % 3)
            fs = 16.25
            signal_data = rng.standard_normal(n) + 0.2 * np.sin(np.arange(n) * 0.34)
            imag = rng.standard_normal(n) * 0.35
            if is_complex:
                signal_data = signal_data + 1j * imag
            edge = ("none", "zeros")[total_cases % 2]
            padded = total_cases % 3 != 0
            if mode == 0:
                scaling = ("spectrum", "psd")[total_cases % 4 > 1]
            else:
                scaling = ("density", "spectrum")[total_cases % 4 > 1]
                edge = "none"
                padded = False
            x = signal_data.astype(np.complex128)
            rows = [
                f"{mode} {n} {nw} {no} {nf} {is_complex} {int(padded)} {fs}",
                f"{edge} {scaling} {trend} {window}",
            ]
            rows.extend(f"{val.real:.17g} {val.imag:.17g}" for val in x)
            proc = subprocess.run(
                [str(executable)], input="\n".join(rows) + "\n", text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=8,
            )
            if proc.returncode != 0:
                raise AssertionError(f"Fortran failed case {total_cases}: {proc.stderr}")
            lines = proc.stdout.strip().splitlines()
            nfreq, nframes = map(int, lines[0].split())
            f = np.array([float(item) for item in lines[1:1+nfreq]])
            t = np.array([float(item) for item in lines[1+nfreq:1+nfreq+nframes]])
            entries = [tuple(map(float, item.split())) for item in lines[1+nfreq+nframes:]]
            if len(entries) != nfreq * nframes:
                raise AssertionError((total_cases, "wrong output length", len(entries)))
            y = np.array([complex(re, im) for re, im in entries]).reshape(nframes, nfreq).T
            scipy_trend = False if trend == "none" else trend
            window_options = {} if window == "default" else {
                "window": ("tukey", 0.5) if window == "tukey" else window,
            }
            if mode == 0:
                f_ref, t_ref, y_ref = signal.stft(
                    signal_data, fs=fs, nperseg=nw, noverlap=no,
                    nfft=nf, detrend=scipy_trend, boundary=None if edge == "none" else edge,
                    padded=padded, return_onesided=not is_complex, scaling=scaling,
                    **window_options,
                )
            else:
                f_ref, t_ref, y_ref = signal.spectrogram(
                    signal_data, fs=fs, nperseg=nw, noverlap=no,
                    nfft=nf, detrend=scipy_trend, return_onesided=not is_complex,
                    scaling=scaling, mode="psd", **window_options,
                )
            if y.shape != y_ref.shape:
                raise AssertionError((total_cases, "shape", y.shape, y_ref.shape))
            if not np.allclose(f, f_ref, atol=1e-12, rtol=0):
                raise AssertionError((total_cases, "frequencies", f, f_ref))
            if not np.allclose(t, t_ref, atol=1e-12, rtol=0):
                raise AssertionError((total_cases, "times", t, t_ref))
            error = float(np.max(np.abs(y - y_ref)))
            if error > 1e-10:
                raise AssertionError((total_cases, mode, nw, no, nf, is_complex, scaling, error))
            maximum_error = max(maximum_error, error)
            total_cases += 1
            total_values += y.size
    print(f"PASS: {total_cases} SciPy cases, {total_values} spectral bins, "
          f"max absolute error {maximum_error:.3g}")


if __name__ == "__main__":
    main()
