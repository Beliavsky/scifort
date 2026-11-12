#!/usr/bin/env python3
"""Independently compare legacy STFT boundary modes with SciPy 1.17.

Usage: python test/compare_signal_stft_boundaries_scipy.py /path/to/driver
The driver is the separate, dependency-free Fortran program in this directory.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy import signal


def run(path: str) -> None:
    completed = subprocess.run([path], check=True, capture_output=True, text=True)
    lines = iter(completed.stdout.splitlines())
    cases = 0
    coefficients = 0
    worst = 0.0
    for line in lines:
        fields = line.split()
        if fields[0] != "CASE":
            raise ValueError(f"Unexpected Fortran output: {line!r}")
        it, n, nw, no, nf, padded, nfreq, nt = map(int, fields[1:9])
        mode, scaling, detrend, window = fields[9:13]
        freqs = np.fromstring(next(lines), sep=" ")
        times = np.fromstring(next(lines), sep=" ")
        if len(freqs) != nfreq or len(times) != nt:
            raise AssertionError(f"Reference-coordinate shape mismatch in case {cases}")
        raw = np.array([
            np.fromstring(next(lines), sep=" ") for _ in range(nt)
        ])
        if raw.shape != (nt, 2 * nfreq):
            raise AssertionError(f"Invalid driver spectrum shape {raw.shape}")
        computed = (raw[:, ::2] + 1j * raw[:, 1::2]).T
        i = np.arange(1, n + 1, dtype=np.float64)
        x = 0.8 * np.cos(0.43 * i) + 0.07 * i - 0.02 * np.sin(0.31 * i * i)
        if it == 2:
            x = x + 1j * (0.23 * np.sin(0.21 * i) - 0.03 * i)
        expected_f, expected_t, expected = signal.stft(
            x, fs=5.5, window=window, nperseg=nw, noverlap=no,
            nfft=nf, detrend=False if detrend == "none" else detrend,
            boundary=None if mode == "none" else mode, padded=bool(padded),
            return_onesided=it == 1, scaling=scaling,
        )
        tag = f"case={cases} kind={it} n={n} nw={nw} no={no} nf={nf} {mode} {scaling} {detrend} {window}"
        if expected.shape != computed.shape:
            raise AssertionError(f"{tag}: {computed.shape} != {expected.shape}")
        np.testing.assert_allclose(freqs, expected_f, atol=2e-13, rtol=2e-13, err_msg=tag)
        np.testing.assert_allclose(times, expected_t, atol=2e-13, rtol=2e-13, err_msg=tag)
        np.testing.assert_allclose(computed, expected, atol=5e-12, rtol=5e-12, err_msg=tag)
        worst = max(worst, float(np.max(np.abs(computed - expected))))
        coefficients += expected.size
        cases += 1
    print(f"STFT SciPy comparisons passed: {cases} cases, {coefficients} coefficients; "
          f"max absolute difference {worst:.3g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("driver", help="compiled stft_boundaries_driver executable")
    args = parser.parse_args()
    run(args.driver)
