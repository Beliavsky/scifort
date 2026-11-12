#!/usr/bin/env python3
"""Compare SciFort's spectral window arguments with SciPy 1.17.0.

Build the Fortran reference driver first (see README), then pass its executable
as an argument. Checks spectra, one-/two-sided frequency grids and frame times.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import warnings

import numpy as np
import scipy
from scipy import signal


def data_for(style: int, method: int):
    if style == 1:
        return ("tukey", 0.32)
    if style == 2:
        return ("kaiser", 6.2)
    if style == 3:
        return ("gaussian", 2.8)
    if style == 4:
        return ("general_gaussian", 1.4, 3.0)
    if style == 5:
        return ("general_cosine", [0.5, 0.4, 0.1])
    if style == 6:
        n = 42 if method == 6 else 10
        i = np.arange(1, n + 1)
        return 0.6 + 0.15 * np.sin(0.31 * i)
    if style == 7:
        return ("tukey", 0.25) if method == 2 else ("boxcar" if method == 6 else "hann")
    if style == 8:
        return "tukey"
    raise ValueError(style)


def parse_line(text: str, prefix: str):
    tokens = text.split()
    assert tokens[0] == prefix, (tokens[0], prefix)
    n = int(tokens[1])
    values = np.asarray([float(t) for t in tokens[2:]], dtype=np.float64)
    return n, values


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("executable", help="compiled Fortran reference driver")
    args = ap.parse_args()
    raw = subprocess.check_output([args.executable], text=True).splitlines()
    assert len(raw) == 96 * 4, f"unexpected output lines: {len(raw)}"
    i = np.arange(1, 43, dtype=np.float64)
    x = np.sin(0.37 * i) + 0.012 * i
    y = np.cos(0.11 * i) + 0.018 * i - 0.15 * x
    z = x + 1j * y
    q = y - 0.25j * x
    total = 0
    coefficients = 0
    largest = (0, None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        warnings.filterwarnings("ignore", message="Input data is complex, switching to return_onesided=False")
        for j in range(0, len(raw), 4):
            method, style, dtype = map(int, raw[j].split()[1:])
            nfreq, freqs = parse_line(raw[j + 1], "FREQ")
            nt, times = parse_line(raw[j + 2], "TIME")
            ndata, payload = parse_line(raw[j + 3], "DATA")
            assert len(freqs) == nfreq
            assert len(times) == nt
            assert len(payload) == 2 * ndata
            actual = payload[0::2] + 1j * payload[1::2]
            real = dtype == 1
            u = x if real else z
            v = y if real else q
            w = data_for(style, method)
            kw = dict(fs=7.25, window=w, nfft=13, nperseg=10,
                      noverlap=4, detrend="linear")
            if method == 1:
                f, t, output = signal.stft(u, boundary="odd", padded=True,
                                           scaling="psd", **kw)
            elif method == 2:
                f, t, output = signal.spectrogram(u, scaling="density", mode="psd", **kw)
            elif method == 3:
                f, output = signal.welch(u, scaling="density", **kw)
                t = np.empty(0)
            elif method == 4:
                f, output = signal.csd(u, v, scaling="density", **kw)
                t = np.empty(0)
            elif method == 5:
                f, output = signal.coherence(u, v, **kw)
                t = np.empty(0)
            elif method == 6:
                f, output = signal.periodogram(u, fs=7.25, window=w, nfft=45,
                                               detrend="linear", scaling="density")
                t = np.empty(0)
            else:
                raise AssertionError(method)
            expected = np.asarray(output).ravel(order="F")
            info = (method, style, dtype)
            if len(expected) != len(actual) or len(f) != nfreq or len(t) != nt:
                raise AssertionError(f"size mismatch at {info}: "
                    f"output expected {len(expected)} actual {len(actual)}, "
                    f"f {len(f)} vs {nfreq}, t {len(t)} vs {nt}")
            ferr = np.max(np.abs(f - freqs))
            terr = np.max(np.abs(t - times)) if nt else 0.0
            err = np.max(np.abs(actual - expected))
            total += 1
            coefficients += len(actual)
            if err > largest[0]:
                largest = (float(err), info)
            if ferr > 1e-11 or terr > 1e-11 or not np.allclose(
                    actual, expected, rtol=3e-11, atol=3e-12, equal_nan=True):
                print(f"FAILED case {info}: err={err:.5g}, "
                      f"ferr={ferr:.3g}, terr={terr:.3g}")
                k = np.argmax(np.abs(actual - expected))
                print(f"  coefficient {k}: Fortran={actual[k]}, SciPy={expected[k]}")
                sys.exit(1)
    print(f"SciPy {scipy.__version__}: {total} parameterized-window cases, "
          f"{coefficients} coefficients PASSED; maximum absolute error "
          f"{largest[0]:.3e} in (method, window style, input type)={largest[1]}")


if __name__ == "__main__":
    main()
