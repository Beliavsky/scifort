#!/usr/bin/env python3
"""Reference comparison of SciFort ShortTimeFFT with scipy.signal.ShortTimeFFT.

Usage: python test/compare_signal_short_time_fft.py /path/to/fortran_reference_executable
Compile validation/reference_drivers/reference_signal_short_time_fft.f90 for the executable.
The Fortran program reads short_fft_cases.txt and writes short_fft_results.txt.
"""
import itertools
import pathlib
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import ShortTimeFFT, get_window


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: compare_signal_short_time_fft.py EXECUTABLE")
    exe = str(pathlib.Path(sys.argv[1]).resolve())
    cases = []
    for k, (nw, hop, nf, mode, scale, pad, phase) in enumerate(itertools.product(
        [5, 6, 7, 8, 9, 12],
        [2, 3],
        [8, 13, 16],
        ["onesided", "onesided2X", "twosided", "centered"],
        ["none", "magnitude", "psd"],
        ["zeros", "edge", "even", "odd"],
        [0, 2],
    )):
        if nf < nw or (mode == "onesided2X" and scale == "none"):
            continue
        # Deliberately exercise odd/even window lengths and FFT sizes, each padding mode,
        # all scaling modes, and real or complex STFT inputs in a deterministic spread.
        if (k * 37 + nw) % 9:
            continue
        n = 19 + (k % 23)
        kind = "complex" if mode in ("twosided", "centered") and k % 3 == 0 else "real"
        fs = 16.0
        sft = ShortTimeFFT(
            get_window("hann", nw), hop, fs, mfft=nf, fft_mode=mode,
            scale_to=None if scale == "none" else scale, phase_shift=phase,
        )
        start = sft.p_min + (k % 2)
        end = sft.p_max(n) - (k % 3)
        if start >= end:
            continue
        offset = (k % 5) - 2
        trend = ["none", "constant", "linear"][k % 3]
        cases.append((n, nw, hop, nf, fs, mode, scale, phase, pad, kind,
                      start, end, offset, trend))
    assert len(cases) >= 50, len(cases)
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        (work / "short_fft_cases.txt").write_text(
            str(len(cases)) + "\n" + "\n".join(" ".join(map(str, row)) for row in cases) + "\n"
        )
        run = subprocess.run([exe], cwd=work, capture_output=True, text=True)
        if run.returncode: raise RuntimeError(f"Fortran reference exit {run.returncode}:\n{run.stdout}\n{run.stderr}")
        lines = iter((work / "short_fft_results.txt").read_text().splitlines())
        largest = 0.0
        nvalues = 0
        for ic, (n, nw, hop, nf, fs, mode, scale, phase, pad, kind, start, end, offset, trend) in enumerate(cases, 1):
            seq = list(map(int, next(lines).split()))
            if seq[0] != ic:
                raise AssertionError(f"Invalid case index: {seq[0]} != {ic}")
            sizef, sizet, p0 = seq[1:]
            fvalues = np.array([
                complex(*map(float, next(lines).split())) for _ in range(sizef * sizet)
            ]).reshape((sizet, sizef)).T
            x = np.cos(0.17 * np.arange(n)) + 0.3 * np.sin(0.32 * np.arange(n))
            if kind == "complex":
                x = x + 0.37j * np.cos(0.09 * np.arange(n))
            sft = ShortTimeFFT(
                get_window("hann", nw), hop, fs,
                mfft=nf, fft_mode=mode, scale_to=None if scale == "none" else scale,
                phase_shift=phase,
            )
            if trend == "none":
                ref = sft.stft(x, p0=start, p1=end, k_offset=offset, padding=pad)
            else:
                ref = sft.stft_detrend(x, trend, p0=start, p1=end, k_offset=offset, padding=pad)
            if ref.shape != (sizef, sizet) or start != p0:
                raise AssertionError(f"Case {ic}: shape/p_min mismatch: "
                                     f"{ref.shape} != {(sizef,sizet)}; {start} != {p0}")
            difference = np.abs(ref - fvalues)
            err = np.max(difference)
            if err > 3.0e-11 * max(1.0, np.max(np.abs(ref))):
                q, p = np.unravel_index(np.argmax(difference), difference.shape)
                raise AssertionError(f"Case {ic}: {mode=} {pad=} {nw=} {hop=} {nf=} {scale=} "
                                     f"{phase=} {kind=} {trend=} {offset=} {start=} {end=} err={err:g}, bin={q}, slice={p}, "
                                     f"expected={ref[q,p]}, got={fvalues[q,p]}")
            largest = max(largest, err)
            nvalues += ref.size
        assert next(lines, None) is None
    print(f"ShortTimeFFT vs SciPy: {len(cases)} cases, {nvalues} complex values, "
          f"maximum absolute difference {largest:.4g}")


if __name__ == "__main__":
    main()
