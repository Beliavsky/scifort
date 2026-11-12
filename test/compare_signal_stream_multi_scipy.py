#!/usr/bin/env python3
"""Compare SciFort multichannel stream output independently with SciPy 1-D channels.

Usage: python test/compare_signal_stream_multi_scipy.py PATH_TO_REFERENCE_EXECUTABLE
"""
from __future__ import annotations

import itertools
from pathlib import Path
import subprocess
import sys
import warnings

import numpy as np
from scipy import signal


def run(executable: Path, header: str, config: str, data: np.ndarray) -> list[str]:
    rows = [header, config]
    rows.extend(f"{z.real:.17g} {z.imag:.17g}" for z in data)
    process = subprocess.run(
        [str(executable)], input="\n".join(rows) + "\n",
        capture_output=True, text=True, timeout=12,
    )
    if process.returncode:
        raise AssertionError(f"Fortran driver failed: {header}: {process.stderr}")
    return process.stdout.strip().splitlines()


def compare_forward(executable: Path):
    rng = np.random.default_rng(122914)
    cases = bins = 0
    worst = 0.0
    for mode, nch, cmplx, nw, overlap, pad_fft, window in itertools.product(
        (0, 1), (1, 3), (0, 1), (5, 8), (0, 1), (0, 1), ("boxcar", "hann")
    ):
        nf = nw + (2 if pad_fft else 0)
        no = nw // 2 if overlap else 0
        n = 3 * nw + 17
        fs = 13.75
        x = rng.normal(size=(n, nch)) + 0.3 * np.sin(np.arange(n)[:, None] * .13)
        imag = 0.4 * rng.normal(size=(n, nch))
        if cmplx:
            x = x + 1j * imag
        boundary = "zeros" if cases % 2 == 0 else "none"
        padded = cases % 3 != 0
        scaling = ("spectrum", "psd")[cases % 4 > 1]
        if mode == 1:
            boundary = "none"
            padded = False
            scaling = ("density", "spectrum")[cases % 4 > 1]
        detrend = ("none", "constant", "linear")[cases % 3]
        header = f"{mode} {nch} {n} {nw} {no} {nf} {cmplx} {int(padded)} {fs}"
        config = f"{boundary} {scaling} {detrend} {window}"
        lines = run(executable, header, config, x.T.reshape(-1))
        nfreq, nframes, channels = map(int, lines[0].split())
        assert channels == nch
        f = np.array([float(v) for v in lines[1:1 + nfreq]])
        t = np.array([float(v) for v in lines[1 + nfreq:1 + nfreq + nframes]])
        v = np.array([complex(*map(float, line.split())) for line in lines[1 + nfreq + nframes:]])
        if v.size != nfreq * nch * nframes:
            raise AssertionError((cases, "output size", v.size))
        y = v.reshape(nframes, nch, nfreq).transpose(2, 1, 0)
        trend = False if detrend == "none" else detrend
        ref_signals = x.T.copy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=DeprecationWarning)
            if mode == 0:
                rf, rt, ry = signal.stft(
                    ref_signals, fs=fs, window=window, nperseg=nw, noverlap=no,
                    nfft=nf, detrend=trend, boundary=None if boundary == "none" else boundary,
                    padded=padded, scaling=scaling, return_onesided=not cmplx,
                )
            else:
                rf, rt, ry = signal.spectrogram(
                    ref_signals, fs=fs, window=window, nperseg=nw, noverlap=no,
                    nfft=nf, detrend=trend, scaling=scaling, return_onesided=not cmplx,
                    mode="psd",
                )
        ry = ry.transpose(1, 0, 2)
        if y.shape != ry.shape:
            raise AssertionError((cases, "shape", y.shape, ry.shape))
        if not np.allclose(f, rf, rtol=0, atol=1e-12):
            raise AssertionError((cases, "frequencies", f, rf))
        if not np.allclose(t, rt, rtol=0, atol=1e-12):
            raise AssertionError((cases, "times", t, rt))
        err = float(np.max(np.abs(y - ry)))
        if err > 1e-10:
            raise AssertionError((cases, mode, window, nw, no, nf, cmplx, err))
        worst = max(worst, err)
        cases += 1
        bins += y.size
    return cases, bins, worst


def compare_inverse(executable: Path):
    rng = np.random.default_rng(71781)
    cases = values = 0
    worst = 0.0
    for nch, onesided, nw, pad_fft, window, boundary, scaling in itertools.product(
        (1, 3), (0, 1), (7, 8), (0, 1), ("boxcar", "hann"), (False, True),
        ("spectrum", "psd")
    ):
        nf = nw + (2 if pad_fft else 0)
        nfreq = nf // 2 + 1 if onesided else nf
        no = nw // 2 if cases % 2 == 0 else nw - 1
        nframes = 5 + cases % 4
        fs = 15.125
        z = rng.normal(size=(nfreq, nch, nframes)) + .5j * rng.normal(
            size=(nfreq, nch, nframes)
        )
        header = f"{2 if onesided else 3} {nch} {nframes} {nw} {no} {nfreq} {nf} {int(boundary)} {fs}"
        lines = run(executable, header, f"none {scaling} none {window}", z.transpose(2, 1, 0).reshape(-1))
        nsamples, channels, third = map(int, lines[0].split())
        assert channels == nch and third == 0
        t = np.array([float(v) for v in lines[1:1 + nsamples]])
        v = np.array([complex(*map(float, line.split())) for line in lines[1 + nsamples:]])
        if v.size != nsamples * nch:
            raise AssertionError((cases, "inverse output size", v.size, nsamples * nch))
        x = v.reshape(nch, nsamples).T
        reference = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for j in range(nch):
                tref, xref = signal.istft(
                    z[:, j, :], fs=fs, window=window, nperseg=nw,
                    noverlap=no, nfft=nf, input_onesided=bool(onesided),
                    boundary=boundary, scaling=scaling,
                )
                if not np.allclose(t, tref, rtol=0, atol=1e-12):
                    raise AssertionError((cases, "inverse times", t, tref))
                reference.append(xref)
        ref = np.stack(reference, axis=1)
        if x.shape != ref.shape:
            raise AssertionError((cases, "inverse shape", x.shape, ref.shape))
        err = float(np.max(np.abs(x - ref)))
        if err > 1e-10:
            raise AssertionError((cases, "inverse", onesided, nw, nf, no, err))
        worst = max(worst, err)
        cases += 1
        values += x.size
    return cases, values, worst


def main():
    if len(sys.argv) != 2:
        raise SystemExit("provide path to reference_signal_stream_multi executable")
    executable = Path(sys.argv[1]).resolve()
    if not executable.is_file():
        raise SystemExit(f"executable not found: {executable}")
    fcases, bins, ferr = compare_forward(executable)
    icases, values, ierr = compare_inverse(executable)
    print(f"PASS: {fcases} SciPy forward/spectrogram cases, {bins} spectral bins, "
          f"max error {ferr:.3g}")
    print(f"PASS: {icases} SciPy inverse cases, {values} samples, max error {ierr:.3g}")


if __name__ == "__main__":
    main()
