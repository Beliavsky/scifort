#!/usr/bin/env python3
"""Numerically compare rank-independent SciFort inverse STFT to SciPy 1.17.

Usage: python test/compare_signal_istft_nd_scipy.py /path/to/reference_signal_istft_nd
"""
import subprocess
import sys
import warnings

import numpy as np
from scipy.signal import istft


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_signal_istft_nd_scipy.py DRIVER")
    rng = np.random.default_rng(20261025)
    executable = sys.argv[1]
    tested, coefficients, worst = 0, 0, 0.0
    for axis in (0, 1):
        for onesided in (True, False):
            for real_output in (True, False):
                if real_output and not onesided:
                    continue
                for nw, overlap, nfft in ((4, 2, 4), (6, 3, 8),
                                          (7, 4, 7), (8, 6, 12), (10, 5, 10)):
                    for window in ("hann", "hamming", "boxcar"):
                        for scaling in ("spectrum", "psd"):
                            for boundary in (True, False):
                                # Unextended Hann or Hamming may have zero
                                # denominator at the first/last sample; test
                                # boxcar only when boundary=False.
                                if not boundary and window != "boxcar":
                                    continue
                                fs = 9.5
                                nfreq = nfft // 2 + 1 if onesided else nfft
                                nframes = 5
                                if axis == 0:
                                    shape = (nfreq, 3, nframes)
                                else:
                                    shape = (3, nfreq, nframes)
                                z = rng.standard_normal(shape) + 1j*rng.standard_normal(shape)
                                nout = nw + (nframes-1)*(nw-overlap)
                                if boundary:
                                    nout -= 2*(nw//2)
                                header = (f"{shape[0]} {shape[1]} {shape[2]} {nout} "
                                          f"{axis+1} {nw} {overlap} {nfft} "
                                          f"{int(real_output)} {int(boundary)} "
                                          f"{int(onesided)} {fs:.17g}\n{window} {scaling}\n")
                                flat = z.ravel(order="F")
                                payload = header + "".join(
                                    f"{c.real:.17g} {c.imag:.17g}\n" for c in flat)
                                proc = subprocess.run([executable], input=payload,
                                                      text=True, capture_output=True,
                                                      check=True)
                                lines = proc.stdout.splitlines()
                                if not lines or int(lines[0]) != 0:
                                    raise AssertionError((axis, shape, lines, proc.stderr))
                                values = np.array([[float(t) for t in line.split()]
                                                   for line in lines[1:]], dtype=float)
                                if values.shape != (3*nout, 2):
                                    raise AssertionError((values.shape, shape, nout))
                                y = values[:, 0]+1j*values[:, 1]
                                if real_output:
                                    y = y.real
                                if axis == 0:
                                    y = y.reshape((nout, 3), order="F")
                                else:
                                    y = y.reshape((3, nout), order="F")
                                with warnings.catch_warnings():
                                    warnings.simplefilter("ignore", UserWarning)
                                    _, expected = istft(z, fs=fs, window=window,
                                                        nperseg=nw, noverlap=overlap,
                                                        nfft=nfft, input_onesided=onesided,
                                                        boundary=boundary, scaling=scaling,
                                                        freq_axis=axis, time_axis=2)
                                # SciPy puts reconstructed time last; SciFort
                                # replaces the selected frequency axis.
                                expected = np.moveaxis(expected, -1, axis)
                                if real_output:
                                    expected = expected.real
                                if y.shape != expected.shape:
                                    raise AssertionError((y.shape, expected.shape))
                                error = float(np.max(np.abs(y-expected)))
                                worst = max(worst, error)
                                if not np.allclose(y, expected, rtol=2e-11, atol=4e-12):
                                    raise AssertionError((axis, shape, window, overlap,
                                                          scaling, boundary, error))
                                tested += 1
                                coefficients += y.size
    print(f"PASS: {tested} SciPy 1.17 inverse STFT cases, "
          f"{coefficients} samples, max absolute error {worst:.3g}")


if __name__ == "__main__":
    main()
