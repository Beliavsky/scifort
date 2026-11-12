#!/usr/bin/env python3
"""Compare equal-dual window construction with SciPy.

Compile validation/reference_drivers/reference_signal_short_time_fft_tight.f90, then pass executable path.
"""
import pathlib
import subprocess
import sys

import numpy as np
from scipy.signal import ShortTimeFFT


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: compare_signal_short_time_fft_tight.py EXE")
    run = subprocess.run([str(pathlib.Path(sys.argv[1]).resolve())],
                         check=True, capture_output=True, text=True)
    output = iter(run.stdout.splitlines())
    number = 0
    worst = 0.0
    while (header := next(output, None)) is not None:
        n, hop, nfft, _, status, scale = header.split()
        n, hop, nfft = int(n), int(hop), int(nfft)
        assert status == "0"
        got = np.array([list(map(float, next(output).split())) for _ in range(n)])
        desired = np.exp(-((np.arange(n) - 0.5 * (n - 1)) / 2) ** 2)
        ref = ShortTimeFFT.from_win_equals_dual(
            desired, hop, fs=8, fft_mode="twosided", mfft=nfft,
            scale_to=None if scale == "none" else scale,
        )
        err = max(np.max(np.abs(got[:, 0] + 1j * got[:, 1] - ref.win)),
                  np.max(np.abs(got[:, 2] + 1j * got[:, 3] - ref.dual_win)))
        assert err < 5.0e-12, (n, hop, nfft, scale, err)
        worst = max(worst, err)
        number += 1
    assert number == 96, number
    print(f"Equal-dual/tight windows: {number} SciPy cases passed, max difference {worst:.4g}")


if __name__ == "__main__":
    main()
