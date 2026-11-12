#!/usr/bin/env python3
"""SciPy comparisons for the cost-aware FFT-convolution padding selector."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.signal import fftconvolve, oaconvolve

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checked", action="store_true")
    args = parser.parse_args()
    rng = np.random.default_rng(20261026)
    a = np.asfortranarray(rng.standard_normal((101, 251)))
    b = np.asfortranarray(rng.standard_normal((29, 59)))
    ca = np.asfortranarray(a + 1j * rng.standard_normal(a.shape))
    cb = np.asfortranarray(b + 1j * rng.standard_normal(b.shape))
    sa = a.reshape((101,) + (1,) * 13 + (251,), order="F")
    sb = b.reshape((29,) + (1,) * 13 + (59,), order="F")
    x = rng.standard_normal(9001)
    h = rng.standard_normal(103)
    da = np.asfortranarray(rng.standard_normal((513, 257)))
    db = np.asfortranarray(rng.standard_normal((11, 7)))
    arrays = [a, b, ca, cb, sa, sb, x, h, da, db]
    expected = [
        fftconvolve(a, b, mode="full"),
        fftconvolve(a, b, mode="same"),
        fftconvolve(a, b, mode="valid"),
        fftconvolve(ca, cb),
        fftconvolve(sa, sb, axes=(0, 14)),
        oaconvolve(x, h),
        oaconvolve(da, db),
    ]
    sources = ["kinds", "signal_dpss", "signal_windows_extra", "signal_spectral", "signal_fftconvolve"]
    with tempfile.TemporaryDirectory(prefix="scifort_fft_padding_") as temp:
        folder = Path(temp)
        flags = ["-std=f2018", "-Werror=implicit-interface", "-J", str(folder), "-I", str(folder)]
        flags += ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if args.checked else ["-O3"]
        cmd = ["gfortran", *flags, *(str(ROOT / "src" / f"{name}.f90") for name in sources),
               str(ROOT / "validation/reference_drivers/signal_fft_padding_reference_driver.f90"), "-o", str(folder / "driver")]
        subprocess.run(cmd, check=True, cwd=folder, capture_output=True, text=True)
        with (folder / "input.bin").open("wb") as stream:
            for array in arrays:
                stream.write(np.asarray(array).tobytes(order="F"))
        subprocess.run([str(folder / "driver"), str(folder / "input.bin"), str(folder / "output.bin")],
                       check=True, cwd=folder, capture_output=True, text=True, timeout=120)
        raw = (folder / "output.bin").read_bytes()
        pos = 0
        maximum = 0.0
        count = 0
        for index, want in enumerate(expected):
            itemsize = np.asarray(want).dtype.itemsize
            needed = want.size * itemsize
            got = np.frombuffer(raw, dtype=want.dtype, offset=pos, count=want.size).reshape(
                want.shape, order="F")
            pos += needed
            error = float(np.max(np.abs(got - want)))
            maximum = max(maximum, error)
            count += want.size
            if not np.allclose(got, want, atol=3e-10, rtol=3e-11):
                raise AssertionError(f"case {index}: largest discrepancy {error:.3g}")
        if pos != len(raw):
            raise AssertionError("Unexpected trailing output from Fortran reference driver")
        mode = "checked" if args.checked else "optimized"
        print(f"PASS {len(expected)} SciPy padding cases ({mode}); "
              f"{count} outputs; maximum absolute error {maximum:.3e}")


if __name__ == "__main__":
    main()
