#!/usr/bin/env python3
"""Compare Fourier-bin weighted SciFort resampling with SciPy 1.17.

Run: python test/compare_signal_resample_frequency_window_scipy.py
The standalone compiler dependency walk is local; no fpm or network required.
Set FC for a different Fortran compiler, or SCIFORT_OPT=-O3 for optimization.
"""
from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
DRIVER = ROOT / "tools/signal_reference/resample_frequency_window_driver.f90"


def compile_driver(build: Path) -> Path:
    modules = {}
    for source in SRC.rglob("*.f90"):
        for name in re.findall(r"^\s*module\s+(?!procedure|function|subroutine)([a-z]\w*)",
                               source.read_text(), re.I | re.M):
            modules[name.lower()] = source
    visited, visiting, order = set(), set(), []

    def visit(name: str):
        if name in visited:
            return
        if name in visiting:
            raise ValueError(f"cycle at {name}")
        visiting.add(name)
        source = modules[name]
        for dependency in re.findall(r"^\s*use\s+(scifort_\w+)",
                                     source.read_text(), re.I | re.M):
            visit(dependency.lower())
        visiting.remove(name)
        visited.add(name)
        order.append(source)

    visit("scifort_signal_multirate_nd_into")
    fc = os.environ.get("FC", "gfortran")
    flags = ["-std=f2018", "-Werror=implicit-interface", "-I", str(build),
             "-J", str(build), os.environ.get("SCIFORT_OPT", "-O0"), "-fcheck=all"]
    objects = []
    for source in order:
        obj = build / (source.stem + ".o")
        subprocess.run([fc, *flags, "-c", str(source), "-o", str(obj)],
                       cwd=build, check=True, capture_output=True, text=True)
        objects.append(str(obj))
    exe = build / "resample_frequency_window_driver"
    subprocess.run([fc, *flags, str(DRIVER), *objects, "-o", str(exe)],
                   cwd=build, check=True, capture_output=True, text=True)
    return exe


def make_cases():
    rng = np.random.default_rng(20261026)
    cases = []
    # Include the one-point case, odd/even lengths, upsampling, downsampling,
    # general asymmetric real windows, complex windows, and frequency-domain input.
    for length in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 19, 20, 23, 29]:
        for kind, weights_kind, from_fft in [(0, 0, 0), (1, 0, 0), (1, 1, 0),
                                             (1, 0, 1), (1, 1, 1),
                                             (0, 2, 0), (1, 2, 0), (1, 2, 1),
                                             (0, 3, 0), (1, 3, 1)]:
            for repeat in range(5):
                channels = int(rng.integers(1, 5))
                count = int(rng.integers(1, 37))
                if repeat == 0:
                    count = length
                elif repeat == 1:
                    count = max(1, length // 2)
                elif repeat == 2:
                    count = length + 5
                x = rng.normal(size=(channels, length))
                if kind:
                    x = x + 1j*rng.normal(size=x.shape)
                weights = rng.normal(loc=0.6, scale=0.5, size=length)
                if weights_kind == 1:
                    weights = weights + 1j*rng.normal(scale=0.2, size=length)
                transformed = np.fft.fft(x, axis=1) if from_fft else x
                cases.append((length, count, channels, kind, weights_kind, from_fft,
                              transformed, weights))
    return cases


def main():
    cases = make_cases()
    data = []
    expected = []
    for n, count, channels, kind, weights_kind, from_fft, x, weights in cases:
        data.append(f"{n} {count} {channels} {kind} {weights_kind} {from_fft}\n")
        for weight in weights:
            data.append(f"{weight.real:.17g} {weight.imag:.17g}\n")
        for entry in x.ravel(order="F"):
            data.append(f"{entry.real:.17g} {entry.imag:.17g}\n")
        window = weights if weights_kind in (0, 1) else "hann" if weights_kind == 2 else None
        expected.append(signal.resample(x, count, axis=1, window=window,
                                        domain="freq" if from_fft else "time"))
    with tempfile.TemporaryDirectory(prefix="scifort_resample_weights_") as directory:
        driver = compile_driver(Path(directory))
        run = subprocess.run([str(driver)], input="".join(data), check=True,
                             capture_output=True, text=True)
        lines = iter(run.stdout.splitlines())
        largest = 0.0
        comparisons = 0
        for case_number, ref in enumerate(expected):
            st, channels, count = map(int, next(lines).split())
            if st != 0 or ref.shape != (channels, count):
                raise AssertionError(f"invalid output at case {case_number}: {st} {channels} {count}")
            received = np.array([complex(*map(float, next(lines).split()))
                                 for _ in range(channels*count)]).reshape((channels,count), order="F")
            error = float(np.max(np.abs(received-ref)))
            largest = max(largest, error)
            np.testing.assert_allclose(received, ref, rtol=6e-12, atol=6e-12,
                                       err_msg=f"case {case_number}: {cases[case_number][:6]}")
            comparisons += ref.size
    print(f"PASS: {len(cases)} SciPy resample array-window cases, {comparisons} values, "
          f"largest absolute difference {largest:.3e}")


if __name__ == "__main__":
    main()
