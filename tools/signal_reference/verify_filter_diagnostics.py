#!/usr/bin/env python3
"""Check SciFort lfiltic and group_delay against SciPy (600 seeded cases).

Requires NumPy, SciPy, and gfortran; no runtime SciPy dependency in SciFort.
    python tools/signal_reference/verify_filter_diagnostics.py
"""
from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261016


def real_line(values: np.ndarray) -> str:
    return " ".join(f"{float(v):.17e}" for v in values)


def complex_line(values: np.ndarray) -> str:
    return " ".join(f"({z.real:.17e},{z.imag:.17e})" for z in values)


def main() -> None:
    rng = np.random.default_rng(SEED)
    cases = []
    lines = []
    for op in (1, 2, 3, 4):
        for _ in range(150):
            is_complex = op in (2, 4)
            nb, na = rng.integers(1, 7, size=2)
            b = np.r_[rng.uniform(.7, 2.0), rng.normal(scale=.16, size=nb-1)]
            a = np.r_[rng.uniform(.7, 2.0), rng.normal(scale=.16, size=na-1)]
            if is_complex:
                b = b + 1j*rng.normal(scale=.15, size=nb)
                a = a + 1j*rng.normal(scale=.15, size=na)
            writer = complex_line if is_complex else real_line
            fs = float(rng.choice((2.0, 2.0*np.pi, 200.0)))
            if op in (1, 2):
                ny = int(rng.integers(0, 9))
                nx = int(rng.integers(-1, 9))
                y = rng.normal(size=ny)
                x = None if nx == -1 else rng.normal(size=nx)
                if is_complex:
                    y = y + 1j*rng.normal(size=ny)
                    if x is not None:
                        x = x + 1j*rng.normal(size=nx)
                expected = signal.lfiltic(b, a, y, x)
                lines.append(f"{op} {nb} {na} {ny} {nx} 0 {fs:.17e}")
                lines.extend((writer(b), writer(a)))
                if ny > 0:
                    lines.append(writer(y))
                if x is not None and nx > 0:
                    lines.append(writer(x))
            else:
                nw = int(rng.integers(2, 18))
                w = rng.uniform(-.3*fs, .8*fs, size=nw)
                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    _, expected = signal.group_delay((b, a), w=w, fs=fs)
                lines.append(f"{op} {nb} {na} 0 0 {nw} {fs:.17e}")
                lines.extend((writer(b), writer(a), real_line(w)))
            cases.append((op, expected))

    with tempfile.TemporaryDirectory(prefix="scifort-filter-diagnostics-") as tmp:
        directory = Path(tmp)
        executable = directory / "filter_diagnostics_driver"
        command = [
            "gfortran", "-std=f2018", "-O2", "-Werror=implicit-interface",
            "-J", str(directory), "-I", str(directory),
            str(ROOT / "src/kinds.f90"),
            str(ROOT / "src/signal_lfilter.f90"),
            str(ROOT / "src/signal_spectral.f90"),
            str(ROOT / "src/signal_analog_design.f90"),
            str(ROOT / "src/signal_design.f90"),
            str(ROOT / "tools/signal_reference/filter_diagnostics_driver.f90"),
            "-o", str(executable),
        ]
        subprocess.run(command, cwd=directory, check=True)
        result = subprocess.run([str(executable)], input="\n".join(lines)+"\n",
                                text=True, capture_output=True, check=True)

    output = iter(result.stdout.splitlines())
    largest_error = 0.0
    for count, (op, expected) in enumerate(cases, start=1):
        status, length = map(int, next(output).split())
        if status != 0 or length != len(expected):
            raise AssertionError(f"case {count}: unexpected status/size {status}/{length}")
        if op == 2:
            received = np.asarray([complex(*map(float, next(output).split()))
                                   for _ in range(length)])
        elif length == 0:
            received = np.array([])
        else:
            received = np.fromstring(next(output), sep=" ")
        if received.shape != expected.shape:
            raise AssertionError(f"case {count}: output shape differs")
        if length:
            difference = np.max(np.abs(received - expected))
            largest_error = max(largest_error, float(difference))
            if not np.allclose(received, expected, rtol=2e-10, atol=3e-11):
                raise AssertionError(f"case {count}: {difference:g} error, op={op}")
    print(f"Verified {len(cases)} cases against SciPy {signal.__version__ if hasattr(signal, '__version__') else ''}; "
          f"maximum absolute error: {largest_error:.3g}")


if __name__ == "__main__":
    main()
