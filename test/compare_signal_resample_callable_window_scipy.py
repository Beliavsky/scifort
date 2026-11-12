#!/usr/bin/env python3
"""Cross-check callable Fourier-bin resampling with SciPy.

Invoke from any directory: python test/compare_signal_resample_callable_window_scipy.py
The default performs its own standalone Fortran module build.
Set SCIFORT_OPT=-O3 to enable optimization in that build. Alternatively,
SCIFORT_REFERENCE_EXE can specify an already built reference driver.
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
DRIVER = ROOT / "tools/signal_reference/resample_callable_window_driver.f90"


def build_driver(build: Path) -> Path:
    fc = os.environ.get("FC", "gfortran")
    flags = ["-std=f2018", "-Werror=implicit-interface", "-ffree-line-length-none",
             os.environ.get("SCIFORT_OPT", "-O0"), "-I", str(build), "-J", str(build)]
    modules = {}
    for source in SRC.rglob("*.f90"):
        for name in re.findall(r"^\s*module\s+(?!procedure|function|subroutine)([a-z]\w*)",
                               source.read_text(), re.I | re.M):
            modules[name.lower()] = source
    visited, active, order = set(), set(), []

    def walk(module):
        if module in visited:
            return
        if module in active:
            raise RuntimeError(f"circular module dependency: {module}")
        active.add(module)
        source = modules[module]
        for dep in re.findall(r"^\s*use\s+(scifort_\w+)", source.read_text(), re.I | re.M):
            walk(dep.lower())
        active.remove(module)
        visited.add(module)
        order.append(source)

    walk("scifort_signal")
    objects = []
    for source in order:
        obj = build / (source.stem + ".o")
        subprocess.run([fc, *flags, "-c", str(source), "-o", str(obj)], cwd=build,
                       check=True, capture_output=True, text=True)
        objects.append(obj)
    exe = build / "resample_callable_window_driver"
    subprocess.run([fc, *flags, str(DRIVER), *map(str, objects), "-o", str(exe)],
                   cwd=build, capture_output=True, text=True, check=True)
    return exe


def make_cases():
    rng = np.random.default_rng(20261026)
    cases = []
    for n in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 19, 20, 23, 29):
        for mode in range(5):
            for repeat in range(5):
                channels = int(rng.integers(1, 5))
                m = [n, max(1, n//2), n+7, 1, int(rng.integers(1, 45))][repeat]
                x = rng.normal(size=(n,channels))
                if mode:
                    x = x+1j*rng.normal(size=x.shape)
                if mode >= 3:
                    x = np.fft.fft(x, axis=0)
                cases.append((n,m,channels,mode,x))
    return cases


def main():
    cases = make_cases()
    lines = []
    ref = []
    for n,m,nc,mode,x in cases:
        lines.append(f"{n} {m} {nc} {mode}\n")
        for v in x.ravel(order="F"):
            lines.append(f"{v.real:.17g} {v.imag:.17g}\n")
        if mode in (0,1,3):
            def window(f):
                return 0.75+0.25*np.sin(11*f)+0.3*f*f
        else:
            def window(f):
                return 1.0+0.2*np.cos(8*f)+1j*(0.4*np.sin(7*f)-0.1*f)
        ref.append(signal.resample(x,m,axis=0,domain="freq" if mode >= 3 else "time",
                                   window=window))
    with tempfile.TemporaryDirectory(prefix="scifort_callable_window_") as tmp:
        exe = os.environ.get("SCIFORT_REFERENCE_EXE") or build_driver(Path(tmp))
        run = subprocess.run([str(exe)],input="".join(lines), capture_output=True,
                             text=True, check=True)
        out = iter(run.stdout.splitlines())
        error_max, total = 0.0, 0
        for index, ((n,m,nc,mode,_),expected) in enumerate(zip(cases,ref)):
            st, nm, nn = map(int, next(out).split())
            if st != 0 or nm != m or nn != nc:
                raise AssertionError(f"case {index} status or shape: {st}, {nm}, {nn}")
            values = np.array([complex(*map(float,next(out).split()))
                               for _ in range(m*nc)]).reshape((m,nc),order="F")
            difference = float(np.max(np.abs(expected-values)))
            error_max = max(error_max,difference)
            np.testing.assert_allclose(values,expected,rtol=8e-12,atol=8e-12,
                err_msg=f"case {index} n={n}, m={m}, channels={nc}, mode={mode}")
            total += values.size
        print(f"PASS: {len(cases)} callable resampling cases, {total} values, "
              f"maximum absolute difference {error_max:.3e}")


if __name__ == "__main__":
    main()
