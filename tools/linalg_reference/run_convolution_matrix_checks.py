#!/usr/bin/env python3
"""Compare SciFort convolution matrices and FFT products with SciPy/NumPy.

The Fortran driver reads cases in one process. A failed comparison exits nonzero.
Use --mode checked or --mode release for different GNU Fortran build options.
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.linalg import convolution_matrix

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    "src/kinds.f90",
    "src/linalg.f90",
    "src/linalg_convolution_matrix.f90",
    "src/signal_windows_extra.f90",
    "src/signal_dpss.f90",
    "src/signal_spectral.f90",
    "src/linalg_matmul_toeplitz.f90",
    "src/linalg_matmul_convolution_matrix.f90",
    "tools/linalg_reference/driver_convolution_matrix.f90",
]


def complex_line(values: np.ndarray) -> str:
    """Format a complex Fortran list-directed input record in column-major order."""
    v = np.asarray(values).ravel(order="F")
    return " ".join(f"({z.real:.17e},{z.imag:.17e})" for z in v) + "\n"


def make_references(seed: int, cases: int) -> bytes:
    """Produce coefficients, right-hand sides, and independent SciPy references."""
    rng = np.random.default_rng(seed)
    lines: list[str] = []
    for idx in range(cases):
        # Force balanced type/mode coverage, including more kernels than samples.
        ctype = 1 + idx % 2
        xtype = 1 + (idx // 2) % 2
        mode = ("full", "same", "valid")[idx % 3]
        m = 1 + int(rng.integers(0, 85))
        n = 1 + int(rng.integers(0, 85))
        if idx % 17 == 0:
            m = 1
        if idx % 19 == 0:
            n = 1
        if idx % 23 == 0:
            m = n
        if idx % 29 == 0:
            m = 2
            n = 7
        if idx % 31 == 0:
            m = 7
            n = 2
        nrhs = 1 + idx % 4
        a = rng.standard_normal(m)
        x = rng.standard_normal((n, nrhs))
        if ctype == 2:
            a = a + 1j * rng.standard_normal(m)
        if xtype == 2:
            x = x + 1j * rng.standard_normal((n, nrhs))
        scale = (1e-80, 1e-20, 1.0, 1e20, 1e80)[idx % 5]
        a *= scale
        expected = convolution_matrix(a, n, mode)
        product = np.column_stack(
            [np.convolve(a, x[:, j], mode=mode) for j in range(nrhs)]
        )
        # The expected matrix and product come from two different SciPy/NumPy APIs.
        assert expected.shape == (len(product), n)
        dense_product = expected @ x
        relative = np.max(np.abs(dense_product - product)) / max(1e-200, np.max(np.abs(product)))
        if relative > 1e-12:
            raise AssertionError(f"Python reference mismatch, case={idx} err={relative}")
        lines.append(f"{m} {n} {nrhs} {mode} {ctype} {xtype}\n")
        lines.extend((complex_line(a), complex_line(x), complex_line(expected), complex_line(product)))
    return "".join(lines).encode()


def run(mode: str, cases: int, seed: int) -> None:
    """Compile focused Fortran modules and process all test cases."""
    flags = (
        ["-std=f2018", "-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
        if mode == "checked"
        else ["-std=f2018", "-O3", "-march=native"]
    )
    with tempfile.TemporaryDirectory(prefix="scifort-convolution-") as temp:
        directory = Path(temp)
        command = ["gfortran", *flags, "-J", str(directory), "-I", str(directory)]
        command.extend(str(ROOT / source) for source in SOURCES)
        command.extend(["-o", str(directory / "driver")])
        subprocess.run(command, cwd=directory, check=True)
        data = make_references(seed, cases)
        result = subprocess.run(
            [str(directory / "driver")], input=data, cwd=directory,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        if result.returncode:
            raise RuntimeError(
                f"Fortran driver failed: {result.stdout.decode()}\n{result.stderr.decode()}"
            )
        print(f"{mode}: {result.stdout.decode().strip()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("checked", "release"), default="checked")
    parser.add_argument("--cases", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=20261109)
    args = parser.parse_args()
    run(args.mode, args.cases, args.seed)
