"""Validate tiled rank-2 oaconvolve_nd against independent SciPy outputs.

Usage: python test/compare_signal_oaconvolve_tiled_scipy.py [gfortran-path]
Compiles the required signal modules in dependency order.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import oaconvolve

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [
    ROOT / "src" / name
    for name in ("kinds.f90", "signal_dpss.f90", "signal_windows_extra.f90", "signal_spectral.f90", "signal_fftconvolve.f90", "signal_convolve.f90")
]
DRIVER = ROOT / "validation" / "reference_drivers" / "signal_oaconvolve_tiled_reference_driver.f90"


def run() -> None:
    compiler = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FC", "gfortran")
    # Each shape covers a distinct block partition, broadcast or axis path.
    shape_cases = [
        ((129, 65), (3, 5), None),
        ((3, 5), (129, 65), None),
        ((129, 3), (5, 131), None),
        ((65, 129), (3, 5), None),
        ((2, 129), (1, 5), (1,)),
        ((1, 5), (2, 129), (1,)),
        ((129, 2), (5, 1), (0,)),
        ((5, 1), (129, 2), (0,)),
        ((129, 1), (3, 2), (0,)),
        ((3, 2), (129, 1), (0,)),
    ]
    rng = np.random.default_rng(20261019)
    with tempfile.TemporaryDirectory(prefix="scifort-oa-") as dirname:
        tmp = Path(dirname)
        exe = tmp / ("reference.exe" if os.name == "nt" else "reference")
        subprocess.run(
            [compiler, "-std=f2018", "-O2", "-fcheck=all", "-J", str(tmp),
             "-I", str(tmp), *map(str, SOURCES), str(DRIVER), "-o", str(exe)],
            cwd=tmp, check=True,
        )
        input_path = tmp / "input.bin"
        output_path = tmp / "output.bin"
        cases = []
        with input_path.open("wb") as writer:
            for shape_a, shape_b, axes in shape_cases:
                for mode in ("full", "same", "valid"):
                    conv_axes = range(2) if axes is None else axes
                    if mode == "valid" and not (
                        all(shape_a[j] >= shape_b[j] for j in conv_axes)
                        or all(shape_a[j] <= shape_b[j] for j in conv_axes)
                    ):
                        continue
                    for typ in range(1, 5):
                        for _ in range(2):
                            a = rng.normal(size=shape_a)
                            b = rng.normal(size=shape_b)
                            if typ in (2, 4):
                                a = a + 1j * rng.normal(size=shape_a)
                            if typ in (2, 3):
                                b = b + 1j * rng.normal(size=shape_b)
                            a = np.asfortranarray(a)
                            b = np.asfortranarray(b)
                            expected = oaconvolve(a, b, mode=mode, axes=axes)
                            axis_code = 0 if axes is None else axes[0] + 1
                            header = np.array(
                                [*shape_a, *shape_b, *expected.shape, axis_code,
                                 ("full", "same", "valid").index(mode) + 1, typ],
                                dtype=np.int32,
                            )
                            writer.write(header.tobytes())
                            writer.write(a.tobytes(order="F"))
                            writer.write(b.tobytes(order="F"))
                            cases.append((expected, typ, shape_a, shape_b, mode, axes))
        subprocess.run([str(exe), str(input_path), str(output_path)], check=True)
        raw = output_path.read_bytes()
        location = 0
        maximum = 0.0
        for i, (expected, typ, shape_a, shape_b, mode, axes) in enumerate(cases):
            status = np.frombuffer(raw, dtype=np.int32, count=1, offset=location)[0]
            location += 4
            if status != 0:
                raise AssertionError(f"Case {i} {shape_a} {shape_b} {mode}: status {status}")
            dtype = np.dtype(np.float64 if typ == 1 else np.complex128)
            number = expected.size
            got = np.frombuffer(raw, dtype=dtype, count=number, offset=location).reshape(
                expected.shape, order="F"
            )
            location += number * dtype.itemsize
            error = float(np.max(np.abs(got - expected)))
            maximum = max(maximum, error)
            if not np.allclose(got, expected, atol=2e-11, rtol=2e-11):
                raise AssertionError(
                    f"Case {i} {shape_a} {shape_b} {mode} {axes} "
                    f"type={typ}: max error={error:.3e}"
                )
        if location != len(raw):
            raise AssertionError("Trailing output from the Fortran driver")
        print(f"SciPy oaconvolve comparisons: {len(cases)} passed; max abs error {maximum:.3e}")


if __name__ == "__main__":
    run()
