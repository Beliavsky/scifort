#!/usr/bin/env python3
"""Compile SciFort Hankel routines and compare their outputs with SciPy.

Tests the actual Fortran implementations against scipy.linalg.hankel @ x,
with a separate hand-built anti-diagonal dense NumPy reference.  Some
cases compare the Fortran Hankel constructor element by element as well.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import hankel

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src"
DRIVER = ROOT / "tools" / "linalg_reference" / "driver_hankel.f90"
MODULES = (
    "kinds.f90", "linalg.f90", "signal_windows_extra.f90",
    "signal_dpss.f90", "signal_spectral.f90",
    "linalg_matmul_toeplitz.f90", "linalg_hankel.f90",
    "linalg_matmul_hankel.f90",
)


def compile_driver(build: Path, mode: str) -> Path:
    build.mkdir(parents=True, exist_ok=True)
    flags = ["-std=f2018"]
    if mode == "checked":
        flags += ["-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
    else:
        flags += ["-O3", "-march=native"]
    objects = []
    for filename in MODULES:
        obj = build / (Path(filename).stem + ".o")
        subprocess.run(["gfortran", *flags, "-J", str(build), "-I", str(build),
                        "-c", str(SOURCE / filename), "-o", str(obj)],
                       check=True, capture_output=True, text=True)
        objects.append(str(obj))
    executable = build / "driver_hankel"
    subprocess.run(["gfortran", *flags, "-I", str(build), *objects,
                    str(DRIVER), "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    return executable


def make_cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = [0, 1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 16, 17, 19, 20, 23, 29,
             31, 32, 37, 40, 47, 48, 53, 59, 64, 67, 73, 96, 97, 101, 128]
    cases = []
    for index in range(count):
        n = sizes[(index + index // 32) % len(sizes)]
        explicit = index % 4 != 1
        m = sizes[(index * 13 + index // 32 + 5) % len(sizes)] if explicit else n
        nrhs = (1, 2, 4)[index % 3]
        kind = (index // 3) % 4
        scale = 10.0 ** (-80, -20, 0, 20, 80)[(index // 11) % 5]
        c = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        r = rng.standard_normal(m) + 1j * rng.standard_normal(m)
        x = rng.standard_normal((m, nrhs)) + 1j * rng.standard_normal((m, nrhs))
        if index % 17 == 0:
            c[:] = 0
        if index % 13 == 0 and explicit and m:
            r[0] = 1234567 - 123456j  # ignored first element of last row
        if kind in (0, 1):
            c = c.real
            r = r.real
        if kind in (0, 2):
            x = x.real
        c *= scale
        r *= scale
        x *= 10.0 ** (-20, 0, 20)[(index // 19) % 3]
        emit_dense = (index % 7 == 0 and n * m <= 4000)
        cases.append((c, r if explicit else None, x, kind, n, m, nrhs, emit_dense))
    return cases


def manual_hankel(c: np.ndarray, r: np.ndarray | None) -> np.ndarray:
    n = c.size
    m = n if r is None else r.size
    if n == 0 or m == 0:
        return np.zeros((n, m), dtype=np.result_type(c, r if r is not None else c))
    stream = np.concatenate([c, np.zeros(m - 1, dtype=c.dtype) if r is None else r[1:]])
    positions = np.arange(n)[:, None] + np.arange(m)[None, :]
    return stream[positions]


def run_cases(executable: Path, cases):
    payload = []
    for c, r, x, kind, n, m, nrhs, emit_dense in cases:
        payload.append(f"{n} {m} {nrhs} {kind} {int(r is not None)} {int(emit_dense)}\n")
        for entry in c:
            payload.append(f"{np.real(entry):.17e} {np.imag(entry):.17e}\n")
        if r is not None:
            for entry in r:
                payload.append(f"{np.real(entry):.17e} {np.imag(entry):.17e}\n")
        for j in range(nrhs):
            for entry in x[:, j]:
                payload.append(f"{np.real(entry):.17e} {np.imag(entry):.17e}\n")
    payload.append("-1 0 0 0 0 0\n")
    proc = subprocess.run([str(executable)], input="".join(payload), text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode:
        raise RuntimeError(f"Fortran driver failed: {proc.stderr[-3000:]}")
    output = proc.stdout.splitlines()
    pos = passed = 0
    maximum_error = maximum_dense_error = maximum_hankel_error = 0.0
    constructor_cases = 0
    for index, (c, r, x, kind, n, m, nrhs, emit_dense) in enumerate(cases):
        label = (f"case {index}: n={n} m={m} rhs={nrhs} kind={kind} "
                 f"row={r is not None}")
        if pos >= len(output):
            raise AssertionError(f"Missing Fortran output: {label}")
        status = int(output[pos])
        pos += 1
        if status != 0:
            raise AssertionError(f"{label}: Fortran status {status}")
        numbers = []
        for _ in range(n * nrhs):
            re, im = (float(value) for value in output[pos].split())
            pos += 1
            numbers.append(re + 1j * im)
        actual = np.array(numbers).reshape(nrhs, n).T
        expected_matrix = hankel(c, r)
        independent_matrix = manual_hankel(c, r)
        if not np.allclose(expected_matrix, independent_matrix, rtol=0, atol=0):
            raise AssertionError(f"SciPy vs independently constructed Hankel: {label}")
        scipy_result = expected_matrix @ x
        dense_result = independent_matrix @ x
        denom = max(np.linalg.norm(scipy_result), 1e-250)
        error = np.linalg.norm(actual-scipy_result) / denom
        dense_error = np.linalg.norm(actual-dense_result) / max(np.linalg.norm(dense_result), 1e-250)
        if not np.isfinite(error) or not np.isfinite(dense_error):
            raise AssertionError(f"Nonfinite product error: {label}")
        if error > 3e-11 or dense_error > 3e-11:
            raise AssertionError(f"{label}: SciPy {error:.3e}, dense {dense_error:.3e}")
        if kind == 0 and np.any(actual.imag != 0.0):
            raise AssertionError(f"{label}: real inputs produced complex output")
        maximum_error = max(maximum_error, error)
        maximum_dense_error = max(maximum_dense_error, dense_error)
        if emit_dense:
            numbers = []
            for _ in range(n * m):
                re, im = (float(value) for value in output[pos].split())
                pos += 1
                numbers.append(re + 1j * im)
            actual_matrix = np.array(numbers).reshape(m, n).T
            matrix_error = np.linalg.norm(actual_matrix-expected_matrix) / max(
                np.linalg.norm(expected_matrix), 1e-250)
            if not np.isfinite(matrix_error) or matrix_error > 2e-14:
                raise AssertionError(f"{label}: Hankel constructor {matrix_error:.3e}")
            maximum_hankel_error = max(maximum_hankel_error, matrix_error)
            constructor_cases += 1
        passed += 1
    if pos != len(output):
        raise AssertionError(f"Fortran emitted {len(output)-pos} unexpected lines")
    return passed, constructor_cases, maximum_error, maximum_dense_error, maximum_hankel_error


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("checked", "release"), default="checked")
    parser.add_argument("--cases", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=20261109)
    parser.add_argument("--build", type=Path)
    args = parser.parse_args()
    if args.cases < 0:
        parser.error("cases must be nonnegative")
    build = args.build or ROOT / "_hankel_reference_build" / args.mode
    exe = compile_driver(build, args.mode)
    cases = make_cases(args.cases, args.seed)
    passed, constructors, max_scipy, max_dense, max_constructor = run_cases(exe, cases)
    print(f"SciPy {scipy.__version__}, NumPy {np.__version__}; {args.mode} build")
    print(f"PASS {passed}/{args.cases} products; {constructors} constructors; "
          f"max relative SciPy difference {max_scipy:.3e}; "
          f"max independent dense difference {max_dense:.3e}; "
          f"max constructor difference {max_constructor:.3e}")


if __name__ == "__main__":
    main()
