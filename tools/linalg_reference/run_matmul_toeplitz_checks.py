#!/usr/bin/env python3
"""Compile and compare SciFort matmul_toeplitz with SciPy and dense NumPy.

Independent randomized tests of the actual Fortran routines, including all
four real/complex input combinations, square and rectangular matrices,
explicit/default rows, and several RHS columns. Run from anywhere.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import matmul_toeplitz, toeplitz

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src"
DRIVER = ROOT / "tools" / "linalg_reference" / "driver_matmul_toeplitz.f90"
MODULES = (
    "kinds.f90", "linalg.f90", "signal_windows_extra.f90",
    "signal_dpss.f90", "signal_spectral.f90", "linalg_matmul_toeplitz.f90",
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
    executable = build / "driver_matmul_toeplitz"
    subprocess.run(["gfortran", *flags, "-I", str(build), *objects,
                    str(DRIVER), "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    return executable


def make_cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = [1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 16, 17, 19, 20, 23, 29,
             31, 32, 37, 40, 47, 48, 53, 59, 64, 67, 73, 96, 97, 101, 128]
    cases = []
    for index in range(count):
        n = sizes[(index + index // 31) % len(sizes)]
        explicit = index % 4 != 1
        m = sizes[(index * 13 + index // 31 + 7) % len(sizes)] if explicit else n
        nrhs = (1, 2, 4)[index % 3]
        kind = (index // 3) % 4
        scale = 10.0 ** (-80, -20, 0, 20, 80)[(index // 11) % 5]
        c = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        r = rng.standard_normal(m) + 1j * rng.standard_normal(m)
        x = rng.standard_normal((m, nrhs)) + 1j * rng.standard_normal((m, nrhs))
        # Include exactly zero coefficients and large first row disagreements.
        if index % 17 == 0:
            c[:] = 0
        if index % 13 == 0 and explicit:
            r[0] = 1234567 - 123456j   # must be ignored, even if not c[0]
        if kind in (0, 1):
            c = c.real
            r = r.real
        if kind in (0, 2):
            x = x.real
        c *= scale
        r *= scale
        x *= 10.0 ** (-20, 0, 20)[(index // 19) % 3]
        cases.append((c, r if explicit else None, x, kind, n, m, nrhs))
    return cases


def run_cases(executable: Path, cases):
    payload = []
    for c, r, x, kind, n, m, nrhs in cases:
        payload.append(f"{n} {m} {nrhs} {kind} {int(r is not None)}\n")
        for entry in c:
            payload.append(f"{np.real(entry):.17e} {np.imag(entry):.17e}\n")
        if r is not None:
            for entry in r:
                payload.append(f"{np.real(entry):.17e} {np.imag(entry):.17e}\n")
        for j in range(nrhs):
            for entry in x[:, j]:
                payload.append(f"{np.real(entry):.17e} {np.imag(entry):.17e}\n")
    payload.append("-1 0 0 0 0\n")
    proc = subprocess.run([str(executable)], input="".join(payload), text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode:
        raise RuntimeError(f"Fortran driver failed: {proc.stderr[-3000:]}")
    output = proc.stdout.splitlines()
    pos = passed = 0
    maximum_error = maximum_dense_error = 0.0
    for index, (c, r, x, kind, n, m, nrhs) in enumerate(cases):
        label = f"case {index}: n={n} m={m} rhs={nrhs} kind={kind} row={r is not None}"
        if pos >= len(output):
            raise AssertionError(f"Missing Fortran output: {label}")
        status = int(output[pos])
        pos += 1
        if status != 0:
            raise AssertionError(f"{label}: Fortran status {status}")
        numbers = []
        for _ in range(n * nrhs):
            real, imag = (float(value) for value in output[pos].split())
            pos += 1
            numbers.append(real + 1j * imag)
        actual = np.array(numbers).reshape(nrhs, n).T
        arg = (c, r) if r is not None else c
        expected = matmul_toeplitz(arg, x)
        independent = toeplitz(c, r) @ x
        # Relative norm measured against SciPy (not the coefficient scale).
        denom = max(np.linalg.norm(expected), 1e-250)
        error = np.linalg.norm(actual-expected) / denom
        dense_error = np.linalg.norm(actual-independent) / max(np.linalg.norm(independent), 1e-250)
        if not np.isfinite(error) or not np.isfinite(dense_error):
            raise AssertionError(f"Nonfinite error: {label}")
        if error > 3e-11 or dense_error > 3e-11:
            raise AssertionError(f"{label}: SciPy {error:.3e}; dense {dense_error:.3e}")
        if kind == 0 and np.any(actual.imag != 0.0):
            raise AssertionError(f"{label}: real input produced complex output")
        maximum_error = max(maximum_error, error)
        maximum_dense_error = max(maximum_dense_error, dense_error)
        passed += 1
    if pos != len(output):
        raise AssertionError(f"Fortran emitted {len(output)-pos} unexpected lines")
    return passed, maximum_error, maximum_dense_error


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("checked", "release"), default="checked")
    parser.add_argument("--cases", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261109)
    parser.add_argument("--build", type=Path)
    args = parser.parse_args()
    if args.cases < 0:
        parser.error("cases must be nonnegative")
    build = args.build or ROOT / "_matmul_toeplitz_reference_build" / args.mode
    executable = compile_driver(build, args.mode)
    cases = make_cases(args.cases, args.seed)
    passed, maximum_error, maximum_dense_error = run_cases(executable, cases)
    print(f"SciPy {scipy.__version__}, NumPy {np.__version__}; {args.mode} build")
    print(f"PASS {passed}/{args.cases}; max relative SciPy difference "
          f"{maximum_error:.3e}; max relative dense difference {maximum_dense_error:.3e}")


if __name__ == "__main__":
    main()
