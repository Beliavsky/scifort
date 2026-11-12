#!/usr/bin/env python3
"""Independent SciPy comparisons for the SciFort circulant solver.

This compiles a dedicated Fortran I/O driver and tests the actual Fortran
implementation against scipy.linalg.solve_circulant, not a Python port.
Run from the project root. No FPM dependency is required.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import LinAlgError, circulant, solve_circulant

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src"
DRIVER = ROOT / "tools" / "linalg_reference" / "driver_circulant.f90"


def compile_driver(build: Path, mode: str) -> Path:
    build.mkdir(parents=True, exist_ok=True)
    flags = ["-std=f2018"]
    if mode == "checked":
        flags += ["-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
    else:
        flags += ["-O3", "-march=native"]
    objects = []
    modules = [
        "kinds.f90", "linalg.f90", "signal_windows_extra.f90",
        "signal_dpss.f90", "signal_spectral.f90", "linalg_circulant.f90",
    ]
    for filename in modules:
        obj = build / (filename.removesuffix(".f90") + ".o")
        subprocess.run(["gfortran", *flags, "-J", str(build), "-I", str(build),
                        "-c", str(SOURCE / filename), "-o", str(obj)],
                       check=True, capture_output=True, text=True)
        objects.append(str(obj))
    exe = build / "driver_circulant"
    subprocess.run(["gfortran", *flags, "-I", str(build), *objects,
                    str(DRIVER), "-o", str(exe)],
                   check=True, capture_output=True, text=True)
    return exe


def make_cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    orders = [1, 2, 3, 4, 5, 7, 8, 11, 12, 16, 17, 19, 20, 31,
              32, 37, 48, 64, 67, 97, 128]
    cases = []
    for index in range(count):
        n = orders[index % len(orders)]
        nrhs = [1, 2, 4][index % 3]
        kind = index % 4
        scenario = (index // 4) % 8
        scale = 10.0 ** [-80, -12, 0, 18, 80][(index // 5) % 5]
        c = (0.12 * rng.standard_normal(n) + 0.12j * rng.standard_normal(n))
        c[0] = (2.0 + 0.05 * n) + 0.05j
        b = rng.standard_normal((n, nrhs)) + 1j * rng.standard_normal((n, nrhs))
        singular = "lstsq" if scenario in (2, 3, 5, 7) else "raise"
        tol = None
        if scenario in (2, 3) and n > 1:
            c = np.ones(n, dtype=complex)  # exactly singular at all but DC
        if scenario in (4, 5) and n > 1:
            c = np.zeros(n, dtype=complex)
            c[0] = 1.0 + 0.0j
            c[1] = 1.0 - 1e-15 + 0.0j
            if n != 2:
                c[2] = 0.8
            tol = 1e-12 * max(abs(np.fft.fft(c)))
        if scenario in (6, 7) and n > 1:
            c = np.zeros(n, dtype=complex)  # zero circulant
        if scenario == 1 and n > 2:
            c[1] = -0.01
        if kind in (0, 2):
            c = c.real
        if kind in (0, 3):
            b = b.real
        c = c * scale
        b = b * scale
        if tol is not None:
            tol *= scale
        cases.append((c, b, kind, singular, tol, n, nrhs))
    return cases


def format_line(arr):
    return " ".join(f"{v:.17e}" for v in arr) + "\n"


def run(cases, exe: Path):
    payload = []
    for c, b, kind, policy, tol, n, nrhs in cases:
        payload.append(f"{n} {nrhs} {kind} {int(policy == 'lstsq')} "
                       f"{tol if tol is not None else -1.0:.17e}\n")
        payload.append(format_line(np.real(c)))
        payload.append(format_line(np.imag(c)))
        for j in range(nrhs):
            payload.append(format_line(np.real(b[:, j])))
            payload.append(format_line(np.imag(b[:, j])))
    payload.append("-1 0 0 0 -1\n")
    proc = subprocess.run([str(exe)], input="".join(payload),
                          text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode:
        raise RuntimeError(f"Fortran driver returned {proc.returncode}: {proc.stderr[-3000:]}")
    lines = proc.stdout.splitlines()
    pos, maximum_error, maximum_residual, passed, singular_count = 0, 0.0, 0.0, 0, 0
    for index, (c, b, kind, policy, tol, n, nrhs) in enumerate(cases):
        label = f"#{index}: n={n} rhs={nrhs} kind={kind} policy={policy} tol={tol}"
        try:
            # SciPy 1.17 has a mask-broadcast issue for singular c with
            # b shaped (n,nrhs). Use its supported trailing-vector axis.
            expected = solve_circulant(c, b.T, singular=policy, tol=tol,
                                       baxis=-1, outaxis=-1).T
            expected_singular = False
        except LinAlgError:
            expected_singular = True
            expected = None
        if pos >= len(lines):
            raise AssertionError(f"Fortran output missing at {label}")
        status = int(lines[pos])
        pos += 1
        if expected_singular:
            if status != 5:
                raise AssertionError(f"{label}: expected singular (5), got {status}")
            singular_count += 1
        else:
            if status != 0:
                raise AssertionError(f"{label}: expected success (0), got {status}")
            raw = []
            for _ in range(n * nrhs):
                if pos >= len(lines):
                    raise AssertionError(f"{label}: truncated Fortran output")
                r, im = (float(s) for s in lines[pos].split())
                raw.append(r + 1j * im)
                pos += 1
            result = np.array(raw).reshape(nrhs, n).T
            # Test with non-dimensional forward error and normalized residual.
            err = np.linalg.norm(result-expected) / max(np.linalg.norm(expected), 1e-20)
            maximum_error = max(maximum_error, err)
            cc = circulant(c)
            residual = np.linalg.norm(cc @ result-b) / max(
                np.linalg.norm(cc)*np.linalg.norm(result) + np.linalg.norm(b), 1e-250)
            if policy == "raise":
                maximum_residual = max(maximum_residual, residual)
            if err > 2e-10:
                raise AssertionError(f"{label}: relative error {err:.3e}")
            if policy == 'raise' and residual > 2e-10:
                raise AssertionError(f"{label}: normalized residual {residual:.3e}")
            if kind == 0 and np.max(abs(result.imag)) > 1e-14:
                raise AssertionError(f"{label}: real-real produced complex solution")
        passed += 1
    if pos != len(lines):
        raise AssertionError(f"Fortran produced {len(lines)-pos} trailing lines")
    return passed, singular_count, maximum_error, maximum_residual


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261109)
    parser.add_argument("--mode", choices=("checked", "release"), default="checked")
    parser.add_argument("--build", type=Path)
    args = parser.parse_args()
    if args.cases < 0:
        parser.error("cases must be nonnegative")
    build = args.build or ROOT / "_circulant_reference_build" / args.mode
    exe = compile_driver(build, args.mode)
    cases = make_cases(args.cases, args.seed)
    passed, singular_count, maxerr, maxres = run(cases, exe)
    print(f"SciPy {scipy.__version__}, NumPy {np.__version__}; {args.mode} build")
    print(f"PASS {passed}/{args.cases}; singular detections {singular_count}; "
          f"max relative solution error {maxerr:.3e}; "
          f"max normalized nonsingular residual {maxres:.3e}")


if __name__ == "__main__":
    main()
