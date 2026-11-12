#!/usr/bin/env python3
"""Validate nonsymmetric Schur-subspace locking against dense SciPy eigenpairs.

Fortran returns actual complex right eigenvectors, so reference checks compute
independent true residuals ||Av-lambda v|| as well as matched eigenvalues.
Eigenvector phases and ordering within ties are not compared.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]


def priority(values: np.ndarray, which: str) -> np.ndarray:
    return {"LM": lambda x: abs(x), "LR": lambda x: x.real,
            "SR": lambda x: -x.real, "LI": lambda x: x.imag,
            "SI": lambda x: -x.imag}[which](values)


def matrix(rng: np.random.Generator, n: int, mode: str, shift: complex) -> np.ndarray:
    a = np.diag(rng.uniform(-0.6, 0.6, n))
    for j in range(n - 1):
        a[j, j+1] = rng.uniform(0.015, 0.085)
        a[j+1, j] = -rng.uniform(0.008, 0.045)
    if shift:
        if shift.imag:
            a[:4, :4] = 0.0
            a[0, 1], a[1, 0] = -5.0, 5.0
            a[2, 3], a[3, 2] = -3.8, 3.8
        else:
            a[0, 0], a[1, 1] = 1.0, 1.4
    elif mode in ("LR", "LM"):
        a[0, 0] = rng.uniform(7.5, 9.0)
        a[1, 1] = rng.uniform(4.5, 6.0)
    elif mode == "SR":
        a[0, 0] = -rng.uniform(7.5, 9.0)
        a[1, 1] = -rng.uniform(4.5, 6.0)
    else:
        a[:4, :4] = 0.0
        a[0, 1], a[1, 0] = -rng.uniform(7, 9), rng.uniform(7, 9)
        a[2, 3], a[3, 2] = -rng.uniform(4, 5.5), rng.uniform(4, 5.5)
    # Non-trivial upper-right coupling is essential: a projected-right-vector
    # scheme can incorrectly discard components in a locked invariant block.
    a[0, 5] = rng.uniform(-4, 4)
    a[1, 7] = rng.uniform(-3, 3)
    return a


def build_driver(build: Path, mode: str) -> Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    cc = ["gfortran", "-std=f2018", *flags, "-J", str(build), "-I", str(build)]
    objects = []
    for name, rel in (("kinds", "src/kinds.f90"),
                      ("sparse", "src/sparse.f90"),
                      ("lsmr", "src/optimize/lsmr.f90"),
                      ("sparse_linalg", "src/sparse_linalg.f90")):
        obj = build / (name + ".o")
        subprocess.run([*cc, "-c", str(ROOT / rel), "-o", str(obj)],
                       cwd=build, check=True)
        objects.append(obj)
    executable = build / "schur_reference"
    subprocess.run([*cc, str(ROOT / "tools/sparse_linalg_reference" /
                                  "sparse_schur_reference_driver.f90"),
                    *(str(obj) for obj in objects), "-o", str(executable)],
                   cwd=build, check=True)
    return executable


def run_one(executable: Path, a: np.ndarray, mode: str, shift: complex,
            ncv: int, k: int = 2, cap: int = 90) -> tuple[int, int]:
    n = len(a)
    args = f"{n} {k} {ncv} {cap} {mode} 1.0e-8 {shift.real} {shift.imag}\n"
    args += " ".join(repr(float(v)) for v in a.flatten(order="F")) + "\n"
    output = subprocess.run([str(executable)], input=args, capture_output=True,
                            text=True, check=True, timeout=35).stdout.splitlines()
    status, cycles, converged, frozen, products = map(int, output[0].split())
    if status != 0:
        raise AssertionError(f"nonconvergence status={status}, mode={mode}, "
                             f"shift={shift}, cycles={cycles}")
    if converged != k or products < ncv or cycles < 1:
        raise AssertionError(f"invalid diagnostics: {output[0]}")
    if len(output) != 1 + k * (n + 1):
        raise AssertionError("unexpected driver output length")
    values = np.empty(k, dtype=complex)
    errors = np.empty(k)
    vectors = np.empty((n, k), dtype=complex)
    for j in range(k):
        base = 1 + j * (n + 1)
        r, im, err = map(float, output[base].split())
        values[j], errors[j] = complex(r, im), err
        vectors[:, j] = [complex(*map(float, line.split()))
                          for line in output[base+1:base+1+n]]
    reference = eigvals(a)
    if shift:
        ranking = abs(1.0 / (reference - shift))
    else:
        ranking = priority(reference, mode)
    expected = reference[np.argsort(-ranking)[:k]]
    rows, cols = linear_sum_assignment(abs(values[:, None] - expected[None, :]))
    delta = max(abs(values[rows] - expected[cols]))
    if delta > 2e-4:
        raise AssertionError(f"incorrect eigenvalues {mode=} {shift=} {delta=} "
                             f"returned={values}, reference={expected}")
    true_errors = np.linalg.norm(a @ vectors - vectors * values[None, :], axis=0)
    if not np.allclose(true_errors, errors, rtol=1e-6, atol=1e-10):
        raise AssertionError(f"incorrect reported true residuals: {true_errors=} {errors=}")
    if max(true_errors) > 1e-6 * max(1, np.linalg.norm(a, ord=np.inf)):
        raise AssertionError(f"poor eigenpair residuals: {true_errors=}")
    if not np.allclose(np.linalg.norm(vectors, axis=0), 1.0, atol=1e-9):
        raise AssertionError("non-unit eigenvectors")
    return cycles, frozen


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    p.add_argument("--cases", type=int, default=70)
    p.add_argument("--seed", type=int, default=20261104)
    args = p.parse_args()
    if args.cases < 1:
        p.error("--cases must be positive")
    for build_mode in (("checked", "optimized") if args.mode == "both" else (args.mode,)):
        rng = np.random.default_rng(args.seed)
        with tempfile.TemporaryDirectory(prefix="scifort-schur-") as tmp:
            exe = build_driver(Path(tmp), build_mode)
            cycles = []
            frozen = []
            cases = [("LR", 0j), ("LM", 0j), ("SR", 0j),
                     ("LI", 0j), ("SI", 0j),
                     ("LM", 1.03+0j), ("LM", 0.0+4.9j)]
            for i in range(args.cases):
                mode, shift = cases[i % len(cases)]
                n = (24, 40, 60)[(i // len(cases)) % 3]
                a = matrix(rng, n, mode, shift)
                steps, locked = run_one(exe, a, mode, shift, 7 if n < 60 else 9)
                cycles.append(steps)
                frozen.append(locked)
            multi_locks = []
            for _ in range(12):
                n = 64
                a = np.diag(0.6 * np.cos(0.24 * np.arange(1, n+1)))
                for j in range(n-1):
                    a[j, j+1] = 0.07 * rng.uniform(0.5, 1.0)
                    a[j+1, j] = -0.025 * rng.uniform(0.5, 1.0)
                a[0, 0], a[1, 1], a[2, 2] = 9.0, 7.0, 5.0
                a[0, 4], a[1, 5] = rng.uniform(0, 2), rng.uniform(0, 2)
                steps, locked = run_one(exe, a, 'LR', 0j, 5, k=3, cap=75)
                if locked < 2:
                    raise AssertionError(f"expected two-dimensional Schur lock, got {locked}")
                multi_locks.append((steps, locked))
            print(f"{build_mode}: {args.cases} SciPy spectral comparisons PASSED; "
                  f"{sum(l > 0 for l in frozen)} locked invariant subspaces; "
                  f"{sum(s > 1 for s in cycles)} restarted cases; "
                  f"maximum {max(cycles)} cycles; "
                  f"12 additional three-eigenpair cases with 2+ Schur locks PASSED", flush=True)


if __name__ == "__main__":
    main()
