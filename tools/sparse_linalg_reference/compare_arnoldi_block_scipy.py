#!/usr/bin/env python3
"""Test block restarted Arnoldi against SciPy dense eigensystems and residuals.

In exact repeated eigenspaces an ARPACK single-start call can return just one
member, so eigenvalue *multiplicity* comes from scipy.linalg.eigvals. This
also independently validates right eigenpairs, unit norms and independence.
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


def build_driver(build: Path, mode: str) -> Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    cc = ["gfortran", "-std=f2018", *flags, "-J", str(build), "-I", str(build)]
    objects = []
    for name, rel in (("kinds", "src/kinds.f90"), ("sparse", "src/sparse.f90"),
                      ("lsmr", "src/optimize/lsmr.f90"),
                      ("sparse_linalg", "src/sparse_linalg.f90")):
        obj = build / f"{name}.o"
        subprocess.run([*cc, "-c", str(ROOT / rel), "-o", str(obj)], cwd=build, check=True)
        objects.append(obj)
    exe = build / "arnoldi_block_reference"
    subprocess.run([*cc, str(ROOT / "tools/sparse_linalg_reference" /
                                   "sparse_arnoldi_block_reference_driver.f90"),
                    *(str(obj) for obj in objects), "-o", str(exe)],
                   cwd=build, check=True)
    return exe


def make_matrix(rng: np.random.Generator, n: int, scenario: int):
    a = np.diag(rng.uniform(-0.25, 0.25, n))
    a[np.arange(n-1), np.arange(1, n)] = rng.uniform(0.015, 0.065, n-1)
    k, b, which, shift = 2, 2, "LR", 0j
    if scenario in (0, 7):
        a[0, 0], a[1, 1], a[2, 2] = 9.0, 9.0, 8.9
        a[0, 1] = 0.0  # Semisimple repeated eigenvalue.
        a[0, 5], a[1, 6] = rng.uniform(0.1, 0.6), -rng.uniform(0.1, 0.6)
        if scenario == 7:
            which = "LM"
    elif scenario == 1:
        a[0, 0], a[1, 1], a[2, 2] = 9.0, 8.999, 8.998
        a[0, 1] = a[1, 2] = 0.0
        k, b = 3, 3
    elif scenario == 2:
        a[0, 0], a[1, 1], a[2, 2] = -9.0, -9.0, -8.9
        a[0, 1] = 0.0
        which = "SR"
    elif scenario in (3, 4, 6):
        a[:4, :4] = 0.0
        a[0, 1], a[1, 0] = -8.0, 8.0
        a[2, 3], a[3, 2] = -8.0, 8.0
        a[0, 6] = rng.uniform(0.1, 0.6)
        a[2, 7] = rng.uniform(0.1, 0.6)
        which = "LI" if scenario != 4 else "SI"
        if scenario == 6:
            which, shift = "LM", 0.0 + 8.05j
    elif scenario == 5:
        a[0, 0], a[1, 1], a[2, 2] = 1.0, 1.0, 1.25
        a[0, 1] = 0.0
        which, shift = "LM", 1.025 + 0.0j
    elif scenario == 8:
        a[0, 0], a[1, 1], a[2, 2] = 9.0, 7.0, 5.0
        a[0, 4], a[1, 5] = rng.uniform(0.5, 2), rng.uniform(-1, -0.2)
        k, b = 3, 2
    else:
        raise AssertionError(f"unexpected scenario {scenario}")
    return a, k, b, which, shift


def priority(vals: np.ndarray, target: str, sigma: complex) -> np.ndarray:
    if sigma:
        return 1.0 / abs(vals - sigma)
    return {"LR": lambda v: v.real, "SR": lambda v: -v.real,
            "LI": lambda v: v.imag, "SI": lambda v: -v.imag,
            "LM": abs}[target](vals)


def test_one(exe: Path, a: np.ndarray, k: int, b: int, target: str,
             sigma: complex, locking: bool) -> tuple[int, int]:
    n = len(a)
    p = max(10, k + b + 5)
    cap = 90
    args = f"{n} {k} {p} {cap} {b} {target} 1e-8 {sigma.real} {sigma.imag} {int(locking)}\n"
    args += " ".join(repr(float(x)) for x in a.flatten(order="F")) + "\n"
    result = subprocess.run([str(exe)], input=args, capture_output=True,
                            text=True, timeout=30)
    if result.returncode != 0:
        raise AssertionError(f"native crash: {result.stderr}")
    output = result.stdout.splitlines()
    status, cycles, wanted, locked, products = map(int, output[0].split())
    if status != 0 or wanted != k:
        raise AssertionError(f"nonconvergence: {status=} {cycles=} {wanted=} "
                             f"{locked=} {k=} {b=} {target=} {sigma=}")
    if len(output) != 1 + k * (n + 1):
        raise AssertionError(f"unexpected output size {len(output)}")
    values = np.zeros(k, dtype=complex)
    vectors = np.zeros((n, k), dtype=complex)
    errors = np.zeros(k)
    for j in range(k):
        first = 1 + j * (n + 1)
        re, im, err = map(float, output[first].split())
        values[j], errors[j] = re + 1j * im, err
        vectors[:, j] = [complex(*map(float, entry.split()))
                         for entry in output[first+1:first+1+n]]
    true_errors = np.linalg.norm(a @ vectors - vectors*values[None, :], axis=0)
    if max(true_errors) > 2e-6 or not np.allclose(true_errors, errors, rtol=1e-5, atol=1e-10):
        raise AssertionError(f"incorrect eigenpair residuals: {true_errors=} {errors=}")
    if not np.allclose(np.linalg.norm(vectors, axis=0), 1, atol=2e-9):
        raise AssertionError("non-unit eigenvectors")
    all_eigenvalues = eigvals(a)
    reference = all_eigenvalues[np.argsort(-priority(all_eigenvalues, target, sigma))[:k]]
    rows, cols = linear_sum_assignment(abs(values[:, None]-reference[None, :]))
    delta = max(abs(values[rows]-reference[cols]))
    if delta > 2e-4:
        raise AssertionError(f"wrong eigenvalues: {delta=} {values=} {reference=} {sigma=}")
    # Semisimple repeated eigenspaces require independent eigenvectors, not
    # just two copies of one accurate eigenvector.
    if k == 2 and abs(reference[0] - reference[1]) < 1e-10:
        gram_determinant = 1.0 - abs(np.vdot(vectors[:, 0], vectors[:, 1]))**2
        if gram_determinant <= 1e-2:
            raise AssertionError(f"repeated eigenspace dimension was lost: {gram_determinant=}")
    return cycles, locked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    parser.add_argument("--cases", type=int, default=54)
    parser.add_argument("--seed", type=int, default=20261104)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error("--cases must be positive")
    for mode in (("checked", "optimized") if args.mode == "both" else (args.mode,)):
        rng = np.random.default_rng(args.seed)
        with tempfile.TemporaryDirectory(prefix="scifort-arnoldi-block-") as tmp:
            exe = build_driver(Path(tmp), mode)
            cycles, locks = [], []
            for i in range(args.cases):
                scenario = i % 9
                n = (24, 40, 56)[(i // 9) % 3]
                a, k, b, target, sigma = make_matrix(rng, n, scenario)
                it, count = test_one(exe, a, k, b, target, sigma, locking=i % 5 != 0)
                cycles.append(it)
                locks.append(count)
            print(f"{mode}: {args.cases} block Arnoldi SciPy dense comparisons PASSED; "
                  f"{sum(x > 1 for x in cycles)} restarted, "
                  f"{sum(x > 0 for x in locks)} Schur-locked, "
                  f"max {max(cycles)} cycles", flush=True)


if __name__ == "__main__":
    main()
