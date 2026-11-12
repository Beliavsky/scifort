#!/usr/bin/env python3
"""Independent NumPy spectra for matrix-free/cached complex SciFort Arnoldi.

Example:
    python tools/sparse_reference/compare_sparse_complex_eigs_operator_numpy.py --mode both --cases 90 --scipy-samples 36
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
DRIVER = ROOT / "tools/sparse_reference/complex_eigs_operator_driver.f90"


def compile_driver(path: Path, mode: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace"]
             if mode == "checked" else ["-O3"])
    objects = []
    for name in ["kinds", "sparse", "sparse_coo", "sparse_coo_nd",
                 "sparse_coo_nd_complex", "sparse_complex",
                 "sparse_linalg_complex_direct", "sparse_linalg_complex_eigs"]:
        obj = path / (name + ".o")
        subprocess.run(["gfortran", "-std=f2018", *flags, "-J", str(path),
                        "-I", str(path), "-c", str(SRC / (name + ".f90")),
                        "-o", str(obj)], cwd=path, check=True, capture_output=True)
        objects.append(str(obj))
    exe = path / "operator_reference"
    subprocess.run(["gfortran", "-std=f2018", *flags, "-J", str(path),
                    "-I", str(path), str(DRIVER), *objects, "-o", str(exe)],
                   cwd=path, check=True, capture_output=True)
    return exe


def make_cases(cases: int, seed: int):
    rng = np.random.default_rng(seed)
    specs = []
    lines = [str(cases)]
    for j in range(cases):
        n = int(rng.integers(6, 16))
        k = int(rng.integers(1, min(n - 2, 4) + 1))
        shifted = int(j % 2 == 1)
        which = "LM" if shifted else ["LM", "SM", "LR", "SR", "LI", "SI"][j % 6]
        # Monotonic real part but variable imaginary component, no ties.
        diagonal = np.arange(1, n + 1, dtype=float) + 0.06 * rng.normal(size=n)
        diagonal = diagonal + 1j * (0.12 * np.arange(n) + rng.normal(size=n) * 0.04)
        upper = 0.1 * (rng.normal(size=n - 1) + 1j * rng.normal(size=n - 1))
        middle = int(rng.integers(1, n-2))
        sigma = diagonal[middle] + (0.13 + 0.05j)
        specs.append((n, k, shifted, which, diagonal, upper, sigma))
        lines.append(f"{n} {k} {shifted} {which}")
        lines.extend(f"{z.real:.17e} {z.imag:.17e}" for z in diagonal)
        lines.extend(f"{z.real:.17e} {z.imag:.17e}" for z in upper)
        lines.append(f"{sigma.real:.17e} {sigma.imag:.17e}")
    return specs, "\n".join(lines) + "\n"


def eigen_priority(values: np.ndarray, which: str):
    return {
        "LM": np.abs(values), "SM": -np.abs(values),
        "LR": values.real, "SR": -values.real,
        "LI": values.imag, "SI": -values.imag,
    }[which]


def check(exe: Path, cases: int, seed: int, mode: str, scipy_samples: int = 0):
    specs, data = make_cases(cases, seed)
    run = subprocess.run([str(exe)], input=data, text=True,
                         capture_output=True, check=True, timeout=150)
    rows = run.stdout.splitlines()
    index = 0
    max_err = 0.0
    compared_scipy = 0
    for case_id, (n, k, shifted, which, diag, upper, sigma) in enumerate(specs, 1):
        assert rows[index].strip() == f"CASE {case_id}", (mode, case_id, rows[index])
        index += 1
        info, cached_info, products, num = map(int, rows[index].split())
        index += 1
        assert info == cached_info == 0 and num == k and products >= n, (
            mode, case_id, info, cached_info, num, products)
        block = np.array([[float(x) for x in rows[index + j].split()] for j in range(k)])
        index += k
        got = block[:, 0] + 1j * block[:, 1]
        got_cache = block[:, 3] + 1j * block[:, 4]
        if shifted:
            transformed = 1.0 / (diag - sigma)
            order = np.argsort(-eigen_priority(transformed, which))
        else:
            order = np.argsort(-eigen_priority(diag, which))
        truth = diag[order[:k]]
        # Sorted eigenvalues and true residuals: independently validated.
        error = max(np.max(np.abs(got - truth)), np.max(np.abs(got_cache - truth)))
        max_err = max(max_err, error)
        assert error < 2.0e-6, (mode, case_id, which, truth, got, got_cache)
        assert np.max(block[:, [2, 5]]) < 2.0e-6, (mode, case_id, block)
        if shifted and compared_scipy < scipy_samples:
            from scipy.sparse import diags
            from scipy.sparse.linalg import eigs
            from scipy.optimize import linear_sum_assignment
            a = diags([diag, upper], [0, 1], shape=(n, n), format="csc")
            scipy_values = eigs(a, k=k, sigma=sigma, which=which, tol=1.0e-12)[0]
            cost = np.abs(scipy_values[:, None] - got[None, :])
            rr, cc = linear_sum_assignment(cost)
            assert np.max(cost[rr, cc]) < 2.0e-6, (
                mode, case_id, "SciPy mismatch", scipy_values, got)
            compared_scipy += 1
    assert index == len(rows), f"Unexpected trailing output: {rows[index:]}"
    print(f"{mode}: {cases}/{cases} matrix-free and cached cases pass; "
          f"max eigenvalue error {max_err:.3e}; "
          f"direct SciPy checks {compared_scipy}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["checked", "optimized", "both"], default="both")
    parser.add_argument("--cases", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261106)
    parser.add_argument("--scipy-samples", type=int, default=0,
                        help="Number of shifted cases to verify directly with SciPy eigs")
    args = parser.parse_args()
    for mode in (["checked", "optimized"] if args.mode == "both" else [args.mode]):
        with tempfile.TemporaryDirectory(prefix="scifort-eigs-operator-") as dirname:
            exe = compile_driver(Path(dirname), mode)
            check(exe, args.cases, args.seed, mode, args.scipy_samples)


if __name__ == "__main__":
    main()
