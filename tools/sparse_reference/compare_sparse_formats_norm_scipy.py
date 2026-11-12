#!/usr/bin/env python3
"""Randomized checks of SciFort DIA/BSR/norm against SciPy 1.17 and NumPy."""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import norm as sparse_norm

ROOT = Path(__file__).resolve().parents[2]


def compile_reference(build: Path, mode: str) -> Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    cc = ["gfortran", "-std=f2018", *flags, "-J", str(build), "-I", str(build)]
    for name, src in (("kinds", "src/kinds.f90"), ("sparse", "src/sparse.f90"),
                      ("lsmr", "src/optimize/lsmr.f90"), ("sparse_linalg", "src/sparse_linalg.f90")):
        subprocess.run([*cc, "-c", str(ROOT / src), "-o", str(build / f"{name}.o")],
                       check=True, cwd=build, capture_output=True, text=True)
    exe = build / "sparse_formats_norm_reference_driver"
    subprocess.run([*cc, str(ROOT / "tools/sparse_reference/sparse_formats_norm_reference_driver.f90"),
                    *(str(build / f"{name}.o") for name in ("kinds", "sparse", "lsmr", "sparse_linalg")),
                    "-o", str(exe)], check=True, cwd=build, capture_output=True, text=True)
    return exe


def compare(mode: str, cases: int, seed: int, prebuilt: Path | None = None) -> None:
    rng = np.random.default_rng(seed)
    inputs = [str(cases)]
    expected = []
    for t in range(cases):
        m = int(rng.integers(2, 10))
        n = int(rng.integers(2, 10))
        br = int(rng.choice([i for i in range(1, m + 1) if m % i == 0]))
        bc = int(rng.choice([i for i in range(1, n + 1) if n % i == 0]))
        nnz = int(rng.integers(0, 3 * (m + n) + 1))
        rows = rng.integers(0, m, size=nnz)
        cols = rng.integers(0, n, size=nnz)
        values = rng.uniform(-5.0, 5.0, size=nnz)
        values[rng.random(nnz) < 0.12] = 0.0
        inputs.append(f"{m} {n} {nnz} {br} {bc}")
        inputs.extend(f"{r+1} {c+1} {v:.17g}" for r, c, v in zip(rows, cols, values))
        mat = sparse.coo_array((values, (rows, cols)), shape=(m, n)).tocsr()
        mat.sum_duplicates()
        dense = mat.toarray()
        expected.append((dense, mat, np.arange(1, n+1) / n, t, br, bc))
    with tempfile.TemporaryDirectory(prefix="scifort-formats-norm-") as temp:
        exe = prebuilt if prebuilt is not None else compile_reference(Path(temp), mode)
        proc = subprocess.run([str(exe)], input="\n".join(inputs) + "\n", check=True,
                              capture_output=True, text=True)
    lines = proc.stdout.splitlines()
    if len(lines) != 7 * cases:
        raise AssertionError(f"Expected {7*cases} output lines, got {len(lines)}; {proc.stderr}")
    count = 0
    max_error = 0.0
    for dense, mat, x, t, br, bc in expected:
        chunk = lines[7*t:7*t+7]
        parsed = {}
        for line, tag in zip(chunk, "NCRDBYZ"):
            pieces = line.split()
            if pieces[0] != tag:
                raise AssertionError(f"case {t}: expected {tag}, got {pieces[0]}")
            parsed[tag] = np.array([float(z) for z in pieces[1:]], dtype=float)
        spectral = np.linalg.norm(dense, ord=2)
        norms = np.array([sparse_norm(mat), sparse_norm(mat, ord=1), sparse_norm(mat, ord=-1),
                          sparse_norm(mat, ord=np.inf), sparse_norm(mat, ord=-np.inf)])
        np.testing.assert_allclose(parsed["N"][:5], norms, rtol=2e-12, atol=2e-12)
        if parsed["N"][5] >= 0:
            np.testing.assert_allclose(parsed["N"][5], spectral, rtol=1e-7, atol=1e-8)
            count += 1
        np.testing.assert_allclose(parsed["C"], np.linalg.norm(dense, axis=0), atol=2e-12)
        np.testing.assert_allclose(parsed["R"], np.linalg.norm(dense, ord=1, axis=1), atol=2e-12)
        for tag in ("D", "B"):
            np.testing.assert_allclose(parsed[tag], dense.ravel(order="F"), rtol=2e-12, atol=2e-12)
        for tag in ("Y", "Z"):
            np.testing.assert_allclose(parsed[tag], mat @ x, rtol=2e-12, atol=2e-12)
        if t < 16:
            # Check SciPy's DIA and BSR conversions explicitly, including odd blocks.
            np.testing.assert_allclose(sparse.dia_matrix(mat).toarray(), dense, atol=2e-12)
            np.testing.assert_allclose(mat.tobsr(blocksize=(br, bc)).toarray(), dense, atol=2e-12)
        max_error = max(max_error, float(np.max(np.abs(parsed["D"] - dense.ravel(order="F")))))
    print(f"{mode}: {cases} randomized CSR/DIA/BSR/norm cases passed; "
          f"{count} successful spectral norm comparisons, max DIA error={max_error:.2g}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("both", "checked", "optimized"), default="both")
    parser.add_argument("--cases", type=int, default=150)
    parser.add_argument("--seed", type=int, default=20261105)
    parser.add_argument("--executable", type=Path, help="Optional prebuilt Fortran comparison driver")
    args = parser.parse_args()
    for mode in (("checked", "optimized") if args.mode == "both" else (args.mode,)):
        compare(mode, args.cases, args.seed, args.executable)


if __name__ == "__main__":
    main()
