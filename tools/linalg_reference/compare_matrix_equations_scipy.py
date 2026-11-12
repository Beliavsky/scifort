#!/usr/bin/env python3
"""Independent SciPy reference checks for SciFort dense matrix equations.

Builds the standalone Fortran driver with gfortran, runs a stream of random
problems, and compares solutions and true matrix-equation residuals.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import subprocess
import tempfile

import numpy as np
from scipy.linalg import (
    solve_continuous_lyapunov,
    solve_discrete_lyapunov,
    solve_sylvester,
)

ROOT = pathlib.Path(__file__).resolve().parents[2]


def compile_driver(work: pathlib.Path, mode: str) -> pathlib.Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    sources = [
        ROOT / "src/kinds.f90",
        ROOT / "src/linalg.f90",
        ROOT / "src/linalg_matrix_functions.f90",
        ROOT / "src/linalg_matrix_equations.f90",
        ROOT / "tools/linalg_reference/driver_matrix_equations.f90",
    ]
    executable = work / f"matrix_equations_{mode}{'.exe' if os.name == 'nt' else ''}"
    subprocess.run(
        ["gfortran", *flags, "-ffree-line-length-none", *map(str, sources), "-o", str(executable)],
        cwd=work,
        check=True,
        capture_output=True,
        text=True,
    )
    return executable


def formatted_values(array: np.ndarray) -> str:
    vals = array.flatten(order="F")
    if np.iscomplexobj(array):
        return " ".join(f"({x.real:.17e},{x.imag:.17e})" for x in vals)
    return " ".join(f"{x:.17e}" for x in vals)


def random_case(rng: np.random.Generator, code: int, cmplx: bool, max_dim: int):
    m = int(rng.integers(1, max_dim + 1))
    n = int(rng.integers(1, max_dim + 1)) if code == 1 else m

    def rand(shape):
        values = rng.normal(size=shape)
        if cmplx:
            values = values + 1j * rng.normal(size=shape)
        return values

    a = rand((m, m)) * 0.18
    # Keep the equations nonsingular, and the discrete-time system stable.
    if code in (1, 2):
        a += (2.0 if code == 1 else -2.0) * np.eye(m)
    b = None
    if code == 1:
        b = rand((n, n)) * 0.18 + 3.0 * np.eye(n)
    q = rand((m, n))
    if code == 1:
        reference = solve_sylvester(a, b, q)
    elif code == 2:
        reference = solve_continuous_lyapunov(a, q)
    else:
        reference = solve_discrete_lyapunov(a, q)
    return m, n, a, b, q, reference


def run_cases(executable: pathlib.Path, rng: np.random.Generator, count: int, max_dim: int) -> tuple[int, float, float]:
    problems = []
    lines = []
    for idx in range(count):
        code = (idx % 3) + 1
        cmplx = (idx // 3) % 2 == 1
        m, n, a, b, q, reference = random_case(rng, code, cmplx, max_dim)
        lines.append(f"{code} {'C' if cmplx else 'R'} {m} {n}")
        lines.append(formatted_values(a))
        if b is not None:
            lines.append(formatted_values(b))
        lines.append(formatted_values(q))
        problems.append((code, m, n, a, b, q, reference, cmplx))
    proc = subprocess.run([str(executable)], input="\n".join(lines) + "\n", capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(f"Fortran driver returned {proc.returncode}: {proc.stderr[-2000:]}")
    results = iter(proc.stdout.splitlines())
    worst_error = 0.0
    worst_residual = 0.0
    for idx, (code, m, n, a, b, q, ref, cmplx) in enumerate(problems):
        status = int(next(results))
        if status != 0:
            raise AssertionError(f"case {idx}: status={status}")
        entries = []
        for _ in range(m * n):
            numbers = [float(value) for value in next(results).split()]
            entries.append(complex(*numbers) if cmplx else numbers[0])
        x = np.asarray(entries).reshape((m, n), order="F")
        error = np.linalg.norm(x - ref) / max(1.0, np.linalg.norm(ref))
        if code == 1:
            residual = a @ x + x @ b - q
        elif code == 2:
            residual = a @ x + x @ a.conj().T - q
        else:
            residual = x - a @ x @ a.conj().T - q
        relative_residual = np.linalg.norm(residual) / max(1.0, np.linalg.norm(q))
        worst_error = max(worst_error, error)
        worst_residual = max(worst_residual, relative_residual)
        if error > 1e-10 or relative_residual > 1e-10:
            raise AssertionError(f"case {idx}: {error=:.4e}, {relative_residual=:.4e}")
    assert next(results, None) is None, "extra Fortran output"
    return count, worst_error, worst_residual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["checked", "optimized", "both"], default="both")
    parser.add_argument("--cases", type=int, default=360)
    parser.add_argument("--seed", type=int, default=20261106)
    parser.add_argument("--max-dim", type=int, default=6, help="maximum matrix size per dimension")
    args = parser.parse_args()
    if args.cases < 1:
        parser.error("--cases must be positive")
    if not (1 <= args.max_dim <= 20):
        parser.error("--max-dim must be between 1 and 20")
    modes = ["checked", "optimized"] if args.mode == "both" else [args.mode]
    with tempfile.TemporaryDirectory(prefix="scifort-matrix-equations-") as d:
        for mode in modes:
            exe = compile_driver(pathlib.Path(d), mode)
            rng = np.random.default_rng(args.seed)
            count, error, residual = run_cases(exe, rng, args.cases, args.max_dim)
            print(f"{mode}: {count}/{count} passed; max relative error={error:.3e}, "
                  f"max relative residual={residual:.3e}")


if __name__ == "__main__":
    main()
