#!/usr/bin/env python3
"""Randomized SciPy CARE/DARE references for SciFort's dense Riccati solvers.

Runs native Fortran in a single stream, checks returned status, compares
stabilizing solutions, and independently evaluates the original generalized
matrix equations.  Uses no SciPy code in the Fortran implementation.
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import tempfile

import numpy as np
from scipy.linalg import solve_continuous_are, solve_discrete_are, LinAlgError

ROOT = pathlib.Path(__file__).resolve().parents[2]


class ReferenceInaccurate(Exception):
    """SciPy returned a matrix not satisfying its documented equation."""


def compile_driver(work: pathlib.Path, mode: str) -> pathlib.Path:
    flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    files = [
        ROOT / "src/kinds.f90", ROOT / "src/linalg.f90",
        ROOT / "src/linalg_matrix_functions.f90",
        ROOT / "src/linalg_matrix_equations.f90", ROOT / "src/linalg_riccati.f90",
        ROOT / "tools/linalg_reference/driver_riccati.f90",
    ]
    output = work / f"riccati_{mode}"
    subprocess.run(["gfortran", *flags, "-ffree-line-length-none", *map(str, files),
                    "-o", str(output)], check=True, cwd=work, capture_output=True, text=True)
    return output


def matrix_text(arr: np.ndarray, complex_: bool) -> str:
    vals = arr.flatten(order="F")
    if complex_:
        return " ".join(f"({float(v.real):.17e},{float(v.imag):.17e})" for v in vals) + "\n"
    return " ".join(f"{float(v):.17e}" for v in vals) + "\n"


def case(rng: np.random.Generator, idx: int, max_n: int):
    n = int(rng.integers(1, max_n + 1))
    m = n if idx % 5 else int(rng.integers(1, n + 1))
    complex_ = idx % 2 == 1
    method = idx // 2 % 2 + 1
    generalized = (idx // 4) % 3 != 0
    cross = (idx // 7) % 3 != 0

    def rand(shape):
        a = rng.normal(size=shape)
        if complex_:
            a = a + 1j * rng.normal(size=shape)
        return a

    a = rand((n, n)) * 0.22
    if method == 1:
        a += [-1.0, 1.0, 0.55, -2.0][(idx // 8) % 4] * np.eye(n)
    else:
        a += [0.25, 1.3, -0.9, 0.75][(idx // 8) % 4] * np.eye(n)
    b = rand((n, m)) * 0.5
    if m == n:
        b += np.eye(n)
    q0 = rand((n, n))
    r0 = rand((m, m))
    q = q0.conj().T @ q0 + 0.45 * np.eye(n)
    r = r0.conj().T @ r0 + 0.9 * np.eye(m)
    e = np.eye(n) + (0.02 * rand((n, n)) if generalized else 0.0)
    s = 0.01 * rand((n, m)) if cross else np.zeros((n, m))
    dtype = np.complex128 if complex_ else np.float64
    a, b, q, r, e, s = [np.asarray(arr, dtype=dtype) for arr in (a, b, q, r, e, s)]
    solve = solve_continuous_are if method == 1 else solve_discrete_are
    # SciPy 1.17 default balancing can give a poor generalized residual;
    # unbalanced QZ is independently checked below.
    ref = solve(a, b, q, r, e=e, s=s, balanced=False)
    record = method, complex_, n, m, a, b, q, r, e, s, ref
    if residual(record, ref) > 2.0e-8:
        raise ReferenceInaccurate("SciPy reference failed independent residual verification")
    return record


def residual(data, x):
    method, _, _, _, a, b, q, r, e, s, _ = data
    if method == 1:
        lhs = e.conj().T @ x @ a + a.conj().T @ x @ e + q
        left = e.conj().T @ x @ b + s
        right = b.conj().T @ x @ e + s.conj().T
        k = np.linalg.solve(r, right)
    else:
        lhs = a.conj().T @ x @ a - e.conj().T @ x @ e + q
        left = a.conj().T @ x @ b + s
        right = b.conj().T @ x @ a + s.conj().T
        k = np.linalg.solve(r + b.conj().T @ x @ b, right)
    return np.linalg.norm(lhs - left @ k) / max(1.0, np.linalg.norm(q), np.linalg.norm(lhs))


def run(exe: pathlib.Path, rng: np.random.Generator, count: int, max_n: int):
    data = []
    skipped = 0
    inaccurate = 0
    attempt = 0
    while len(data) < count:
        try:
            data.append(case(rng, attempt, max_n))
        except LinAlgError:
            # Hamiltonian/symplectic pencil not isolated by SciPy.
            skipped += 1
        except ReferenceInaccurate:
            inaccurate += 1
        if skipped + inaccurate > max(20, count):
            raise AssertionError("Too many SciPy reference failures")
        attempt += 1
    inp = []
    for method, complex_, n, m, a, b, q, r, e, s, _ in data:
        inp.append(f"{method} {'C' if complex_ else 'R'} {n} {m}\n")
        inp.extend(matrix_text(arr, complex_) for arr in (a, b, q, r, e, s))
    proc = subprocess.run([str(exe)], input="".join(inp), text=True,
                          capture_output=True, timeout=180)
    if proc.returncode:
        raise RuntimeError(proc.stderr[-5000:])
    lines = iter(proc.stdout.splitlines())
    max_error = 0.0
    max_res = 0.0
    for i, record in enumerate(data):
        status, steps, reported_res = next(lines).split()
        n = record[2]
        got = np.asarray([complex(*map(float, next(lines).split()))
                          for _ in range(n*n)]).reshape((n, n), order="F")
        ref = record[-1]
        if int(status) != 0:
            raise AssertionError(f"case {i} {record[:4]} status={status}, steps={steps}, residual={reported_res}")
        err = np.linalg.norm(got - ref) / max(1.0, np.linalg.norm(ref))
        res = residual(record, got)
        max_error = max(max_error, err)
        max_res = max(max_res, res)
        if err > 2.0e-7 or res > 2.0e-8 or abs(float(reported_res) - res) > 2.0e-8:
            raise AssertionError(f"case {i} {record[:4]} error={err:.4e}, residual={res:.4e}, steps={steps}")
    if next(lines, None) is not None:
        raise AssertionError("Unexpected extra output")
    return max_error, max_res, skipped, inaccurate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["both", "checked", "optimized"], default="both")
    parser.add_argument("--cases", type=int, default=360)
    parser.add_argument("--max-n", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261106)
    opts = parser.parse_args()
    for mode in (["checked", "optimized"] if opts.mode == "both" else [opts.mode]):
        with tempfile.TemporaryDirectory(prefix="riccati_") as path:
            exe = compile_driver(pathlib.Path(path), mode)
            maxerr, maxres, skipped, inaccurate = run(exe, np.random.default_rng(opts.seed), opts.cases, opts.max_n)
            print(f"{mode}: {opts.cases} CARE/DARE comparisons passed; "
                  f"max solution error {maxerr:.3e}; max residual {maxres:.3e}; "
                  f"SciPy raised {skipped}; SciPy residual failures {inaccurate}")


if __name__ == "__main__":
    main()
