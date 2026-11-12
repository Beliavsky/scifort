#!/usr/bin/env python3
"""Compare Schur-based SciFort Sylvester/Lyapunov solvers against SciPy.

Fortran driver receives a streamed collection of problems, avoiding per-case
process startup costs. Test families include random, triangular, clustered,
normal and nonnormal coefficient matrices. No BLAS/LAPACK linked in Fortran.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import os
import subprocess
import tempfile

# The randomized matrices are at most order 40. Threaded BLAS overhead can
# dominate their tiny reference solves on high-core-count machines, so use a
# single thread for deterministic, bounded SciPy reference-test execution.
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_variable] = "1"

import numpy as np
from scipy.linalg import solve_sylvester, solve_continuous_lyapunov, solve_discrete_lyapunov

ROOT = Path(__file__).resolve().parents[2]


def compile_driver(work: Path, mode: str) -> Path:
    flags = ["-O1", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
    source_paths = [
        "src/kinds.f90", "src/linalg.f90", "src/linalg_matrix_functions.f90",
        "src/linalg_matrix_equations.f90",
        "tools/linalg_reference/driver_matrix_equations_schur.f90",
    ]
    exe = work / ("driver_schur_" + mode + (".exe" if os.name == "nt" else ""))
    proc = subprocess.run(["gfortran", "-std=f2018", *flags, *map(lambda s: str(ROOT / s), source_paths),
                            "-o", str(exe)], cwd=work, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(proc.stderr[-4000:])
    return exe


def write_values(a: np.ndarray) -> str:
    if np.iscomplexobj(a):
        return " ".join(f"({z.real:.17e},{z.imag:.17e})" for z in a.ravel(order="F"))
    return " ".join(f"{z:.17e}" for z in a.ravel(order="F"))


def make_case(rng: np.random.Generator, idx: int, max_order: int, large: bool):
    code = idx % 3 + 1
    complex_input = (idx // 3) % 2 == 1
    family = (idx // 6) % 5
    low = max(21, max_order // 2) if large and max_order > 20 else 1
    m = int(rng.integers(low, max_order + 1))
    n = int(rng.integers(low, max_order + 1)) if code == 1 else m

    def sample(shape):
        z = rng.normal(size=shape)
        if complex_input:
            z = z + 1j*rng.normal(size=shape)
        return z

    def matrix(order: int, shift: float, magnitude: float):
        r = sample((order, order))
        if family == 0:  # shifted Gaussian (general, mostly well conditioned)
            return shift*np.eye(order) + magnitude*r/max(1, np.sqrt(order))
        if family == 1:  # upper triangular, repeated diagonal
            result = magnitude*np.triu(r, 1)/max(1, np.sqrt(order))
            result[np.diag_indices(order)] = shift
            return result
        if family == 2:  # bidiagonal Jordan-like (moderate nonnormality)
            result = shift*np.eye(order, dtype=r.dtype)
            result += np.diag(np.full(order-1, magnitude), 1)
            return result
        if family == 3:  # normal with randomly rotated spectrum
            h, _ = np.linalg.qr(sample((order, order)))
            diag = shift + magnitude*rng.uniform(-0.2, 0.2, size=order)
            return (h * diag) @ h.conj().T
        # shifted, modestly nonnormal triangular matrix
        return shift*np.eye(order) + magnitude*np.triu(r)/max(1, np.sqrt(order))

    shift = 2.0 if code == 1 else (-2.0 if code == 2 else 0.05)
    magnitude = 0.4 if code != 3 else 0.2
    a = matrix(m, shift, magnitude)
    b = matrix(n, 3.0, magnitude) if code == 1 else None
    q = sample((m, n))
    if code == 1:
        reference = solve_sylvester(a, b, q)
    elif code == 2:
        reference = solve_continuous_lyapunov(a, q)
    else:
        reference = solve_discrete_lyapunov(a, q)
    return (code, complex_input, m, n, a, b, q, reference, family)


def test_mode(exe: Path, rng: np.random.Generator, cases: int, max_order: int, large: bool):
    records = [make_case(rng, i, max_order, large) for i in range(cases)]
    lines = []
    for code, complex_input, m, n, a, b, q, _, _ in records:
        lines.append(f"{code} {'C' if complex_input else 'R'} {m} {n}")
        lines.append(write_values(a))
        if b is not None:
            lines.append(write_values(b))
        lines.append(write_values(q))
    proc = subprocess.run([str(exe)], input="\n".join(lines) + "\n", capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(f"driver failed with {proc.returncode}: {proc.stderr[-2000:]}")
    iterator = iter(proc.stdout.splitlines())
    worst_diff = 0.0
    worst_backward = 0.0
    for idx, (code, complex_input, m, n, a, b, q, reference, family) in enumerate(records):
        try:
            status = int(next(iterator))
            entries = [tuple(map(float, next(iterator).split())) for _ in range(m*n)]
        except (StopIteration, ValueError) as exc:
            raise AssertionError(f"case {idx}: invalid driver output") from exc
        if status:
            raise AssertionError(f"case {idx}: nonzero Fortran status {status}, code={code}, shape={m}x{n}, family={family}")
        values = [complex(*v) if complex_input else v[0] for v in entries]
        x = np.array(values).reshape((m, n), order="F")
        diff = np.linalg.norm(x-reference)/max(1.0, np.linalg.norm(reference))
        if code == 1:
            residual = a @ x + x @ b - q
            denom = (np.linalg.norm(a) + np.linalg.norm(b))*np.linalg.norm(x) + np.linalg.norm(q)
        elif code == 2:
            residual = a @ x + x @ a.conj().T - q
            denom = 2*np.linalg.norm(a)*np.linalg.norm(x) + np.linalg.norm(q)
        else:
            residual = x - a @ x @ a.conj().T - q
            denom = (1.0 + np.linalg.norm(a)**2)*np.linalg.norm(x) + np.linalg.norm(q)
        backward = np.linalg.norm(residual)/max(np.finfo(float).tiny, denom)
        worst_diff = max(worst_diff, diff)
        worst_backward = max(worst_backward, backward)
        if not (np.isfinite(diff) and np.isfinite(backward)) or diff > 3e-8 or backward > 3e-9:
            raise AssertionError(f"case {idx}: code={code}, {m}x{n}, family={family}, diff={diff:.3e}, residual={backward:.3e}")
    if next(iterator, None) is not None:
        raise AssertionError("extra driver output")
    return worst_diff, worst_backward


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    p.add_argument("--cases", type=int, default=240)
    p.add_argument("--max-order", type=int, default=32)
    p.add_argument("--seed", type=int, default=711)
    p.add_argument("--large", action="store_true", help="concentrate cases at orders 21 through max-order")
    args = p.parse_args()
    if args.cases < 1 or not 1 <= args.max_order <= 40:
        p.error("requires cases >= 1 and 1 <= max-order <= 40")
    modes = ["checked", "optimized"] if args.mode == "both" else [args.mode]
    with tempfile.TemporaryDirectory(prefix="scifort-schur-") as d:
        for mode in modes:
            exe = compile_driver(Path(d), mode)
            diff, backward = test_mode(exe, np.random.default_rng(args.seed),
                                       args.cases, args.max_order, args.large)
            print(f"{mode}: {args.cases}/{args.cases} passed; max rel diff={diff:.3e}; max backward error={backward:.3e}")


if __name__ == "__main__":
    main()
