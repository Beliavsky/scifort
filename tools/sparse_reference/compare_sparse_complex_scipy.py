#!/usr/bin/env python3
"""Randomized independent SciPy/NumPy checks for complex CSR/CSC + Krylov.

Usage: python compare_sparse_complex_scipy.py path/to/complex_sparse_reference_driver
The executable processes cases from stdin, using a stable one-based COO format.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sla


@dataclass
class Case:
    dense: np.ndarray
    rows: np.ndarray
    cols: np.ndarray
    vals: np.ndarray
    b: np.ndarray
    x0: np.ndarray
    solver: int
    fmt: int
    jacobi: bool
    restart: int
    cap: int


def draw_complex(rng, size):
    return rng.normal(size=size) + 1j * rng.normal(size=size)


def make_case(rng, solver, n, fmt, jacobi, restart):
    if solver == 0:
        m = int(rng.integers(1, 16))
        raw = draw_complex(rng, (n, m))
        raw[rng.random((n, m)) < 0.65] = 0
        cap = 0
    elif solver == 1:
        m = n
        # B = D+R, A=B^H B + constant I is Hermitian SPD.
        raw = draw_complex(rng, (n, n)) * (rng.random((n, n)) < 0.13)
        raw = raw + (2.0 + rng.uniform(0.5, 2.0, n)) * np.eye(n)
        raw = raw.conj().T @ raw + np.eye(n)
        cap = max(30, n * 7)
    else:
        m = n
        raw = draw_complex(rng, (n, n)) * (rng.random((n, n)) < 0.18)
        raw += np.diag(n + 3.0 + rng.uniform(0.5, 1.5, n) + 0.5j)
        cap = max(50, n * 10)

    # Split some entries into duplicates and retain explicitly stored zeros.
    rows, cols = np.nonzero(raw)
    rows, cols = rows.tolist(), cols.tolist()
    vals = [raw[i, j] for i, j in zip(rows, cols)]
    if vals:
        for _ in range(min(8, len(vals))):
            k = int(rng.integers(0, len(vals)))
            val = vals[k]
            vals[k] = val * 0.4
            rows.append(rows[k])
            cols.append(cols[k])
            vals.append(val * 0.6)
    for _ in range(int(rng.integers(0, 4))):
        rows.append(int(rng.integers(n)))
        cols.append(int(rng.integers(m)))
        vals.append(0j)
    order = rng.permutation(len(vals))
    rows, cols = np.asarray(rows, dtype=np.int64)[order], np.asarray(cols, dtype=np.int64)[order]
    vals = np.asarray(vals, dtype=np.complex128)[order]
    sp = sparse.coo_array((vals, (rows, cols)), shape=(n, m)).tocsr()
    d = sp.toarray()
    # Explicit x0 must be nonzero in about half the cases.
    x0 = draw_complex(rng, m) if rng.random() < 0.6 else np.zeros(m, complex)
    b = draw_complex(rng, n) if solver == 0 else sp @ draw_complex(rng, m)
    return Case(d, rows, cols, vals, b, x0, solver, fmt, bool(jacobi), restart, cap)


def pack(case):
    n, m = case.dense.shape
    lines = [f"{n} {m} {len(case.vals)} {case.solver} {case.fmt} {int(case.jacobi)} {case.restart} {case.cap}"]
    for row, col, value in zip(case.rows, case.cols, case.vals):
        lines.append(f"{row+1} {col+1} {value.real:.17e} {value.imag:.17e}")
    for array in (case.b, case.x0):
        for value in array:
            lines.append(f"{value.real:.17e} {value.imag:.17e}")
    return "\n".join(lines) + "\n"


def compare(executable, seed, per_solver):
    rng = np.random.default_rng(seed)
    cases = []
    for solver in range(4):
        for i in range(per_solver):
            n = int(rng.integers(1, 21))
            fmt = 1 + i % 2
            jacobi = bool(i % 3)
            restart = min(n, max(1, int(rng.integers(1, 7))))
            cases.append(make_case(rng, solver, n, fmt, jacobi, restart))
    proc = subprocess.run([executable], input="".join(map(pack, cases)), text=True,
                          capture_output=True, check=True, timeout=120)
    data = iter(proc.stdout.splitlines())
    worst = 0.0
    for index, case in enumerate(cases):
        try:
            status_s, steps_s, residual_s = next(data).split()
        except (ValueError, StopIteration) as exc:
            raise AssertionError(f"bad output header at case {index}: {proc.stdout[-500:]}") from exc
        status = int(status_s)
        steps = int(steps_s)
        residual = float(residual_s)
        count = case.dense.shape[0]
        solution = np.array([complex(*map(float, next(data).split())) for _ in range(count)])
        if status != 0:
            raise AssertionError(f"case {index}: solver {case.solver} returned status {status}, steps {steps}, "
                                 f"residual {residual}, shape {case.dense.shape}, restart={case.restart}")
        if case.solver == 0:
            expected = sparse.coo_array((case.vals, (case.rows, case.cols)),
                                        shape=case.dense.shape) @ case.x0
        else:
            expected = np.linalg.solve(case.dense, case.b)
            a = sparse.csr_array(case.dense)
            if case.solver == 1:
                reference, rc = sla.cg(a, case.b, x0=case.x0, atol=1e-12, rtol=1e-10,
                                       maxiter=case.cap)
            elif case.solver == 2:
                reference, rc = sla.bicgstab(a, case.b, x0=case.x0, atol=1e-12, rtol=1e-10,
                                             maxiter=case.cap)
            else:
                reference, rc = sla.gmres(a, case.b, x0=case.x0, atol=1e-12, rtol=1e-10,
                                          restart=case.restart, maxiter=case.cap)
            if rc != 0 or not np.allclose(reference, expected, rtol=5e-9, atol=2e-9):
                raise AssertionError(f"case {index}: SciPy reference failed: {rc}")
        err = float(np.linalg.norm(solution - expected))
        worst = max(worst, err)
        if not np.allclose(solution, expected, rtol=5e-9, atol=2e-9):
            raise AssertionError(f"case {index}: wrong result, error {err:.3g}, status={status}, "
                                 f"shape={case.dense.shape}, fmt={case.fmt}, solver={case.solver}")
        true_residual = (float(np.linalg.norm(case.b - case.dense @ solution))
                         if case.solver > 0 else 0.0)
        if case.solver > 0 and not np.isclose(residual, true_residual, rtol=3e-8, atol=2e-11):
            raise AssertionError(f"case {index}: reported residual {residual}, true {true_residual}")
    if next(data, None) is not None:
        raise AssertionError("unexpected trailing Fortran output")
    print(f"PASS {len(cases)} randomized complex sparse comparisons (seed={seed}), "
          f"max Euclidean error={worst:.3g}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('--seed', type=int, default=20261105)
    parser.add_argument('--per-solver', type=int, default=100)
    args = parser.parse_args()
    compare(args.executable, args.seed, args.per_solver)


if __name__ == '__main__':
    main()
