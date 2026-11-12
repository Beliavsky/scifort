#!/usr/bin/env python3
"""Independent SciPy reference tests for scifort_linalg_qr_updates.

Run with a compiled driver::

    python validation/compare_linalg_qr_updates.py --driver /path/to/driver

Alternatively omit --driver to build a checked GNU Fortran driver locally.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import scipy
from scipy import linalg


@dataclass
class Case:
    kind: int
    operation: int
    original: np.ndarray
    q: np.ndarray
    r: np.ndarray
    u: np.ndarray
    v: np.ndarray
    k: int  # one-based
    expected: np.ndarray
    scipy_shapes: tuple[tuple[int, int], tuple[int, int]]


def compiler_driver(directory: Path, optimized: bool = False) -> Path:
    """Compile only the SciFort QR dependency closure; no external packages."""
    root = Path(__file__).resolve().parents[1]
    directory.mkdir(parents=True, exist_ok=True)
    names = ["kinds.f90", "linalg.f90", "linalg_svd.f90", "linalg_qr.f90", "linalg_qr_updates.f90"]
    flags = ["-std=f2008", "-O3"] if optimized else ["-std=f2008", "-O0", "-g", "-fcheck=all", "-fbacktrace"]
    driver = directory / "linalg_qr_updates_driver"
    cmd = ["gfortran", *flags, *(str(root / "src" / n) for n in names),
           str(root / "reference" / "linalg_qr_updates_driver.f90"), "-o", str(driver)]
    subprocess.run(cmd, cwd=directory, check=True, capture_output=True, text=True)
    return driver


def generate(seed: int = 20261109, reps: int = 2) -> list[Case]:
    rng = np.random.default_rng(seed)
    cases = []
    for dtype in (0, 1):
        for m, n in ((1, 1), (1, 4), (4, 1), (2, 3), (3, 2), (4, 6),
                     (6, 4), (5, 5), (9, 6), (7, 11), (12, 8)):
            for full in (False, True):
                for op in range(1, 9):
                    for rep in range(reps):
                        if (op == 7 and m <= 1) or (op == 8 and n <= 1):
                            continue  # SciPy has nonorthonormal empty factor layouts
                        def draw(shape):
                            x = rng.standard_normal(shape)
                            if dtype:
                                x = x + 1j * rng.standard_normal(shape)
                            return x
                        a = draw((m, n))
                        q, r = linalg.qr(a, mode="full" if full else "economic")
                        p = int(rng.integers(1, min(4, max(m, n)) + 1))
                        if op == 1:
                            p = 1
                        if op == 2:
                            p = min(p,max(1,min(m,n)-1))
                        if op in (3, 5):
                            p = 1
                        if op == 7:
                            p = min(p, m - 1)
                        if op == 8:
                            p = min(p, n - 1)
                        k = 1
                        u = np.empty((0, 0), dtype=a.dtype)
                        v = np.empty((0, 0), dtype=a.dtype)
                        try:
                            if op in (1, 2):
                                u = draw((m, p))
                                v = draw((n, p))
                                if rep == 1 and op == 2:  # include updates in span(Q)
                                    u[:, 0] = q[:, 0]
                                uv = np.outer(u[:, 0], v[:, 0].conj()) if op == 1 else u @ v.conj().T
                                expected = a + uv
                                qs, rs = linalg.qr_update(q, r, u[:, 0], v[:, 0]) if op == 1 else linalg.qr_update(q, r, u, v)
                            elif op in (3, 4):
                                u = draw((p, n))
                                k = int(rng.integers(1, m + 2))
                                # np.insert with 2-D data has different axis semantics: use explicit concatenation.
                                expected = np.concatenate((a[:k-1], u, a[k-1:]), axis=0)
                                qs, rs = linalg.qr_insert(q, r, u[0] if op == 3 else u, k-1, which="row")
                            elif op in (5, 6):
                                u = draw((m, p))
                                k = int(rng.integers(1, n + 2))
                                expected = np.concatenate((a[:, :k-1], u, a[:, k-1:]), axis=1)
                                qs, rs = linalg.qr_insert(q, r, u[:, 0] if op == 5 else u, k-1, which="col")
                            elif op == 7:
                                k = int(rng.integers(1, m-p+2))
                                expected = np.concatenate((a[:k-1], a[k-1+p:]), axis=0)
                                qs, rs = linalg.qr_delete(q, r, k-1, p=p, which="row")
                            else:
                                k = int(rng.integers(1, n-p+2))
                                expected = np.concatenate((a[:, :k-1], a[:, k-1+p:]), axis=1)
                                qs, rs = linalg.qr_delete(q, r, k-1, p=p, which="col")
                        except (ValueError, linalg.LinAlgError) as e:
                            # SciPy may reject a thin-column insertion whose new vector is
                            # linearly dependent on Q; see the documented limitation.
                            print(f"skipping SciPy-rejected case m={m} n={n} op={op}: {str(e)[:72]}")
                            continue
                        cases.append(Case(dtype, op, a, q, r, u, v, k, expected, (qs.shape,rs.shape)))
    return cases


def to_text(case: Case) -> list[str]:
    m, n = case.original.shape
    h = case.q.shape[1]
    p = case.u.shape[1] if case.operation in (1, 2, 5, 6) else case.u.shape[0]
    if case.operation in (7, 8):
        p = (m - case.expected.shape[0]) if case.operation == 7 else (n - case.expected.shape[1])
    lines = [f"{case.kind} {case.operation} {m} {n} {h} {p} {case.k}"]
    for arr in (case.q, case.r, case.u, case.v):
        for value in arr.ravel(order="F"):
            lines.append(f"{value.real:.17e} {value.imag:.17e}")
    return lines


def compare(driver: Path, cases: list[Case]) -> dict:
    lines = [str(len(cases))]
    for case in cases:
        lines.extend(to_text(case))
    proc = subprocess.run([str(driver)], input="\n".join(lines)+"\n", text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, timeout=180)
    cursor = iter(proc.stdout.splitlines())
    worst_reconstruction = 0.0
    worst_unitary = 0.0
    for caseidx, case in enumerate(cases):
        status, m2, n2, h2 = map(int, next(cursor).split())
        if status != 0:
            raise AssertionError(f"case {caseidx}: Fortran status={status}, {case}")
        assert (m2, n2) == case.expected.shape
        assert ((m2, h2), (h2, n2)) == case.scipy_shapes, (caseidx, case.scipy_shapes, ((m2,h2),(h2,n2)))
        def read_matrix(rows: int, cols: int) -> np.ndarray:
            values = [complex(*map(float, next(cursor).split())) for _ in range(rows*cols)]
            return np.array(values).reshape((rows,cols),order="F")
        q = read_matrix(m2, h2)
        r = read_matrix(h2, n2)
        if case.kind == 0:
            assert np.max(abs(q.imag),initial=0) == 0
            assert np.max(abs(r.imag),initial=0) == 0
        if not np.all(np.isfinite(q)) or not np.all(np.isfinite(r)):
            raise AssertionError(f"case {caseidx}: nonfinite factor")
        scale = max(1.0, np.max(abs(case.expected), initial=0))
        error = np.max(abs(q@r - case.expected), initial=0)/scale
        unitary = np.max(abs(q.conj().T @ q - np.eye(h2)), initial=0)
        triangular = np.max(abs(np.tril(r, -1)), initial=0)
        worst_reconstruction = max(worst_reconstruction, error)
        worst_unitary = max(worst_unitary, unitary)
        if error > 5e-10 or unitary > 5e-10 or triangular > 5e-10:
            raise AssertionError(f"case {caseidx}: type={case.kind},op={case.operation}, "
                                 f"shape={case.original.shape}, bad QR "
                                 f"reconstruction={error}, unitary={unitary}, triangular={triangular}")
    return {"cases": len(cases), "max_reconstruction": worst_reconstruction,
            "max_orthogonality": worst_unitary, "scipy": scipy.__version__}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, help="previously compiled Fortran reference driver")
    parser.add_argument("--optimized", action="store_true", help="use -O3 for self-built driver")
    parser.add_argument("--reps", type=int, default=2, help="repetitions per shape/type/operation")
    args = parser.parse_args()
    cases = generate(reps=args.reps)
    if args.driver:
        result = compare(args.driver.resolve(), cases)
    else:
        with tempfile.TemporaryDirectory(prefix="scifort_qr_update_") as tmp:
            result = compare(compiler_driver(Path(tmp), args.optimized), cases)
    print("PASS:", result)


if __name__ == "__main__":
    main()
