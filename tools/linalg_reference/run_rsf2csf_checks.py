#!/usr/bin/env python3
"""Compare SciFort rsf2csf with scipy.linalg.rsf2csf and Schur invariants.

Tests input real Schur factors computed by SciPy and, separately, factors
computed by SciFort's public real Schur routine.  The complex triangular
factor is *not* compared elementwise (phases and eigenvalue order vary).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.optimize import linear_sum_assignment
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_matrix_functions.f90",
           "src/linalg_schur.f90"]
ORDERS = (0, 1, 2, 3, 4, 6, 8, 12, 20, 32, 40, 64, 96, 128)


def build(work: Path, mode: str, main: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / (Path(main).stem + "_" + mode)
    cmd = ["gfortran", "-std=f2018", *flags,
           *(str(ROOT / source) for source in SOURCES), str(ROOT / main), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=work, text=True, capture_output=True, timeout=120)
    if p.returncode:
        raise RuntimeError(f"{mode} compilation failed:\n{p.stderr[-4000:]}")
    return exe


def cases(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = tuple(n for n in ORDERS if n <= max_order)
    for case in range(count):
        n = sizes[case % len(sizes)]
        family = (case // len(sizes)) % 12
        a = rng.standard_normal((n, n))
        if family == 1:
            a = a + a.T
        elif family == 2:
            a = a - a.T
        elif family == 3:
            a = np.triu(a)
        elif family == 4:
            a = np.tril(a)
        elif family == 5:
            a = np.eye(n) * 2.0
            if n > 1:
                a += np.triu(rng.standard_normal((n, n)), 1) * 0.03
        elif family == 6:
            a = np.triu(a, 1)
        elif family == 7:
            a *= 1e130
        elif family == 8:
            a *= 1e-130
        elif family == 9:
            a *= 1e-260
        elif family == 10:
            # A chosen block with an explicit conjugate eigenvalue pair.
            if n > 1:
                a[:2, :2] = np.array([[0., -5.], [5., 0.]])
        elif family == 11:
            # Nearly real conjugate eigenvalues, with an isolated subspace.
            a = np.diag(np.linspace(0.1, 0.9, n))
            if n > 1:
                a[:2, :2] = np.array([[0.3, -1e-6], [1e-6, 0.3]])
        # The real Schur engine is bounded to order 40 and may be less robust
        # with strongly nonnormal/underflow-scale cases; those are still
        # tested extensively via SciPy-generated real Schur factors.
        source = ("scifort" if case % 4 == 1 and n <= 40 and
                  family in (0, 1, 2, 3, 4, 5, 10) else "scipy")
        with threadpool_limits(limits=1):
            t, z = linalg.schur(a, output="real")
        yield a, t, z, source


def submit(exe: Path, tests, timeout: int):
    lines = [str(len(tests))]
    for a, t, z, source in tests:
        n = len(a)
        lines.append(f"{n} {1 if source == 'scifort' else 0}")
        if source == 'scifort':
            lines.extend(f"{v:.17e}" for v in a.ravel(order="F"))
        else:
            lines.extend(f"{v:.17e}" for v in t.ravel(order="F"))
            lines.extend(f"{v:.17e}" for v in z.ravel(order="F"))
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
    p = subprocess.run([str(exe)], input="\n".join(lines) + "\n", text=True,
                       capture_output=True, timeout=timeout, env=env)
    if p.returncode:
        raise RuntimeError(f"Fortran driver failed:\n{p.stdout[-1000:]}\n{p.stderr[-4000:]}")
    return iter(p.stdout.splitlines())


def check(exe: Path, tests, tolerance: float, timeout: int):
    output = submit(exe, tests, timeout)
    worst = {"reconstruction": 0., "unitarity": 0., "triangularity": 0., "eigenvalues": 0.}
    count_internal = 0
    for case, (a, t, z, source) in enumerate(tests):
        n = len(a)
        try:
            status = int(next(output))
        except (StopIteration, ValueError) as exc:
            raise RuntimeError(f"case {case}: truncated or malformed driver output") from exc
        if status:
            raise AssertionError(f"case {case}: n={n}, source={source}, status={status}")
        if n == 0:
            continue
        raw = np.fromiter((float(v) for _ in range(n*n) for v in next(output).split()),
                          dtype=float, count=4*n*n).reshape(n*n, 4)
        tc = (raw[:, 0]+1j*raw[:, 1]).reshape((n, n), order="F")
        zc = (raw[:, 2]+1j*raw[:, 3]).reshape((n, n), order="F")
        scale = max(float(np.max(np.abs(a))), 1e-300)
        norm_a = max(float(linalg.norm(a/scale)), 1.0)
        reconstruction = float(linalg.norm(zc@(tc/scale)@zc.conj().T - a/scale))/norm_a
        unitarity = float(linalg.norm(zc.conj().T@zc-np.eye(n)))/max(1, n)
        triangularity = float(linalg.norm(np.tril(tc/scale, -1)))/norm_a
        # Normalize *before* SciPy's conversion.  At ~1e-260, SciPy's
        # rsf2csf/eigvals on unscaled data may lose imaginary parts due to
        # its underlying tiny-matrix eigenvalue solver; this is a known
        # reference limitation, not a SciFort conversion failure.
        with threadpool_limits(limits=1):
            ref_t, _ = linalg.rsf2csf(t / scale, z)
        ref = np.diag(ref_t)
        eig = np.diag(tc)
        cost = np.abs(eig[:, None]/scale - ref[None, :])
        rr, cc = linear_sum_assignment(cost)
        eigenvalues = float(np.max(cost[rr, cc]))/norm_a
        results = {"reconstruction": reconstruction, "unitarity": unitarity,
                   "triangularity": triangularity, "eigenvalues": eigenvalues}
        for key, val in results.items():
            worst[key] = max(worst[key], val)
            if not np.isfinite(val) or val > tolerance:
                raise AssertionError(f"case {case}, n={n}, source={source}, "
                                     f"{key}={val:.3e} > {tolerance:.3e}")
        count_internal += (source == "scifort")
    return worst, count_internal


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    parser.add_argument("--cases", type=int, default=240)
    parser.add_argument("--max-order", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--tolerance", type=float, default=3e-8)
    parser.add_argument("--timeout", type=int, default=150)
    args = parser.parse_args()
    if not 2 <= args.max_order <= 128 or args.cases < 1:
        parser.error("max-order must be 2..128 and cases must be positive")
    test_cases = list(cases(args.cases, args.max_order, args.seed))
    with tempfile.TemporaryDirectory(prefix="scifort_rsf2csf_") as directory:
        work = Path(directory)
        for mode in (["checked", "optimized"] if args.mode == "both" else [args.mode]):
            native = build(work, mode, "test/test_linalg_rsf2csf.f90")
            p = subprocess.run([str(native)], cwd=work, capture_output=True, text=True, timeout=40)
            if p.returncode:
                raise RuntimeError(f"{mode} native tests failed:\n{p.stdout}\n{p.stderr}")
            exe = build(work, mode, "tools/linalg_reference/driver_rsf2csf.f90")
            worst, internal = check(exe, test_cases, args.tolerance, args.timeout)
            print(f"{mode}: native tests passed; {len(test_cases)} SciPy comparisons passed "
                  f"({internal} via SciFort real Schur); " +
                  ", ".join(f"{k}={v:.3e}" for k,v in worst.items()), flush=True)


if __name__ == "__main__":
    main()
