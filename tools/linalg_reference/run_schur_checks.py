#!/usr/bin/env python3
"""Validate SciFort complex Schur and sorted complex Schur against SciPy.

All comparisons are based on reconstruction, unitarity, eigenvalue multisets,
and selected-subspace dimensions; individual Schur matrices are not unique.
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

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_matrix_functions.f90",
           "src/linalg_schur.f90"]
RULES = ["none", "lhp", "rhp", "iuc", "ouc", "positive-imag"]


def matrices(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = [0, 1, 2, 3, 5, 8, 12, max_order]
    for k in range(count):
        n = sizes[k % len(sizes)]
        complex_input = (k // len(sizes)) % 2 != 0
        a = rng.standard_normal((n, n))
        if complex_input:
            a = a + 1j * rng.standard_normal((n, n))
        family = (k // (2 * len(sizes))) % 12
        if family == 1:
            a = np.triu(a)
        elif family == 2:
            a = np.tril(a)
        elif family == 3:
            a = np.diag(np.linspace(-2.5, 2.5, n)).astype(a.dtype)
        elif family == 4:
            a = a + a.conj().T
        elif family == 5:
            a = a - a.conj().T
        elif family == 6:
            a *= 1e110
        elif family == 7:
            a *= 1e-110
        elif family == 8:
            a = np.eye(n, dtype=a.dtype) * 2.3
        elif family == 9:
            a = np.triu(a, 1)  # defective/zero eigenvalues
        elif family == 10:
            a = np.diag(np.linspace(-3, 3, n).astype(a.dtype))
            if n > 1:
                a[0, -1] = 0.25
        elif family == 11:
            # Regular upper triangular matrix to exercise stable reordering.
            a = np.diag(np.linspace(-4.3, 4.3, n).astype(a.dtype))
            if n > 1:
                a += np.triu(rng.standard_normal((n, n)), 1) * 0.15
        # With eigenvalues on a sorting boundary, SciPy/Schur rounding can
        # affect region membership. Test these in unsorted mode instead.
        sort = RULES[(k // len(sizes)) % len(RULES)]
        if family in (3, 5, 9, 10) and sort != "none":
            sort = "none"
        if n and sort != "none":
            # Near selection boundaries, a tiny Schur rounding change can
            # change membership without changing the decomposition accuracy.
            eig = np.linalg.eigvals(a)
            if sort in ("lhp", "rhp") and np.any(
                    np.abs(eig.real) <= 5e-12 * max(float(np.max(np.abs(eig))), 1e-300)):
                sort = "none"
            elif sort in ("iuc", "ouc") and np.any(np.abs(np.abs(eig) - 1) < 5e-12):
                sort = "none"
            elif sort == "positive-imag" and np.any(
                    np.abs(eig.imag) <= 5e-12 * max(float(np.max(np.abs(eig))), 1e-300)):
                sort = "none"
        yield a, sort


def native(work: Path, mode: str) -> None:
    exe = compile_source(work, mode, "test/test_linalg_schur.f90", "native")
    p = subprocess.run([str(exe)], cwd=work, capture_output=True, text=True, timeout=60)
    if p.returncode:
        raise RuntimeError(f"{mode}: native execution failed:\n{p.stdout}\n{p.stderr}")
    print(f"{mode}: {p.stdout.strip()}")


def compile_source(work: Path, mode: str, driver: str, suffix: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / (f"schur_{suffix}_{mode}")
    cmd = ["gfortran", "-std=f2018", *flags, *(str(ROOT / src) for src in SOURCES),
           str(ROOT / driver), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=90)
    if p.returncode:
        raise RuntimeError(f"{mode}: compile failed:\n{p.stderr}")
    return exe


def selected(w: np.ndarray, rule: str) -> np.ndarray:
    if rule == "lhp":
        return w.real < 0.0
    if rule == "rhp":
        return w.real >= 0.0
    if rule == "iuc":
        return np.abs(w) <= 1.0
    if rule == "ouc":
        return np.abs(w) > 1.0
    if rule == "positive-imag":
        return w.imag > 0.0
    return np.zeros(w.shape, dtype=bool)


def compare(exe: Path, cases, tolerance: float):
    buf = [str(len(cases))]
    for a, rule in cases:
        n = a.shape[0]
        buf.append(f"{2 if np.iscomplexobj(a) else 1} {n} {RULES.index(rule)}")
        buf.extend(f"{v.real:.17e} {v.imag:.17e}" for v in a.ravel(order="F"))
    p = subprocess.run([str(exe)], input="\n".join(buf) + "\n", text=True,
                       capture_output=True, timeout=180,
                       env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"})
    if p.returncode:
        raise RuntimeError(f"driver crashed:\n{p.stderr[-3000:]}")
    lines = iter(p.stdout.splitlines())
    worst = dict(reconstruction=0., orthogonality=0., triangularity=0.,
                 eigenvalues=0., invariant_subspace=0.)
    sorted_count = 0
    for k, (a, rule) in enumerate(cases):
        n = a.shape[0]
        status, sdim = [int(x) for x in next(lines).split()]
        if status != 0:
            raise AssertionError(f"case {k} n={n} rule={rule} status={status}")
        if n == 0:
            assert sdim == 0
            continue
        flat = np.asarray([[float(v) for v in next(lines).split()] for _ in range(n * n)])
        t = (flat[:, 0] + 1j * flat[:, 1]).reshape((n, n), order="F")
        z = (flat[:, 2] + 1j * flat[:, 3]).reshape((n, n), order="F")
        scale = max(float(np.max(np.abs(a))), 1e-300)
        denom = max(float(linalg.norm(a / scale)), 1.0)
        recon = float(linalg.norm(z @ (t / scale) @ z.conj().T - a / scale)) / denom
        orth = float(linalg.norm(z.conj().T @ z - np.eye(n))) / n
        tri = float(linalg.norm(np.tril(t / scale, -1))) / denom
        ref_eig = linalg.eigvals(a)
        eigen = np.diag(t)
        errors = np.abs(eigen[:, None] / scale - ref_eig[None, :] / scale)
        rr, cc = linear_sum_assignment(errors)
        spectra = float(np.max(errors[rr, cc])) / denom
        if rule != "none":
            # A Schur diagonal must have the selected eigenvalues leading.
            matches = selected(eigen, rule)
            if int(matches.sum()) != sdim or not np.all(matches[:sdim]) or np.any(matches[sdim:]):
                raise AssertionError(f"case {k} n={n} invalid sorted diagonal for {rule}")
            # Compare selected dimension to the independent SciPy sorted Schur
            # when it can be robustly determined.
            reference_selector = ((lambda w: w.imag > 0.0) if rule == "positive-imag" else rule)
            ref_sorted_t, _, ref_sdim = linalg.schur(a, output="complex", sort=reference_selector)
            if sdim != ref_sdim:
                raise AssertionError(f"case {k} n={n} rule={rule} selected {sdim}, SciPy {ref_sdim}")
            sorted_count += 1
        subspace = 0.0
        if rule != "none" and sdim:
            subspace = float(linalg.norm((a / scale) @ z[:, :sdim] -
                         z[:, :sdim] @ (t[:sdim, :sdim] / scale))) / denom
        for name, val in dict(reconstruction=recon, orthogonality=orth,
                              triangularity=tri, eigenvalues=spectra,
                              invariant_subspace=subspace).items():
            worst[name] = max(worst[name], val)
            if not np.isfinite(val) or val > tolerance:
                raise AssertionError(f"case {k} n={n} rule={rule}: {name}={val:.3g} > {tolerance}")
    return worst, sorted_count


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["checked", "optimized", "both"], default="both")
    ap.add_argument("--cases", type=int, default=180)
    ap.add_argument("--max-order", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20261108)
    ap.add_argument("--tolerance", type=float, default=2e-8)
    args = ap.parse_args()
    cases = list(matrices(args.cases, args.max_order, args.seed))
    with tempfile.TemporaryDirectory(prefix="scifort_schur_") as directory:
        for mode in (["checked", "optimized"] if args.mode == "both" else [args.mode]):
            work = Path(directory)
            native(work, mode)
            exe = compile_source(work, mode, "tools/linalg_reference/driver_schur.f90", "driver")
            worst, sorted_count = compare(exe, cases, args.tolerance)
            print(f"{mode}: {len(cases)} SciPy comparisons passed ({sorted_count} sorted); " +
                  ", ".join(f"{key} {value:.3e}" for key, value in worst.items()))


if __name__ == "__main__":
    main()
