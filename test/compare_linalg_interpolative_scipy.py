#!/usr/bin/env python3
"""Independent deterministic SciPy ID comparisons (no vendored SciPy code).

Run: python test/compare_linalg_interpolative_scipy.py ./driver_interpolative
The driver is compiled from validation/reference_drivers/driver_linalg_interpolative.f90 alongside
kinds, linalg, linalg_svd, linalg_qr, and linalg_interpolative.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
import scipy
from scipy.linalg import interpolative as sli


def build_cases():
    rng = np.random.default_rng(20261109)
    cases = []
    for dtype in (np.float64, np.complex128):
        for m in range(1, 13):
            for n in range(1, 15):
                for variant in range(4):
                    a = rng.standard_normal((m, n))
                    if dtype is np.complex128:
                        a = a + 1j * rng.standard_normal((m, n))
                    if variant == 1:  # rank-one (important for tolerance mode)
                        a = np.outer(a[:, 0], a[0, :])
                    if variant == 2:  # controllable low rank
                        r = max(1, min(m, n) // 2)
                        left = rng.standard_normal((m, r))
                        right = rng.standard_normal((r, n))
                        if dtype is np.complex128:
                            left = left + 1j * rng.standard_normal((m, r))
                            right = right + 1j * rng.standard_normal((r, n))
                        a = left @ right
                    if variant == 3:  # rescaled columns, with well-separated pivots
                        a = a * 10.0 ** np.linspace(-3, 3, n)
                    a = a.astype(dtype)
                    if variant == 3 and (m + n) % 2:
                        a = a * 1e-65
                    if variant == 2:
                        fixed_rank = max(1, min(m, n) // 2)
                    elif variant == 1:
                        fixed_rank = 1
                    else:
                        fixed_rank = max(1, min(m, n) - (m + n) % min(m, n))
                    cases.append((a, 0, fixed_rank, 0.0, variant))
                    if (m + n + variant) % 2 == 0:
                        cases.append((a, 1, 0, 1e-8, variant))
    return cases


def run(driver):
    cases = build_cases()
    parts = [str(len(cases))]
    for a, mode, rank, eps, variant in cases:
        m, n = a.shape
        parts.append(f"{int(np.iscomplexobj(a))} {m} {n} {mode} {rank} {eps:.18e}")
        for val in a.flat:
            parts.append(f"{float(np.real(val)):.18e} {float(np.imag(val)):.18e}")
    completed = subprocess.run(
        [driver], input="\n".join(parts) + "\n", text=True, capture_output=True,
        check=True, timeout=180,
    )
    lines = iter(completed.stdout.splitlines())
    worst_approx = 0.0
    worst_scipy_ratio = 0.0
    worst_exact_projection = 0.0
    same_pivots = 0
    pivot_comparable = 0
    fixed_rank_tests = 0
    tol_tests = 0
    comparisons = 0
    svd_comparisons = 0
    worst_svd_difference = 0.0
    for number, (a, mode, requested_rank, eps, variant) in enumerate(cases, 1):
        header = next(lines).split()
        if len(header) != 4 or header[0] != "CASE" or int(header[1]) != number:
            raise AssertionError((number, header))
        status, k = map(int, header[2:])
        if status != 0:
            raise AssertionError((number, a.shape, variant, mode, status))
        idx_line = next(lines).split()
        idx = np.asarray([int(t) - 1 for t in idx_line], dtype=int)
        if len(idx) != a.shape[1] or not np.array_equal(np.sort(idx), np.arange(a.shape[1])):
            raise AssertionError((number, "invalid permutation", idx))
        if mode == 0 and k != requested_rank:
            raise AssertionError((number, "wrong rank", k, requested_rank))
        proj = np.empty((k, a.shape[1] - k), dtype=a.dtype)
        for j in range(a.shape[1] - k):
            for i in range(k):
                re, im = map(float, next(lines).split())
                proj[i, j] = complex(re, im) if np.iscomplexobj(a) else re
        b = sli.reconstruct_skel_matrix(a, k, idx)
        p = sli.reconstruct_interp_matrix(idx, proj)
        approx = sli.reconstruct_matrix_from_id(b, idx, proj)
        if not np.allclose(approx, b @ p, rtol=1e-13, atol=1e-75):
            raise AssertionError((number, "reconstruction helpers disagree"))
        norm_a = np.linalg.norm(a)
        err = np.linalg.norm(a - approx) / norm_a
        worst_approx = max(worst_approx, err)
        if mode == 0:
            fixed_rank_tests += 1
            expected_idx, expected_proj = sli.interp_decomp(a, k, rand=False)
            expected_a = sli.reconstruct_matrix_from_id(
                a[:, expected_idx[:k]], expected_idx, expected_proj
            )
            expected_err = np.linalg.norm(a - expected_a) / norm_a
            if err > max(25.0 * expected_err, 1.0e-10):
                raise AssertionError((number, "approximation worse than SciPy ID", err, expected_err))
            worst_scipy_ratio = max(worst_scipy_ratio, err / max(expected_err, 1e-14))
            pivot_comparable += 1
            if np.array_equal(idx, expected_idx):
                same_pivots += 1
                if proj.size:
                    diff = np.linalg.norm(proj - expected_proj)
                    worst_exact_projection = max(
                        worst_exact_projection,
                        diff / max(1, np.linalg.norm(expected_proj)),
                    )
            else:
                # Different valid pivot orders need not have matching T.
                if not np.all(np.isfinite(proj)):
                    raise AssertionError((number, "nonfinite projection"))
        else:
            tol_tests += 1
            if err > max(3.0 * eps, 1.0e-12):
                raise AssertionError((number, "tolerance miss", err, eps, k))
            # Compare the SciPy deterministic rank (not necessarily identical:
            # SciFort bounds Frobenius error, SciPy primarily estimates 2-norm).
            k_scipy, _, _ = sli.interp_decomp(a, eps, rand=False)
            if not (0 <= k_scipy <= min(a.shape)):
                raise AssertionError((number, "bad SciPy rank", k_scipy))
        if k > 0 and number % 7 == 0:
            svd_header = next(lines).split()
            if svd_header != ["SVD", "0"]:
                raise AssertionError((number, "id_to_svd failed", svd_header))
            singular_values = np.asarray([float(t) for t in next(lines).split()])
            if singular_values.shape != (k,):
                raise AssertionError((number, "id_to_svd singular value count"))
            reference_u, reference_s, reference_v = sli.id_to_svd(b,idx,proj)
            scale_s = max(np.linalg.norm(reference_s), 1e-90)
            diff = np.linalg.norm(reference_s - singular_values) / scale_s
            worst_svd_difference = max(worst_svd_difference, diff)
            if diff > 2.0e-8:
                raise AssertionError((number, "id_to_svd singular value mismatch", diff))
            svd_comparisons += 1
        comparisons += 1
    if next(lines, None) is not None:
        raise AssertionError("Unexpected trailing output from Fortran")
    print(f"SciPy {scipy.__version__}; {comparisons} ID reference cases passed")
    print(f"Fixed rank: {fixed_rank_tests}; tolerance: {tol_tests}")
    print(f"Matching pivot orders: {same_pivots}/{pivot_comparable}")
    print(f"Max relative reconstruction error (including truncated rank): {worst_approx:.3e}")
    print(f"Max approximation-error ratio to SciPy: {worst_scipy_ratio:.3e}")
    print(f"Max relative projection difference on matching pivots: {worst_exact_projection:.3e}")
    print(f"ID-to-SVD singular-value comparisons: {svd_comparisons}")
    print(f"Max relative ID-to-SVD singular-value difference: {worst_svd_difference:.3e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: compare_linalg_interpolative_scipy.py FORTRAN_DRIVER")
    run(sys.argv[1])
