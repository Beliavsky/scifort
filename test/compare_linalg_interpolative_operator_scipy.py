#!/usr/bin/env python3
"""Independent SciPy LinearOperator comparison for matrix-free Fortran ID.

Usage: python test/compare_linalg_interpolative_operator_scipy.py PATH_TO_DRIVER

Checks genuine forward/adjoint operator usage, selected original columns,
low-rank reconstruction, and adaptive residuals. Random pivots need not
match SciPy's independent random generator or randomized ID algorithm.
"""

from __future__ import annotations

import sys
import subprocess
import warnings
from pathlib import Path

import numpy as np
from scipy.linalg import interpolative
from scipy.sparse.linalg import aslinearoperator


def run_fortran(exe: Path, a: np.ndarray, target: int | float, seed: int) -> tuple:
    m, n = a.shape
    typ = int(np.iscomplexobj(a))
    is_tolerance = not isinstance(target, int)
    data = [f"{m} {n} {typ} {int(is_tolerance)} {target} {seed}\n"]
    for i in range(m):
        if typ:
            data.append(" ".join(f"({z.real:.17e},{z.imag:.17e})" for z in a[i, :]) + "\n")
        else:
            data.append(" ".join(f"{x:.17e}" for x in a[i, :]) + "\n")
    result = subprocess.run(
        [str(exe)], input="".join(data), text=True, capture_output=True,
        timeout=30, check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Fortran error: {result.stderr}\n{result.stdout[:2000]}")
    lines = result.stdout.splitlines()
    head = lines[0].split()
    status, k, est = int(head[1]), int(head[2]), float(head[3])
    if status != 0:
        raise RuntimeError(f"Fortran status {status}: {result.stdout[:2000]}")
    assert lines[2] == 'INDICES'
    assert lines[4] == 'PROJ'
    idx = np.fromstring(lines[3], sep=" ", dtype=np.int64) - 1
    calls = list(map(int, lines[1].split()[1:]))
    assert sorted(idx) == list(range(n))
    p = np.zeros((k, n - k), dtype=np.complex128 if typ else np.float64)
    for row in lines[5:]:
        if not row.strip():
            continue
        cols = row.split()
        i, j = int(cols[0]) - 1, int(cols[1]) - 1
        p[i, j] = float(cols[2]) + 1j*float(cols[3]) if typ else float(cols[2])
    assert len(lines[5:]) == k*(n-k)
    assert k <= min(m, n)
    assert calls[0] >= k and calls[1] >= min(k, n - k)
    b = a[:, idx[:k]]
    approx = np.zeros_like(a)
    approx[:, idx[:k]] = b
    approx[:, idx[k:]] = b @ p
    scale = float(np.max(abs(a))) if a.size else 0.0
    error = (
        float(np.linalg.norm((a - approx)/scale) / np.linalg.norm(a/scale))
        if scale else 0.0
    )
    return k, error, est, idx, p, calls


def test(exe: Path):
    rng = np.random.default_rng(81382)
    count = 0
    scipy_count = 0
    scipy_failures = 0
    maximum_adaptive_error = 0.0
    maximum_ratio_to_optimal = 0.0
    maximum_scipy_ratio = 0.0
    max_rank_gap = 0
    total_forward = 0
    total_adjoint = 0
    for m, n in [(3,5),(5,3),(6,8),(8,6),(7,7),(12,19),(19,12),(18,26),(26,18),
                 (23,37),(37,23),(36,45)]:
        for is_complex in (False, True):
            for variant in range(4):
                dtype = np.complex128 if is_complex else np.float64
                x = rng.normal(size=(m,n))
                if is_complex:
                    x = x + 1j * rng.normal(size=(m,n))
                if variant == 0:
                    a = x
                elif variant == 1:
                    rank = min(m,n,3)
                    u = rng.normal(size=(m,rank))
                    v = rng.normal(size=(rank,n))
                    if is_complex:
                        u = u + 1j * rng.normal(size=u.shape)
                        v = v + 1j * rng.normal(size=v.shape)
                    a = u @ v
                elif variant == 2:
                    a = x * 1.0e-80
                else:
                    a = x * 1.0e80
                assert a.dtype == dtype
                targets: list[int | float] = [1, min(m,n,3), min(m,n)//2, 0.25, 1.0e-3, 1.0e-7]
                for target in targets:
                    k, err, estimated, idx, proj, calls = run_fortran(exe,a,target,seed=22033)
                    total_forward += calls[0]
                    total_adjoint += calls[1]
                    assert np.isfinite(err)
                    if isinstance(target,int):
                        assert k == target
                        s = np.linalg.svd(a/np.max(abs(a)),compute_uv=False)
                        optimal = np.linalg.norm(s[k:]) / np.linalg.norm(s)
                        ratio = err / max(1.0e-12,optimal)
                        maximum_ratio_to_optimal = max(maximum_ratio_to_optimal,ratio)
                        if err > max(1.0e-9,3.0*optimal):
                            raise AssertionError(('fixed reconstruction',m,n,is_complex,variant,k,err,optimal,idx))
                        if target == 1 and variant in (0,1) and n > 10:
                            assert calls[0] == 1, ('unexpected matrix materialization',calls)
                    else:
                        maximum_adaptive_error = max(maximum_adaptive_error,err)
                        if err > max(2.0e-9,1.05*target):
                            raise AssertionError(('adaptive tolerance',m,n,is_complex,variant,k,err,target,estimated))
                        assert estimated <= max(1e-10,0.85*target), (estimated,target)
                    count += 1
                    # Independently exercise SciPy's actual LinearOperator pathway.
                    # Different randomized pivots, ranks and coefficients are valid.
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore')
                            result = interpolative.interp_decomp(
                                aslinearoperator(a),target,rng=np.random.default_rng(3312)
                            )
                        if isinstance(target,int):
                            sidx, sp = result
                            sk = target
                        else:
                            sk, sidx, sp = result
                        sb = a[:,sidx[:sk]]
                        ref_a = np.zeros_like(a)
                        ref_a[:,sidx[:sk]] = sb
                        ref_a[:,sidx[sk:]] = sb @ sp
                        scerr = np.linalg.norm((ref_a-a)/np.max(abs(a))) / np.linalg.norm(a/np.max(abs(a)))
                        assert np.isfinite(scerr)
                        scipy_count += 1
                        maximum_scipy_ratio = max(maximum_scipy_ratio,err/max(1e-12,scerr))
                        max_rank_gap = max(max_rank_gap,abs(k-sk))
                    except (ValueError, IndexError, np.linalg.LinAlgError, AssertionError, TypeError, RuntimeError):
                        scipy_failures += 1
    print(f"PASS: {count} SciFort matrix-free ID comparisons")
    print(f"SciPy LinearOperator calls: {scipy_count} succeeded; {scipy_failures} raised or were excluded")
    print(f"Maximum fixed-rank error/best-SVD-error ratio: {maximum_ratio_to_optimal:.5g}")
    print(f"Maximum relative adaptive reconstruction error: {maximum_adaptive_error:.5g}")
    print(f"Maximum Fortran error / SciPy error: {maximum_scipy_ratio:.5g}")
    print(f"Maximum difference in adaptive numerical ranks: {max_rank_gap}")
    print(f"Total Fortran callback calls: {total_forward} forward, {total_adjoint} adjoint")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: compare_linalg_interpolative_operator_scipy.py FORTRAN_DRIVER')
    test(Path(sys.argv[1]).resolve())
