#!/usr/bin/env python3
"""Validate native 3D LinearND against Qhull-backed SciPy.

Builds only kinds, KD-tree, scattered interpolation and the streaming driver;
no FPM or third-party Fortran dependencies are needed. Run two independent seeds
for reproducible stochastic regression coverage.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import tempfile
import numpy as np
from scipy.interpolate import LinearNDInterpolator, griddata
from compare_scattered_scipy import build, run_case


def check(exe: Path, rng: np.random.Generator, case_id: int) -> tuple[int, float]:
    n = (4, 6, 9, 16, 30, 55, 85, 120)[case_id % 8]
    channels = 1 + (case_id % 4)
    is_complex = case_id % 3 == 0
    rescale = case_id % 2 == 0
    use_griddata = (case_id // 2) % 2 == 0
    x = rng.normal(size=(n, 3))
    if case_id % 7 == 0:
        x *= np.array([8.0, 0.6, 0.08])
    if case_id % 13 == 0:
        x += rng.uniform(-3000, 3000, size=3)
    if case_id % 17 == 0:
        x *= np.array([1.0, 0.01, 0.2])
    y = rng.normal(size=(n, channels))
    if is_complex:
        y = y + 1j * rng.normal(size=y.shape)
    nq = 24
    # Interpolate among randomly selected tetrahedron vertices to exercise
    # close-to-boundary and interior queries, plus exact sample-site queries.
    w = rng.dirichlet(np.ones(4), size=nq)
    indices = np.array([rng.choice(n, size=4, replace=False) for _ in range(nq)])
    qi = np.einsum("qi,qij->qj", w, x[indices])
    outside = x.mean(axis=0) + 10.0 * (x.max(axis=0) - x.min(axis=0))
    q = np.concatenate((x[:min(n, 5)], qi, outside[None, :]))
    fill = -17.0 + (2.0j if is_complex else 0.0)
    mode = 3 if use_griddata else 1
    actual = run_case(exe, mode, x, y, q, fill, rescale, is_complex)
    if use_griddata:
        ref = griddata(x, y, q, method='linear', fill_value=fill, rescale=rescale)
    else:
        ref = LinearNDInterpolator(x, y, fill_value=fill, rescale=rescale)(q)
    if not is_complex:
        actual = actual.real
    if ref.shape != actual.shape:
        raise AssertionError((case_id, ref.shape, actual.shape))
    error = float(np.max(np.abs(ref - actual)))
    scale = max(1.0, float(np.max(np.abs(ref))))
    if error > 2e-10 * scale:
        row, channel = np.unravel_index(np.argmax(np.abs(ref-actual)), ref.shape)
        raise AssertionError(f"case {case_id}: {error:.5g} at query {row}, channel {channel}; "
                             f"native={actual[row,channel]}, scipy={ref[row,channel]}")
    return int(ref.size), error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=20261031)
    parser.add_argument('--cases', type=int, default=100)
    parser.add_argument('--only', choices=['checked', 'optimized', 'both'], default='both')
    args = parser.parse_args()
    for name, flags in (('checked', ['-O0', '-g', '-fcheck=all', '-fbacktrace']),
                        ('optimized', ['-O3'])):
        if args.only != 'both' and args.only != name:
            continue
        rng = np.random.default_rng(args.seed)
        total, worst = 0, 0.0
        with tempfile.TemporaryDirectory() as tmp:
            exe = build(Path(tmp), flags)
            for case_id in range(args.cases):
                count, error = check(exe, rng, case_id)
                total += count
                worst = max(worst, error)
        print(f'{name}: {args.cases} 3D SciPy comparisons PASS, '
              f'{total} numerical values; max absolute error {worst:.3g}', flush=True)


if __name__ == '__main__':
    main()
