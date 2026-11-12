#!/usr/bin/env python3
"""Check the genuine FITPACK smoothing adapter against independent SciPy tools.

Two fits can have distinct adaptive knot sets while satisfying the same
weighted residual target. We test the SciFort adapter's TCK representation
against scipy.interpolate.BSpline constructed from those SAME knots, then
compare independently fitted SciPy UnivariateSpline solutions separately.

Use --executable PATH for a local FITPACK build without FPM, or omit it to
build the optional dependency through FPM (which requires network access).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import numpy as np
from scipy.interpolate import BSpline, UnivariateSpline

ROOT = Path(__file__).resolve().parents[2]
OPTIONAL = ROOT / 'optional' / 'fitpack_backend'


def reference_executable(supplied: Path | None) -> Path:
    if supplied is not None:
        exe = supplied.resolve()
        if not exe.is_file():
            raise FileNotFoundError(exe)
        return exe
    subprocess.run(['fpm', 'build', '--compiler', 'gfortran'], cwd=OPTIONAL, check=True)
    candidates = [f for f in (OPTIONAL / 'build').rglob('fitpack_reference_driver*')
                  if f.is_file() and f.name in ('fitpack_reference_driver', 'fitpack_reference_driver.exe')]
    if len(candidates) != 1:
        raise RuntimeError(f'Expected one FPM executable, found: {candidates}')
    return candidates[0]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable', type=Path, help='Already linked real FITPACK reference driver')
    p.add_argument('--cases', type=int, default=140)
    p.add_argument('--seed', type=int, default=20261101)
    args = p.parse_args()
    exe = reference_executable(args.executable)
    rng = np.random.default_rng(args.seed)
    max_tck_error = 0.0
    max_residual_error = 0.0
    max_agreeing_knots_error = 0.0
    max_different_knots_error = 0.0
    different_knots = 0

    for i in range(args.cases):
        n = 25 + (i * 13) % 75
        k = 1 + i % 5
        x = np.cumsum(rng.uniform(.1, .8, size=n))
        x = (x - x[0]) / (x[-1] - x[0]) * 6.0
        y = np.sin(x * (1 + i % 4) * .5) + .15 * np.cos(3 * x) + rng.normal(scale=.035, size=n)
        w = np.exp(rng.uniform(-.5, .5, size=n))
        s = [0.0, .005, .05, .2, 1.0, 3.0][i % 6]
        q = np.linspace(x[0], x[-1], 71)
        lines = [f'{n} {k} {len(q)} {s:.17g}']
        lines.extend(' '.join(f'{v:.17g}' for v in arr) for arr in (x, y, w, q))
        result = subprocess.run([str(exe)], input='\n'.join(lines) + '\n', text=True,
                                capture_output=True, check=True)
        output = result.stdout.splitlines()
        if not output or int(output[0]) != 0:
            raise AssertionError(f'case {i}: FITPACK returned {output[:2]}')
        residual = float(output[1])
        nt, nc, degree = map(int, output[2].split())
        knots = np.fromstring(output[3], sep=' ')
        coeff = np.fromstring(output[4], sep=' ')
        values = np.fromstring(output[5], sep=' ')
        if nt != nc or nt != len(knots) or nc != len(coeff) or degree != k:
            raise AssertionError(f'case {i}: malformed TCK')

        # Independent check of the SciFort TCK evaluation and FITPACK residual.
        same_tck = BSpline(knots, coeff, k)
        tck_values = same_tck(q)
        tck_error = float(np.max(np.abs(values - tck_values)))
        max_tck_error = max(max_tck_error, tck_error)
        np.testing.assert_allclose(values, tck_values, rtol=1.e-10, atol=1.e-10,
                                   err_msg=f'case {i}: SciFort TCK evaluation')
        computed_residual = float(np.sum((w * (y - same_tck(x))) ** 2))
        residual_error = abs(computed_residual - residual)
        max_residual_error = max(max_residual_error, residual_error)
        np.testing.assert_allclose(residual, computed_residual, rtol=1.e-9, atol=1.e-9,
                                   err_msg=f'case {i}: FITPACK reported residual')
        if s == 0.0:
            np.testing.assert_allclose(same_tck(x), y, rtol=1.e-8, atol=1.e-8,
                                       err_msg=f'case {i}: exact interpolation')
        elif residual > s * 1.001 + 1.e-9:
            raise AssertionError(f'case {i}: residual {residual} exceeds s={s}')

        ref = UnivariateSpline(x, y, w=w, k=k, s=s)
        ref_t = ref._eval_args[0]
        has_same_knots = len(knots) == len(ref_t) and np.allclose(
            knots, ref_t, rtol=1.e-9, atol=1.e-9)
        discrepancy = float(np.max(np.abs(values - ref(q))))
        if has_same_knots:
            max_agreeing_knots_error = max(max_agreeing_knots_error, discrepancy)
            np.testing.assert_allclose(values, ref(q), rtol=1.e-7, atol=1.e-7,
                                       err_msg=f'case {i}: equal-knot SciPy evaluation')
        else:
            different_knots += 1
            max_different_knots_error = max(max_different_knots_error, discrepancy)

    print(f'PASS: {args.cases} genuine FITPACK/SciFort TCK and weighted-residual cases')
    print(f'max same-TCK evaluation error: {max_tck_error:.3g}')
    print(f'max reported-versus-computed residual error: {max_residual_error:.3g}')
    print(f'knot sets differing from SciPy: {different_knots}/{args.cases}')
    print(f'max SciPy evaluation difference for agreeing knots: {max_agreeing_knots_error:.3g}')
    print(f'max SciPy evaluation difference for different knots: {max_different_knots_error:.3g}')
    print('Different adaptive knots are reported, not treated as equal-fit validation.')


if __name__ == '__main__':
    main()
