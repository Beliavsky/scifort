#!/usr/bin/env python3
"""Independent SciPy comparisons of real/complex dense triangular and LU solves."""
from __future__ import annotations

from dataclasses import dataclass
import subprocess
from pathlib import Path

import numpy as np
from scipy import linalg


@dataclass
class Case:
    variant: int  # 1 real triangular, 2 complex triangular, 3 real LU, 4 complex LU
    a: np.ndarray
    b: np.ndarray
    op: str
    lower: bool
    unit: bool
    reference: np.ndarray


def build_cases(count: int, seed: int, max_order: int, stress: bool = False) -> list[Case]:
    rng = np.random.default_rng(seed)
    cases: list[Case] = []
    sizes = [1, 2, 3, 4, 7, 11, max_order]
    if stress:
        sizes.extend([16, 32, 64, 120, 180, 256])
    for idx in range(count):
        variant = 1 + idx % 4
        complex_ = variant in (2, 4)
        triangular = variant <= 2
        lower = bool(rng.integers(2))
        unit = triangular and bool(rng.integers(2))
        op = 'NTC'[idx % 3]
        n = int(rng.choice(sizes if idx % 3 else [int(rng.integers(1, max_order + 1))]))
        n = max(n, 1)
        nrhs = int(rng.choice([1, 2, 3, 4]))
        a = rng.standard_normal((n, n)) / np.sqrt(n)
        if complex_:
            a = a + 1j * rng.standard_normal((n, n)) / np.sqrt(n)
        if triangular:
            # Unused half intentionally remains nonzero to test triangle selection.
            a[np.diag_indices(n)] = 0 if unit else rng.choice([-1.0, 1.0], size=n) * (
                1.5 + rng.random(n))
            if idx % 9 == 0:
                a = a * 0.001
                if unit:
                    a[np.diag_indices(n)] = 0
                else:
                    a[np.diag_indices(n)] = rng.choice([-1.0, 1.0], size=n)
            chosen = np.tril(a) if lower else np.triu(a)
            if unit:
                chosen = chosen.copy()
                chosen[np.diag_indices(n)] = 1
        else:
            a[np.diag_indices(n)] += 2.0
            if n >= 2 and idx % 5 == 0:
                # Swap the first row at elimination step one.
                a[0, 0] = 0
                a[n - 1, 0] = 5.0
            if idx % 11 == 0:
                a = a * 10.0
            chosen = a
        if op == 'T':
            chosen = chosen.T
        elif op == 'C':
            chosen = chosen.conj().T
        truth = rng.standard_normal((n, nrhs))
        if complex_:
            truth = truth + 1j * rng.standard_normal((n, nrhs))
        b = chosen @ truth
        if triangular:
            ref = linalg.solve_triangular(a, b, lower=lower, trans=op,
                                          unit_diagonal=unit, check_finite=True)
        else:
            factors, piv = linalg.lu_factor(a)
            ref = linalg.lu_solve((factors, piv), b, trans='NTC'.index(op))
        cases.append(Case(variant, a, b, op, lower, unit, np.asarray(ref)))
    return cases


def run(executable: Path, cases: list[Case], tolerance: float) -> dict[str, float]:
    lines = [str(len(cases))]
    for case in cases:
        n, nrhs = case.b.shape
        lines.append(f'{case.variant} {n} {nrhs} {case.op} {int(case.lower)} {int(case.unit)}')
        for array in (case.a, case.b):
            for v in np.asarray(array).ravel(order='F'):
                lines.append(f'{np.real(v):.17e} {np.imag(v):.17e}')
    result = subprocess.run([str(executable)], input='\n'.join(lines)+'\n', text=True,
                            capture_output=True, check=True, timeout=180)
    outputs = result.stdout.splitlines()
    pos = 0
    worst = 0.0
    worst_backward = 0.0
    for idx, case in enumerate(cases):
        n, nrhs = case.b.shape
        if pos == len(outputs):
            raise AssertionError(f'Fortran output ended early at case {idx}')
        status = int(outputs[pos].strip())
        pos += 1
        if status != 0:
            raise AssertionError(f'Fortran status {status} at case {idx}, '
                                 f'variant={case.variant}, n={n}, op={case.op}')
        if pos + n*nrhs > len(outputs):
            raise AssertionError(f'Fortran solution truncated at case {idx}')
        values = np.array([complex(*map(float, ln.split()))
                          for ln in outputs[pos:pos+n*nrhs]])
        pos += n*nrhs
        obtained = values.reshape((n, nrhs), order='F')
        if case.variant in (1, 3):
            if np.max(np.abs(obtained.imag)) != 0:
                raise AssertionError(f'Unexpected imaginary output at case {idx}')
            obtained = obtained.real
        error = np.linalg.norm(obtained - case.reference) / max(1, np.linalg.norm(case.reference))
        worst = max(worst, float(error))
        chosen = np.tril(case.a) if case.lower else np.triu(case.a)
        if case.variant >= 3:
            chosen = case.a
        elif case.unit:
            chosen = chosen.copy()
            chosen[np.diag_indices(n)] = 1
        if case.op == 'T':
            chosen = chosen.T
        elif case.op == 'C':
            chosen = chosen.conj().T
        residual = np.linalg.norm(chosen @ obtained - case.b)
        denom = np.linalg.norm(chosen)*np.linalg.norm(obtained) + np.linalg.norm(case.b)
        backward = float(residual / max(denom, 1e-300))
        worst_backward = max(worst_backward, backward)
        if not np.isfinite(error) or error > tolerance or backward > tolerance:
            raise AssertionError(f'SciPy mismatch case={idx} variant={case.variant} n={n} '
                                 f'error={error:g} backward={backward:g}')
    if pos != len(outputs):
        raise AssertionError(f'Unexpected trailing output lines: {len(outputs)-pos}')
    print(f'{len(cases)} SciPy comparisons passed; max relative solution delta '
          f'{worst:.3e}, max normwise backward residual {worst_backward:.3e}')
    return {'cases': len(cases), 'max_relative_error': worst, 'max_backward_residual': worst_backward}
