#!/usr/bin/env python3
"""Compare SciFort CSD with scipy.linalg.cossin on deterministic test matrices.

Usage: python tools/linalg_reference/compare_cossin_scipy.py <driver_cossin>
No external dependencies are needed by the compiled Fortran implementation.
"""
from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass

import numpy as np
from scipy.linalg import cossin

RNG = np.random.default_rng(20261109)


@dataclass
class Case:
    x: np.ndarray
    p: int
    q: int
    swap: bool
    invalid: bool = False


def make_cases() -> list[Case]:
    cases: list[Case] = []
    for complex_mode in (False, True):
        for m in range(2, 19):
            for iteration in range(45):
                p = int(RNG.integers(1, m))
                q = int(RNG.integers(1, m))
                raw = RNG.normal(size=(m, m))
                if complex_mode:
                    raw = raw + 1j * RNG.normal(size=(m, m))
                x, _ = np.linalg.qr(raw)
                if iteration % 12 == 0:
                    x = np.eye(m, dtype=complex if complex_mode else float)
                elif iteration % 12 == 1:
                    x = x[:, ::-1].copy()
                elif iteration % 12 == 2:
                    x = x[RNG.permutation(m), :].copy()
                elif iteration % 12 == 3:
                    # A block rotation with a tiny sine can be numerically
                    # more challenging than a Haar-distributed random matrix.
                    x = np.eye(m, dtype=complex if complex_mode else float)
                    t = 1e-12
                    x[0, 0] = x[1, 1] = np.cos(t)
                    x[0, 1] = -np.sin(t)
                    x[1, 0] = np.sin(t)
                for swap in (False, True):
                    cases.append(Case(x, p, q, swap))
        # A handful of expected shape and unitary-input failures.
        m = 6
        for p, q in ((0, 3), (6, 2), (2, 0), (2, 6)):
            cases.append(Case(np.eye(m,dtype=complex if complex_mode else float),p,q,False,True))
        x = np.eye(m,dtype=complex if complex_mode else float)
        x[1, 1] = 2
        cases.append(Case(x,2,2,False,True))
    return cases


def run(exe: str) -> None:
    cases = make_cases()
    lines = [str(len(cases))]
    for case in cases:
        x = case.x
        m = x.shape[0]
        kind = int(np.iscomplexobj(x))
        lines.append(f'{kind} {m} {case.p} {case.q} {int(case.swap)}')
        for value in x.ravel(order='F'):
            lines.append(f'{value.real:.17e} {value.imag:.17e}')
    proc = subprocess.run(
        [exe], input='\n'.join(lines)+'\n', capture_output=True, text=True, check=True
    )
    output = iter(proc.stdout.splitlines())
    max_matrix = 0.0
    max_angle = 0.0
    max_reconstruction = 0.0
    max_unitary = 0.0
    checked = 0
    for index, case in enumerate(cases):
        stat, err, ue, ve = next(output).split()
        status = int(stat)
        if case.invalid:
            if status == 0:
                raise AssertionError(f'case {index}: expected invalid status')
            checked += 1
            continue
        if status != 0:
            raise AssertionError(f'case {index}: Fortran status {status}')
        th = np.fromstring(next(output),sep=' ')
        mat = np.fromstring(next(output),sep=' ').reshape(case.x.shape,order='F')
        (_, _), angles, (_, _) = cossin(case.x, p=case.p, q=case.q,
                                       separate=True, swap_sign=case.swap)
        _, cs_ref, _ = cossin(case.x,p=case.p,q=case.q,swap_sign=case.swap)
        aerr = float(np.max(np.abs(angles-th)))
        merr = float(np.max(np.abs(mat-cs_ref)))
        rec = float(err)
        unitarity = max(float(ue),float(ve))
        if aerr > 5e-9 or merr > 5e-9 or rec > 5e-9 or unitarity > 5e-9:
            raise AssertionError(f'case {index} m={len(case.x)} '
                                 f'p={case.p} q={case.q} swap={case.swap}: '
                                 f'angles={aerr:.4e} cs={merr:.4e} '
                                 f'reconstruction={rec:.4e} unitary={unitarity:.4e}')
        max_matrix=max(max_matrix,merr)
        max_angle=max(max_angle,aerr)
        max_reconstruction=max(max_reconstruction,rec)
        max_unitary=max(max_unitary,unitarity)
        checked += 1
    if next(output,None) is not None:
        raise AssertionError('unexpected additional Fortran output')
    print(f'PASS {checked} SciPy CSD cases (including {sum(c.invalid for c in cases)} negative tests)')
    print(f'max |CS-SciPy CS| = {max_matrix:.3e}')
    print(f'max |theta-SciPy theta| = {max_angle:.3e}')
    print(f'max elementwise reconstruction error = {max_reconstruction:.3e}')
    print(f'max unitary factor defect = {max_unitary:.3e}')


if __name__ == '__main__':
    if len(sys.argv)!=2:
        raise SystemExit('Usage: compare_cossin_scipy.py <compiled-driver>')
    run(sys.argv[1])
