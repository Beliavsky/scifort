#!/usr/bin/env python3
"""Compare the SciFort matrix exponential with SciPy on independent cases.

Usage:
    python tools/linalg_reference/compare_expm_scipy.py PATH/TO/driver_expm

The driver must be compiled with scifort_kinds, scifort_linalg,
scifort_linalg_dense, and scifort_linalg_expm. No custom Python dependencies
except NumPy and SciPy are needed.
"""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm


@dataclass
class MatrixCase:
    name: str
    a: np.ndarray


def cases() -> list[MatrixCase]:
    rng = np.random.default_rng(20261110)
    output: list[MatrixCase] = []
    for n in (0, 1, 2, 3, 4, 7, 12, 20, 32):
        for complex_input in (False, True):
            dtype = np.complex128 if complex_input else np.float64
            for version in range(11):
                for scale in (1e-8, 0.01, 0.5, 1.0, 4.0, 12.0, 30.0, 60.0):
                    a = rng.normal(size=(n,n)) / max(1.0, np.sqrt(n))
                    if complex_input:
                        a = a + 1j*rng.normal(size=(n,n))/max(1.0,np.sqrt(n))
                    if version == 1:
                        a = (a+a.T.conj())/2
                    elif version == 2:
                        a = (a-a.T.conj())/2
                    elif version == 3:
                        a = np.triu(a)
                    elif version == 4:
                        a = np.tril(a)
                    elif version == 5:
                        a = np.diag(np.diag(a))
                    elif version == 6:
                        a = np.eye(n, dtype=dtype) + a/10
                    elif version == 7:
                        a = np.triu(a, k=1)
                    elif version == 8:
                        if n:
                            a[:, -1] = a[:, 0]  # dependent columns
                    elif version == 9:
                        a = a.T
                    elif version == 10:
                        a = a + 0.4*np.eye(n)
                    a = np.asarray(a*scale, dtype=dtype)
                    output.append(MatrixCase(f'n={n},complex={complex_input},v={version},s={scale}',a))
    for n in (48, 64):
        for complex_input in (False, True):
            a = rng.normal(size=(n,n)) / np.sqrt(n)
            if complex_input:
                a = a + 1j*rng.normal(size=(n,n))/np.sqrt(n)
            output.append(MatrixCase(f'medium dense n={n}, complex={complex_input}',a))
    # Structured, nonnormal and high-frequency cases.
    for t in (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0):
        output.append(MatrixCase(f'rotation {t}', np.array([[0., -t],[t,0.]])))
        output.append(MatrixCase(f'complex oscillation {t}',
                                 np.array([[0., -t],[t,0.]],dtype=complex)*(1+min(0.1,1.0/t)*1j)))
    for t in (0.1, 1.0, 10., 100., 1e4, 1e8):
        output.append(MatrixCase(f'Jordan {t}',
                                 np.array([[0., t, 0.],[0.,0.,t],[0.,0.,0.]])))
    for diagonal in (np.array([-700., 0., 300.]),
                     np.array([-10., 5., 100.]),
                     np.array([0., 1e-100, -1e-100])):
        output.append(MatrixCase(f'large diagonal {diagonal}',np.diag(diagonal)))
    return output


def driver_input(cases_: list[MatrixCase]) -> str:
    lines = [str(len(cases_))]
    for c in cases_:
        n = c.a.shape[0]
        lines.append(f'{2 if np.iscomplexobj(c.a) else 1} {n}')
        for val in c.a.ravel():
            lines.append(f'{val.real:.17e} {val.imag:.17e}')
    return '\n'.join(lines)+'\n'


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver', help='compiled Fortran executable')
    args = parser.parse_args()
    data = cases()
    process = subprocess.run([args.driver],input=driver_input(data),
                             text=True, capture_output=True, check=True)
    tokens = iter(process.stdout.split())
    max_error = 0.0
    max_label = ''
    checks = 0
    elements = 0
    with_scaling = 0
    for idx, case in enumerate(data):
        try:
            status, n, scalings = int(next(tokens)), int(next(tokens)), int(next(tokens))
        except StopIteration as exc:
            raise AssertionError(f'Missing Fortran output for case {idx}') from exc
        assert n == case.a.shape[0], (idx,n,case.a.shape)
        assert status == 0, (case.name, status)
        got = np.empty((n,n), dtype=np.complex128)
        for i in range(n):
            for j in range(n):
                got[i,j] = float(next(tokens)) + 1j*float(next(tokens))
        expected = expm(case.a)
        if expected.size:
            if not np.all(np.isfinite(got)):
                raise AssertionError(('nonfinite result',case.name))
            error = float(np.linalg.norm(got-expected, ord='fro')/
                          max(1.0,np.linalg.norm(expected, ord='fro')))
        else:
            error = 0.0
        if error > max_error:
            max_error = error
            max_label = case.name
        # Complex nonnormal matrices can amplify roundoff; this is deliberately
        # looser than the expected observed error for most cases.
        assert error < 1e-8, (case.name,error,max_error)
        with_scaling += int(scalings > 0)
        checks += 1
        elements += n*n
    assert next(tokens,None) is None, 'unexpected trailing output'
    print(f'PASS: {checks} independent SciPy expm cases; {elements} matrix entries')
    print(f'Cases requiring scaling and squaring: {with_scaling}')
    print(f'Worst normalized Frobenius difference: {max_error:.3e} ({max_label})')


if __name__ == '__main__':
    main()
