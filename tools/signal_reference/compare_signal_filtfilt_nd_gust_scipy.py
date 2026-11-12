#!/usr/bin/env python3
"""Check rank-1/2/5/15 Gustafsson filtfilt against SciPy 1.17.

Usage: python compare_signal_filtfilt_nd_gust_scipy.py DRIVER

The Fortran driver emits samples in column-major order. Complex-coefficient
cases use an independent complex least-squares reference rather than SciPy's
real-only observability matrix for complex input filters in SciPy 1.17.
"""
from __future__ import annotations

import sys
import subprocess

import numpy as np
from scipy import signal, linalg

LENGTHS = [1, 2, 5, 17, 47]
NCASES = 6 * 8 * 4 * 5


def coefficients(choice):
    if choice == 0:
        b, a = np.array([0.37]), np.array([1.2])
    elif choice == 1:
        b = np.array([.31, .12, .03])
        a = np.array([1., -.23, .07])
    elif choice == 2:
        b = np.array([.25, -.11])
        a = np.array([1.1, -.29, .07, -.014])
    else:
        b = np.array([.17, .07, .05, -.008])
        a = np.array([.9, -.21])
    cb, ca = b.astype(complex), a.astype(complex)
    if len(cb) > 1:
        cb[1] += .032j
    if len(ca) > 1:
        ca[1] -= .019j
    else:
        ca[0] += .07j
    return b, a, cb, ca


def complex_reference(b, a, x, irlen):
    """Solve Gustafsson's complex initial-state system in double precision."""
    n = len(x)
    order = max(len(b), len(a)) - 1
    if not order:
        return (b[0] / a[0]) ** 2 * x
    m = min(n, irlen) if n > 2 * irlen else n
    obs = np.empty((m, order), complex)
    for k in range(order):
        zi = np.zeros(order, complex)
        zi[k] = 1
        obs[:, k] = signal.lfilter(b, a, np.zeros(m, complex), zi=zi)[0]
    s = signal.lfilter(b, a, obs[::-1, :], axis=0)
    sr, obsr = s[::-1, :], obs[::-1, :]
    y0 = signal.lfilter(b, a, signal.lfilter(b, a, x)[::-1])[::-1]
    y1 = signal.lfilter(b, a, signal.lfilter(b, a, x[::-1])[::-1])
    if m == n:
        matrix = np.column_stack((sr - obs, obsr - s))
        correction = np.column_stack((sr, obsr))
        delta = y1 - y0
    else:
        matrix = np.zeros((2 * m, 2 * order), complex)
        correction = np.zeros_like(matrix)
        matrix[:m, :order] = sr - obs
        matrix[m:, order:] = obsr - s
        correction[:m, :order] = sr
        correction[m:, order:] = obsr
        delta = np.r_[y1[:m] - y0[:m], y1[-m:] - y0[-m:]]
    ic = linalg.lstsq(matrix, delta)[0]
    shift = correction @ ic
    result = y0.copy()
    if m == n:
        result += shift
    else:
        result[:m] += shift[:m]
        result[-m:] += shift[m:]
    return result


def expected(cid):
    rk = cid % 6
    typ = (cid // 6) % 8
    family = (cid // 48) % 4
    ircase = (cid // 192) % 5
    n = LENGTHS[(cid // 6 + cid // 48 + 7 * (cid // 192)) % 5]
    irlen = [10000, 1, 2, 7, 28][ircase]
    shapes = [(n,), (2, n), (n, 2), (2, 2, n, 2, 1),
              (n, 1, 2, *([1] * 11), 2), (2, 1, 2, *([1] * 11), n)]
    axes = [0, 1, 0, 2, 0, 14]
    shape, axis = shapes[rk], axes[rk]
    assert len(shape) in (1, 2, 5, 15)
    b, a, cb, ca = coefficients(family)
    coeffs = {0: (b, a, False), 1: (cb, ca, True), 2: (b, a, True),
              3: (b, ca, False), 4: (b, ca, True), 5: (cb, a, False),
              6: (cb, a, True), 7: (cb, ca, False)}
    b, a, complex_x = coeffs[typ]
    idx = np.arange(1, np.prod(shape) + 1)
    xr = (np.sin(.087 * idx + .03 * cid) + .003 * idx).reshape(shape, order='F')
    imag = .1 * np.cos(.12 * idx + .013 * cid).reshape(shape, order='F')
    x = xr + 1j * imag if complex_x else xr
    if np.iscomplexobj(b) or np.iscomplexobj(a):
        output = np.apply_along_axis(
            lambda line: complex_reference(b, a, line, irlen), axis, x)
    else:
        output = signal.filtfilt(b, a, x, axis=axis, method='gust', irlen=irlen)
    return output.ravel(order='F'), not (np.iscomplexobj(b) or np.iscomplexobj(a))


def main():
    driver = sys.argv[1]
    requested = [str(k) for k in range(NCASES)]
    completed = subprocess.run([driver], input='\n'.join(requested) + '\n',
                               capture_output=True, text=True, check=True)
    lines = iter(completed.stdout.splitlines())
    max_error = 0.
    max_scaled = 0.
    compared = 0
    direct_scipy = 0
    for cid in range(NCASES):
        header = next(lines).split()
        if len(header) != 3 or header[0] != 'CASE':
            raise AssertionError((cid, header))
        n, code = map(int, header[1:])
        assert code == 0, (cid, code)
        got = np.array([complex(*map(float, next(lines).split()))
                        for _ in range(n)])
        ref, direct = expected(cid)
        assert ref.size == n, (cid, n, ref.size)
        diff = np.max(np.abs(got - ref))
        scaled = diff / (1. + np.max(np.abs(ref)))
        if scaled >= 2e-10:
            raise AssertionError(f'case {cid}: error {diff:.3g}, scaled {scaled:.3g}')
        max_error = max(max_error, diff)
        max_scaled = max(max_scaled, scaled)
        compared += n
        direct_scipy += int(direct)
    assert next(lines, None) is None, 'extra driver output'
    print(f'PASS: {NCASES} rank-1/2/5/15 Gustafsson cases; '
          f'{compared} samples; {direct_scipy} direct SciPy comparisons; '
          f'{NCASES-direct_scipy} independent complex-coefficient references; '
          f'max abs error {max_error:.3g}; max scaled error {max_scaled:.3g}')


if __name__ == '__main__':
    main()
