#!/usr/bin/env python3
"""Compare SciFort continuous Gramians and balancing with independent SciPy.

Usage: python verify_clti_gramians.py /path/to/clti_gramians_driver
"""
from __future__ import annotations

import subprocess
import sys
import numpy as np
from scipy import linalg

FREQUENCIES = [0.0, 0.10, 0.30, 0.60, 1.0, 1.6, 2.5, 4.0]


def complex_random(rng, shape, scale=1.0):
    return scale * (rng.normal(size=shape) + 1j * rng.normal(size=shape)) / np.sqrt(2.0)


def make_case(rng, i):
    n = i % 6
    nu, ny = max(1, n), max(1, n)
    r = (i * 7 // 3) % (n + 1)
    if n:
        q, _ = linalg.qr(complex_random(rng, (n, n)))
        # Stable, mixed real/imaginary poles; add mild nonnormal similarity.
        rates = rng.uniform(0.40, 2.5, n)
        frequencies = rng.uniform(-0.9, 0.9, n)
        basis = q @ (np.eye(n) + complex_random(rng, (n, n), 0.05))
        a = basis @ np.diag(-rates + 1j * frequencies) @ linalg.inv(basis)
        b = complex_random(rng, (n, nu), 0.4) + np.eye(n, nu)
        c = complex_random(rng, (ny, n), 0.4) + np.eye(ny, n)
    else:
        a = np.zeros((0, 0), complex)
        b = np.zeros((0, nu), complex)
        c = np.zeros((ny, 0), complex)
    d = complex_random(rng, (ny, nu), 0.04) if i % 3 == 0 else np.zeros((ny, nu), complex)
    return (n, nu, ny, r), (a, b, c, d)


def read_matrix(lines, index, shape):
    count = int(np.prod(shape))
    vals = [complex(*map(float, lines[j].split())) for j in range(index, index + count)]
    return np.array(vals, dtype=complex).reshape(shape, order='F'), index + count


def gramian_reference(a, b, c):
    if a.shape[0] == 0:
        return np.zeros((0, 0), complex), np.zeros((0, 0), complex)
    wc = linalg.solve_continuous_lyapunov(a, -(b @ b.conj().T))
    wo = linalg.solve_continuous_lyapunov(a.conj().T, -(c.conj().T @ c))
    return wc, wo


def balance_reference(a, b, c, d, wc, wo, r):
    n = a.shape[0]
    if n == 0:
        return np.empty(0), (a, b, c, d)
    lc = linalg.cholesky(wc, lower=True)
    lo = linalg.cholesky(wo, lower=True)
    u, hsv, vh = linalg.svd(lo.conj().T @ lc)
    if r == 0:
        return hsv, (a[:0, :0], b[:0, :], c[:, :0], d)
    t = lc @ vh.conj().T[:, :r] / np.sqrt(hsv[:r])
    ti = (u.conj().T[:r, :] @ lo.conj().T) / np.sqrt(hsv[:r, None])
    return hsv, (ti @ a @ t, ti @ b, c @ t, d)


def response(matrices, w):
    a, b, c, d = matrices
    n = a.shape[0]
    if not n:
        return d
    return c @ linalg.solve(1j * w * np.eye(n) - a, b, assume_a='gen') + d


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    rng = np.random.default_rng(20261025)
    cases = [make_case(rng, i) for i in range(180)]
    tokens = [str(len(cases))]
    for (n, nu, ny, r), arrays in cases:
        tokens.append(f'{n} {nu} {ny} {r}')
        for mat in arrays:
            for x in mat.ravel(order='F'):
                tokens.append(f'({x.real:.17g},{x.imag:.17g})')
    result = subprocess.run([sys.argv[1]], input='\n'.join(tokens) + '\n',
                            capture_output=True, text=True, check=True)
    lines = result.stdout.splitlines()
    index = 0
    max_diff = 0.0
    values = 0
    for i, ((n, nu, ny, r), matrices) in enumerate(cases):
        header = list(map(int, lines[index].split()))
        index += 1
        if header != [i + 1, 0, 0]:
            raise AssertionError((i, 'statuses', header))
        a, b, c, d = matrices
        wc, index = read_matrix(lines, index, (n, n))
        wo, index = read_matrix(lines, index, (n, n))
        h2 = float(lines[index]); index += 1
        hsv = np.array([float(x) for x in lines[index:index+n]])
        index += n
        error_bound = float(lines[index]); index += 1
        actual = []
        for w in FREQUENCIES:
            mat, index = read_matrix(lines, index, (ny, nu))
            actual.append(mat)
        wc_ref, wo_ref = gramian_reference(a, b, c)
        hsv_ref, matrices_ref = balance_reference(a, b, c, d, wc_ref, wo_ref, r)
        max_err = max(np.max(np.abs(wc-wc_ref), initial=0),
                      np.max(np.abs(wo-wo_ref), initial=0),
                      np.max(np.abs(hsv-hsv_ref), initial=0),
                      abs(error_bound - 2 * sum(hsv_ref[r:])))
        if max_err > 3e-8:
            raise AssertionError((i, 'Gramians/Hankel values', max_err))
        if np.count_nonzero(d):
            if not np.isinf(h2):
                raise AssertionError((i, 'H2 feedthrough', h2))
        else:
            reference_h2 = np.sqrt(max(0.0, np.trace(c @ wc_ref @ c.conj().T).real))
            max_err = max(max_err, abs(h2-reference_h2))
            if abs(h2-reference_h2) > 3e-8:
                raise AssertionError((i, 'H2 strictly proper', h2, reference_h2))
        for k, w in enumerate(FREQUENCIES):
            expected = response(matrices_ref, w)
            delta = np.max(abs(actual[k]-expected))
            max_err = max(max_err, delta)
            if delta > 2e-7:
                raise AssertionError((i, 'balanced frequency', w, delta))
            original = response(matrices, w)
            spectral_error = linalg.svdvals(original-actual[k])[0]
            if spectral_error > error_bound + 2e-9:
                raise AssertionError((i, 'H-infinity bound', spectral_error, error_bound))
        max_diff = max(max_diff, max_err)
        values += 2*n*n+n+1+len(FREQUENCIES)*ny*nu
    if index != len(lines):
        raise AssertionError(('unread output', len(lines)-index))
    print(f'PASS: {len(cases)} continuous MIMO Gramian/reduction cases, '
          f'{values} values, maximum discrepancy {max_diff:.3e}')


if __name__ == '__main__':
    main()
