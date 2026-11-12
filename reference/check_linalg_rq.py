#!/usr/bin/env python3
"""Compare SciFort real/complex RQ with scipy.linalg.rq (deterministic).

Usage: python reference/check_linalg_rq.py /path/to/linalg_rq_driver

The driver reads one multi-case text stream from standard input and writes R/Q
factors. Factor sign/phase is convention-dependent, so reconstruction,
row-orthonormality, R trapezoidal structure, and gauge-independent absolute R
entries are checked rather than direct R/Q equality on every case.
"""
from __future__ import annotations

import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy import linalg


def make_cases():
    rng = np.random.default_rng(20261109)
    cases = []
    # Every small matrix shape, including both empty dimensions.
    for m in range(13):
        for n in range(13):
            for iscomplex in (False, True):
                a = rng.standard_normal((m, n))
                if iscomplex:
                    a = a + 1j * rng.standard_normal((m, n))
                for mode in ('full', 'economic', 'r'):
                    cases.append((a, mode))
    # Rectangular, square, rank-deficient, zero, and extreme-scale examples.
    for _ in range(220):
        m = int(rng.integers(1, 28))
        n = int(rng.integers(1, 28))
        iscomplex = bool(rng.integers(2))
        a = rng.standard_normal((m, n))
        if iscomplex:
            a = a + 1j * rng.standard_normal((m, n))
        variant = int(rng.integers(6))
        if variant == 1:
            a[:, 0] = 0
        elif variant == 2:
            a[-1, :] = a[0, :]
        elif variant == 3:
            a *= 1e-50
        elif variant == 4:
            a *= 1e50
        elif variant == 5:
            a[:] = 0
        for mode in ('full', 'economic', 'r'):
            cases.append((a, mode))
    return cases


def input_stream(cases):
    stream = io.StringIO()
    stream.write(f'{len(cases)}\n')
    for a, mode in cases:
        m, n = a.shape
        ctype = np.iscomplexobj(a)
        stream.write(f'{2 if ctype else 1} {mode} {m} {n}\n')
        for v in a.ravel(order='F'):
            if ctype:
                stream.write(f'{v.real:.17e} {v.imag:.17e}\n')
            else:
                stream.write(f'{float(v):.17e}\n')
    return stream.getvalue()


def read_matrix(lines, rows, cols, ctype):
    if not rows or not cols:
        return np.empty((rows, cols), dtype=np.complex128 if ctype else np.float64)
    vals = []
    for _ in range(rows * cols):
        chunks = next(lines).split()
        if ctype:
            assert len(chunks) == 2, chunks
            vals.append(float(chunks[0]) + 1j * float(chunks[1]))
        else:
            assert len(chunks) == 1, chunks
            vals.append(float(chunks[0]))
    return np.asarray(vals).reshape((rows, cols), order='F')


def verify(cases, output):
    lines = iter(output.splitlines())
    maximum_reconstruction = 0.0
    maximum_orthogonality = 0.0
    maximum_absolute_r_difference = 0.0
    counts = {'full': 0, 'economic': 0, 'r': 0}
    for case_id, (a, mode) in enumerate(cases, start=1):
        header = next(lines).split()
        assert header[0] == 'CASE' and int(header[1]) == case_id, header
        assert int(header[2]) == 0, (case_id, mode, header)
        m, n = a.shape
        k = min(m, n)
        rows = k if mode == 'economic' else n
        assert (int(header[3]), int(header[4]), int(header[5])) == (m, rows, 0 if mode == 'r' else rows)
        r = read_matrix(lines, m, rows, np.iscomplexobj(a))
        q = None
        if mode != 'r':
            q = read_matrix(lines, rows, n, np.iscomplexobj(a))
        reference = linalg.rq(a, mode='r' if mode == 'r' else mode)
        reference_r = reference if mode == 'r' else reference[0]
        scale = max(1e-90, float(np.max(np.abs(a))) if a.size else 0.0)
        scale_r = max(1e-90, float(np.max(np.abs(reference_r))) if reference_r.size else 0.0)
        assert r.shape == reference_r.shape
        # Verify prescribed RQ trapezoidal structure exactly (up to sign-zero).
        for i in range(m):
            for j in range(rows):
                if j < i + rows - m:
                    assert r[i, j] == 0, (case_id, 'R structure', i, j)
        if q is not None:
            if m and n:
                rec_err = np.max(np.abs(r @ q - a)) / scale
                assert rec_err < 7e-12, (case_id, mode, 'reconstruction', rec_err)
                maximum_reconstruction = max(maximum_reconstruction, float(rec_err))
            if rows:
                gram = q @ q.conj().T
                orth_err = float(np.max(np.abs(gram - np.eye(rows))))
                assert orth_err < 7e-12, (case_id, mode, 'unitary', orth_err)
                maximum_orthogonality = max(maximum_orthogonality, orth_err)
            if m == 0 or n == 0:
                if mode == 'full':
                    np.testing.assert_array_equal(q, reference[1])
        # Magnitudes of entries are unaffected by row-phase choices in R,
        # on full-rank matrices (for rank-deficient cases compare invariants).
        if min(m,n) and np.linalg.matrix_rank(a / scale) == min(m,n):
            diff = float(np.max(np.abs(np.abs(r) - np.abs(reference_r)))) / max(scale, scale_r)
            assert diff < 7e-12, (case_id, mode, 'SciPy abs(R)', diff)
            maximum_absolute_r_difference = max(maximum_absolute_r_difference, diff)
        counts[mode] += 1
    assert next(lines, None) is None, 'unexpected driver output'
    return counts, maximum_reconstruction, maximum_orthogonality, maximum_absolute_r_difference


def main():
    if len(sys.argv) != 2:
        raise SystemExit('usage: check_linalg_rq.py /path/to/linalg_rq_driver')
    cases = make_cases()
    proc = subprocess.run([str(Path(sys.argv[1]).resolve())], input=input_stream(cases),
                          text=True, capture_output=True, check=True)
    counts, rec, orth, r = verify(cases, proc.stdout)
    print(f'PASS scipy={scipy.__version__} cases={len(cases)} full={counts["full"]} '
          f'economic={counts["economic"]} r_only={counts["r"]}')
    print(f'max_scaled_reconstruction={rec:.4e} max_unitarity={orth:.4e} '
          f'max_scaled_abs_R_difference={r:.4e}')


if __name__ == '__main__':
    main()
