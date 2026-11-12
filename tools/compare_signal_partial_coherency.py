#!/usr/bin/env python3
"""Compare complex partial coherency and shrinkage with independent NumPy/SciPy.

Usage: python tools/compare_signal_partial_coherency.py /path/to/fortran-driver
SciPy is used only for independent multitaper reference generation.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np

from compare_signal_tf_partial_coherence import reference_spectrogram


def reference_partial(s: np.ndarray, ridge: float, alpha: float):
    """Invert trace-loaded, diagonal-shrunk S using independent NumPy LAPACK."""
    nc = len(s)
    diagonal = np.diag(np.diag(s))
    regularized = (1.0 - alpha) * s + alpha * diagonal
    regularized = regularized + np.eye(nc) * ridge * np.trace(s).real / nc
    p = np.linalg.inv(regularized)
    d = np.diag(p).real
    c = -p / np.sqrt(d[:, None]) / np.sqrt(d[None, :])
    np.fill_diagonal(c, 1.0)
    return np.minimum(1.0, np.abs(c)**2), c, p


def build_cases():
    rng = np.random.default_rng(264794)
    cases = []
    for nc in (2, 3, 4):
        for alpha in (0.0, 0.07, 0.45, 1.0):
            for ridge in (0.0, 0.025):
                for trial in range(6):
                    mode = trial % 2
                    nf, nt = 2 + trial % 2, 1 if mode == 0 else 2 + trial % 2
                    s = np.empty((nf, nt, nc, nc), dtype=np.complex128)
                    singular = np.zeros((nf, nt), dtype=bool)
                    for f in range(nf):
                        for t in range(nt):
                            z = rng.normal(size=(nc + 3, nc)) + 1j * rng.normal(size=(nc + 3, nc))
                            if trial == 0 and f == 0 and t == 0:
                                z[:, -1] = z[:, -2]
                                singular[f, t] = alpha == 0.0 and ridge == 0.0
                            scales = 10.0**rng.uniform(-0.75, 0.75, size=nc)
                            z *= scales[None, :]
                            s[f, t] = z.conj().T @ z
                    cases.append((mode, ridge, alpha, s, singular))
    for mode in (2, 3):
        for nc in (2, 3, 4):
            for nperseg in (24, 32):
                for weighting in ('equal', 'eigen'):
                    for trend in ('none', 'linear'):
                        for alpha in (0.0, 0.10, 0.65):
                            n = 2*nperseg + 7
                            x = rng.normal(size=(n, nc))
                            if mode == 3:
                                x = x + 1j*rng.normal(size=(n, nc))
                            noverlap, nw, kmax, nfft, fs = nperseg//2, 3.5, 6, nperseg+3, 8.0
                            ridge = 0.035
                            axis = 1 if len(cases) % 2 else 2
                            params = (x, nw, nperseg, noverlap, kmax, nfft, fs,
                                      weighting, trend, axis)
                            cases.append((mode, ridge, alpha, params, None))
    return cases


def format_input(cases):
    stream = []
    for mode, ridge, alpha, payload, _ in cases:
        if mode < 2:
            s = payload
            nf, nt, nc, _ = s.shape
            stream.append(f'{mode}\n{nf} {nt} {nc} {ridge:.17g} {alpha:.17g}\n')
            for t in range(nt):
                for f in range(nf):
                    for j in range(nc):
                        for i in range(nc):
                            v = s[f, t, i, j]
                            stream.append(f'{v.real:.17g} {v.imag:.17g}\n')
        else:
            x, nw, nperseg, noverlap, kmax, nfft, fs, weighting, trend, axis = payload
            n, nc = x.shape
            stream.append(f'{mode}\n{n} {nc} {nperseg} {noverlap} {nw:.17g} '
                          f'{kmax} {nfft} {fs:.17g} {ridge:.17g} {alpha:.17g} '
                          f'{axis} {weighting} {trend}\n')
            for j in range(nc):
                for i in range(n):
                    v = x[i, j]
                    stream.append(f'{v.real:.17g} {v.imag:.17g}\n')
    return ''.join(stream)


def check(executable):
    cases = build_cases()
    p = subprocess.run([executable], input=format_input(cases),
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       check=True)
    lines = iter(p.stdout.splitlines())
    valid_values = 0
    singular_cells = 0
    max_coherency_error = 0.0
    max_precision_relative_error = 0.0
    max_spectrum_error = 0.0
    for case_no, (mode, ridge, alpha, payload, singular) in enumerate(cases):
        status, nf, nt, nc = map(int, next(lines).split())
        if mode < 2:
            expected = payload
            freqs = np.arange(nf)
            times = np.arange(nt)
        else:
            x, nw, nperseg, noverlap, kmax, nfft, fs, weighting, trend, axis = payload
            expected, freqs, times = reference_spectrogram(
                x, nw, nperseg, noverlap, kmax, nfft, fs, weighting,
                trend, mode == 2)
        if expected.shape != (nf, nt, nc, nc):
            raise AssertionError(f'case {case_no}: unexpected output shape {(nf,nt,nc)}')
        if status not in (0, 2):
            raise AssertionError(f'case {case_no}: status={status}')
        for t in range(nt):
            for f in range(nf):
                f_out, t_out, is_valid = next(lines).split()
                np.testing.assert_allclose(float(f_out), freqs[f], rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(float(t_out), times[t], rtol=1e-12, atol=1e-12)
                is_valid = bool(int(is_valid))
                actual = np.empty((nc, nc), complex)
                coh = np.empty((nc, nc), float)
                c = np.empty((nc, nc), complex)
                precision = np.empty((nc, nc), complex)
                for j in range(nc):
                    for i in range(nc):
                        values = np.fromstring(next(lines), sep=' ')
                        if len(values) != 7:
                            raise AssertionError(f'case {case_no}: malformed output')
                        coh[i, j] = values[0]
                        c[i, j] = values[1] + 1j*values[2]
                        precision[i, j] = values[3] + 1j*values[4]
                        actual[i, j] = values[5] + 1j*values[6]
                max_spectrum_error = max(max_spectrum_error,
                                         np.max(np.abs(actual - expected[f,t])))
                if mode >= 2:
                    np.testing.assert_allclose(actual, expected[f,t], rtol=2e-9, atol=3e-11)
                if singular is not None and singular[f,t]:
                    if is_valid or not np.isnan(c).all() or not np.isnan(coh).all():
                        raise AssertionError(f'case {case_no}: singular cell not masked')
                    singular_cells += 1
                    continue
                if not is_valid:
                    raise AssertionError(f'case {case_no}: unexpected singular cell {(f,t)}')
                rc, rcomplex, rp = reference_partial(expected[f,t],ridge,alpha)
                np.testing.assert_allclose(coh,rc,atol=2e-9,rtol=2e-8)
                np.testing.assert_allclose(c,rcomplex,atol=2e-9,rtol=2e-8)
                np.testing.assert_allclose(precision,rp,atol=1e-7,rtol=2e-8)
                np.testing.assert_allclose(coh,np.abs(c)**2,atol=2e-13,rtol=2e-13)
                np.testing.assert_allclose(c,c.conj().T,atol=2e-13,rtol=2e-13)
                max_coherency_error = max(max_coherency_error, np.max(np.abs(c-rcomplex)))
                max_precision_relative_error = max(max_precision_relative_error,
                    np.max(np.abs(precision-rp)) / max(1.0,np.max(np.abs(rp))))
                valid_values += nc*nc
    if next(lines, None) is not None:
        raise AssertionError('unexpected trailing output')
    print(f'PASS {len(cases)} cases, {valid_values} complex coherency entries, '
          f'{singular_cells} expected singular cells; '
          f'max coherency error {max_coherency_error:.4g}; '
          f'max relative precision error {max_precision_relative_error:.4g}; '
          f'max spectral-matrix error {max_spectrum_error:.4g}', flush=True)


if __name__ == '__main__':
    check(sys.argv[1])
