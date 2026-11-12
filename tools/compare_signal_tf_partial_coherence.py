#!/usr/bin/env python3
"""Independent time-frequency partial-coherence checks using SciPy DPSS and NumPy.

Usage: python tools/compare_signal_tf_partial_coherence.py /path/to/fortran_driver
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy.signal.windows import dpss


def reference_spectrogram(x, nw, nperseg, noverlap, kmax, nfft, fs,
                          weighting, trend, is_real):
    """Compute sliding multitaper CSD matrices without using SciFort results."""
    n, nc = x.shape
    tapers, eigenvalues = dpss(nperseg, nw, Kmax=kmax, sym=False,
                              norm=2, return_ratios=True)
    norm = np.sqrt(fs * np.sum(tapers ** 2, axis=1))
    q = np.ones(kmax) if weighting == 'equal' else eigenvalues.copy()
    q /= q.sum()
    hop = nperseg - noverlap
    starts = range(0, n - nperseg + 1, hop)
    frames = []
    for start in starts:
        seg = x[start:start + nperseg, :].copy()
        if trend in ('constant', 'linear'):
            seg -= seg.mean(axis=0)
        if trend == 'linear':
            grid = np.arange(nperseg, dtype=float) - (nperseg - 1) / 2
            slope = np.sum(grid[:, None] * seg, axis=0) / np.sum(grid ** 2)
            seg -= grid[:, None] * slope
        tapered = tapers[:, :, None] * seg[None, :, :]
        fft = np.fft.rfft if is_real else np.fft.fft
        coeff = fft(tapered, n=nfft, axis=1) / norm[:, None, None]
        mat = np.einsum('kfc,kfd,k->fcd', coeff.conj(), coeff, q,
                        optimize=True)
        if is_real:
            if nfft % 2 == 0:
                mat[1:-1, :, :] *= 2
            else:
                mat[1:, :, :] *= 2
        frames.append(mat)
    matrices = np.stack(frames, axis=1)   # frequency,time,channel,channel
    freq = (np.fft.rfftfreq if is_real else np.fft.fftfreq)(nfft, 1 / fs)
    times = (np.array(list(starts)) + nperseg // 2) / fs
    return matrices, freq, times


def reference_coherence(s, ridge):
    """Directly invert loaded matrices independently at each time/frequency."""
    nc = s.shape[0]
    regularized = s + np.eye(nc) * (ridge * np.trace(s).real / nc)
    precision = np.linalg.inv(regularized)
    diagonal = np.diag(precision).real
    coherence = np.abs(precision) ** 2 / (diagonal[:, None] * diagonal[None, :])
    np.fill_diagonal(coherence, 1.0)
    return coherence, precision


def build_cases():
    rng = np.random.default_rng(66713)
    cases = []
    for nc in (2, 3, 4, 5):
        for ridge in (0.0, 0.015, 0.15):
            for trial in range(6):
                nf, nt = 2 + trial % 2, 2 + trial % 3
                matrix = np.empty((nf, nt, nc, nc), dtype=np.complex128)
                for f in range(nf):
                    for t in range(nt):
                        nobs = nc + 3
                        z = (rng.normal(size=(nobs, nc))
                             + 1j * rng.normal(size=(nobs, nc)))
                        if trial == 2 and (f + t) % 3 == 0:
                            z[:, -1] = z[:, -2]
                        if trial == 4 and f == 0 and t == 0:
                            z[:, :] = z[:, :1]
                        s = z.conj().T @ z
                        matrix[f, t] = s * 10.0 ** int(rng.integers(-3, 4))
                cases.append(('matrix', ridge, matrix))
    for is_real in (True, False):
        for nc in (2, 3, 4):
            for nperseg in (24, 32):
                for weighting in ('equal', 'eigen'):
                    for trend in ('none', 'linear'):
                        n = 2 * nperseg + 7
                        x = rng.normal(size=(n, nc))
                        if not is_real:
                            x = x + 1j * rng.normal(size=(n, nc))
                        noverlap = nperseg // 2
                        nw = 3.5
                        kmax = 7
                        nfft = nperseg + 3
                        fs = 4.0
                        ridge = 0.03
                        axis = 1 if (len(cases) % 2) else 2
                        cases.append(('signal', ridge, (x, nw, nperseg, noverlap,
                                     kmax, nfft, fs, weighting, trend, is_real, axis)))
    return cases


def run(executable):
    cases = build_cases()
    stream = []
    for kind, ridge, content in cases:
        if kind == 'matrix':
            mat = content
            nf, nt, nc, _ = mat.shape
            stream.append(f'0\n{nf} {nt} {nc} {ridge:.17g}\n')
            for t in range(nt):
                for f in range(nf):
                    for j in range(nc):
                        for i in range(nc):
                            v = mat[f, t, i, j]
                            stream.append(f'{v.real:.17g} {v.imag:.17g}\n')
        else:
            (x, nw, nperseg, noverlap, kmax, nfft, fs,
             weighting, trend, is_real, axis) = content
            n, nc = x.shape
            stream.append(f'1\n{n} {nc} {nperseg} {noverlap} {nw:.17g} '
                          f'{kmax} {nfft} {fs:.17g} {ridge:.17g} '
                          f'{int(is_real)} {axis} {weighting} {trend}\n')
            for j in range(nc):
                for i in range(n):
                    v = x[i, j]
                    stream.append(f'{v.real:.17g} {v.imag:.17g}\n')
    p = subprocess.run([executable], input=''.join(stream), text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    lines = iter(p.stdout.splitlines())
    max_coherence_error = 0.0
    max_spectrum_error = 0.0
    values = 0
    singular = 0
    for case_num, (kind, ridge, content) in enumerate(cases):
        status, nf, nt, nc = map(int, next(lines).split())
        if kind == 'matrix':
            expected = content
            freq = np.arange(expected.shape[0])
            times = np.arange(expected.shape[1])
        else:
            (x, nw, nperseg, noverlap, kmax, nfft, fs,
             weighting, trend, is_real, axis) = content
            expected, freq, times = reference_spectrogram(x, nw, nperseg,
                noverlap, kmax, nfft, fs, weighting, trend, is_real)
        if (nf, nt, nc) != (*expected.shape[:2], expected.shape[2]):
            raise AssertionError(f'case {case_num} wrong shapes: {nf,nt,nc}, {expected.shape}')
        if status not in (0, 2):
            raise AssertionError(f'case {case_num} unexpected status {status}')
        for t in range(nt):
            for f in range(nf):
                freq_out, time_out, valid_text = next(lines).split()
                np.testing.assert_allclose(float(freq_out), freq[f], atol=1e-12, rtol=1e-13)
                np.testing.assert_allclose(float(time_out), times[t], atol=1e-12, rtol=1e-13)
                valid = bool(int(valid_text))
                coh = np.empty((nc,nc), float)
                precision = np.empty((nc,nc), complex)
                actual = np.empty((nc,nc), complex)
                for j in range(nc):
                    for i in range(nc):
                        numbers = np.fromstring(next(lines), sep=' ')
                        if len(numbers) != 5:
                            raise AssertionError(f'case {case_num} malformed output')
                        coh[i,j] = numbers[0]
                        precision[i,j] = numbers[1] + 1j*numbers[2]
                        actual[i,j] = numbers[3] + 1j*numbers[4]
                max_spectrum_error = max(max_spectrum_error,
                    float(np.max(np.abs(actual - expected[f,t]))))
                if kind == 'signal':
                    np.testing.assert_allclose(actual,expected[f,t],atol=3e-11,rtol=1e-9)
                if not valid:
                    singular += 1
                    if not np.isnan(coh).all():
                        raise AssertionError('singular time-frequency bins must contain NaNs')
                    continue
                ref_coh, ref_precision = reference_coherence(expected[f,t],ridge)
                max_coherence_error = max(max_coherence_error,
                    float(np.max(np.abs(coh-ref_coh))))
                np.testing.assert_allclose(coh,ref_coh,atol=1e-9,rtol=3e-8)
                np.testing.assert_allclose(precision,ref_precision,atol=4e-7,rtol=3e-8)
                values += nc*nc
    if next(lines,None) is not None:
        raise AssertionError('unexpected trailing output')
    print(f'PASS {len(cases)} cases, {values} valid conditional-coherence values, '
          f'{singular} singular time-frequency cells; '
          f'max coherence error {max_coherence_error:.4g}; '
          f'max CSD matrix error {max_spectrum_error:.4g}', flush=True)


if __name__ == '__main__':
    run(sys.argv[1])
