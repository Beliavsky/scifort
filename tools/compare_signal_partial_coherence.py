#!/usr/bin/env python3
"""Independent NumPy and SciPy references for spectral partial coherence.

Usage: python tools/compare_signal_partial_coherence.py /path/to/reference_partial
"""
from __future__ import annotations

import sys
import subprocess
import numpy as np
from scipy.signal.windows import dpss


def reference_matrix(s, ridge):
    """Inverse-spectral precision and squared conditional coherence."""
    c = s.shape[0]
    a = s + np.eye(c) * (ridge * np.trace(s).real / c)
    precision = np.linalg.inv(a)
    d = np.real(np.diag(precision))
    coherence = np.abs(precision) ** 2 / (d[:, None] * d[None, :])
    np.fill_diagonal(coherence, 1.0)
    return precision, coherence


def reference_multitaper(x, nw, kmax, nfft, fs, weighting, trend, is_real):
    """Compute the CSD matrix directly from SciPy's DPSS and NumPy FFTs."""
    n, nc = x.shape
    tapers, ratios = dpss(n, nw, Kmax=kmax, sym=False, norm=2, return_ratios=True)
    data = x.copy()
    if trend == 'constant':
        data = data - data.mean(axis=0)
    elif trend == 'linear':
        grid = np.arange(n, dtype=float) - (n - 1) / 2
        data = data - data.mean(axis=0)
        data = data - grid[:, None] * np.sum(grid[:, None] * data, axis=0) / np.sum(grid**2)
    if is_real:
        spectrum = np.fft.rfft(tapers[:, :, None] * data[None, :, :], n=nfft, axis=1)
    else:
        spectrum = np.fft.fft(tapers[:, :, None] * data[None, :, :], n=nfft, axis=1)
    if weighting == 'equal':
        w = np.ones(kmax)
    else:
        w = ratios
    cross = np.einsum('kfc,kfd,k->fcd', np.conj(spectrum), spectrum,
                      w / (fs * np.sum(tapers**2, axis=1) * w.sum()), optimize=True)
    if is_real:
        if nfft % 2 == 0:
            cross[1:-1] *= 2
        else:
            cross[1:] *= 2
    return cross


def build_cases():
    rng = np.random.default_rng(4762)
    cases = []
    for nc in [2, 3, 4, 6]:
        for ridge in [0.0, 1e-3, 0.05, 0.3]:
            for trial in range(9):
                nf = int(rng.integers(1, 5))
                matrix = []
                for f in range(nf):
                    nobs = nc + 3 if trial % 4 else max(1, nc - 2)
                    z = rng.normal(size=(nobs, nc)) + 1j * rng.normal(size=(nobs, nc))
                    if trial == 4:
                        # Exactly rank-deficient spectrum, with identical channels.
                        z[:, -1] = z[:, -2]
                    s = z.conj().T @ z
                    scale = 10.0 ** int(rng.integers(-5, 6))
                    matrix.append(s * scale)
                cases.append(('matrix', ridge, np.array(matrix)))
    for is_real in [True, False]:
        for nc in [2,3,4]:
            for nfft in [32,33,48]:
                for weighting in ['equal','eigen']:
                    for trend in ['none','constant','linear']:
                        n = 32
                        x = rng.normal(size=(n,nc))
                        if not is_real:
                            x = x + 1j*rng.normal(size=(n,nc))
                        cases.append(('signal', .05, (x, 3.5, 7, nfft, 4., weighting, trend, is_real)))
    return cases


def run(executable):
    cases = build_cases()
    stream = []
    for kind, ridge, content in cases:
        if kind == 'matrix':
            s = content
            nf, nc, _ = s.shape
            stream.append(f'0\n{nf} {nc} {ridge:.17g}\n')
            for f in range(nf):
                for j in range(nc):
                    for i in range(nc):
                        z = s[f,i,j]
                        stream.append(f'{z.real:.17g} {z.imag:.17g}\n')
        else:
            x,nw,kmax,nfft,fs,weighting,trend,is_real = content
            n,nc=x.shape
            stream.append(f'1\n{n} {nc} {nw} {kmax} {nfft} {fs} {ridge} {int(is_real)} {weighting} {trend}\n')
            for j in range(nc):
                for i in range(n):
                    z = x[i,j]
                    stream.append(f'{z.real:.17g} {z.imag:.17g}\n')
    process = subprocess.run([executable], input=''.join(stream), text=True,
                             capture_output=True, check=True)
    lines = iter(process.stdout.splitlines())
    max_abs = 0.0
    max_spectral = 0.0
    checked_values = 0
    singular_bins = 0
    for kind,ridge,content in cases:
        status, nf, nc = [int(v) for v in next(lines).split()]
        if kind == 'matrix':
            expected = content
        else:
            expected = reference_multitaper(*content)
        if nf != len(expected) or nc != expected.shape[1]:
            raise AssertionError(f'{kind}: output shape {nf} x {nc} vs {expected.shape}')
        if status not in (0,2):
            raise AssertionError(f'Invalid result: {status}')
        for f in range(nf):
            freq, valid = next(lines).split()
            valid = bool(int(valid))
            coh = np.empty((nc,nc))
            precision = np.empty((nc,nc),complex)
            actual = np.empty((nc,nc),complex)
            for j in range(nc):
                for i in range(nc):
                    values = np.fromstring(next(lines),sep=' ')
                    if len(values) != 5:
                        raise AssertionError('Invalid output row')
                    coh[i,j] = values[0]
                    precision[i,j] = values[1]+1j*values[2]
                    actual[i,j] = values[3]+1j*values[4]
            if kind=='signal':
                max_spectral=max(max_spectral,float(np.max(np.abs(actual-expected[f]))))
                np.testing.assert_allclose(actual,expected[f],atol=1e-11,rtol=4e-10)
            # Compare conditional coherence against independently inverted
            # NumPy matrices; skip ill-conditioned, unregularized singular bins.
            if not valid:
                singular_bins += 1
                assert np.all(np.isnan(coh)), 'invalid frequency must be NaN'
                continue
            ref_p,ref_c=reference_matrix(expected[f],ridge)
            abs_err=float(np.max(np.abs(coh-ref_c)))
            max_abs=max(max_abs,abs_err)
            np.testing.assert_allclose(coh,ref_c,atol=3e-9,rtol=2e-8)
            np.testing.assert_allclose(precision,ref_p,atol=5e-7,rtol=3e-8)
            checked_values+=nc*nc
    if next(lines,None) is not None:
        raise AssertionError('Trailing unexpected driver output')
    print(f'PASS {len(cases)} independently checked matrices/data cases; '
          f'{checked_values} partial-coherence entries, {singular_bins} singular frequency bins; '
          f'max coherence error {max_abs:.3g}; max SciPy CSD matrix error {max_spectral:.3g}')


if __name__ == '__main__':
    run(sys.argv[1])
