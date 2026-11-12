"""Independent SciPy DPSS / NumPy FFT validation for SciFort's spectral matrices.

Run from any directory with: python test/compare_signal_multitaper_channels.py
The comparison requires NumPy, SciPy and a local Fortran compiler only for testing.
"""
from __future__ import annotations

from pathlib import Path
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal.windows import dpss

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / name for name in (
    'kinds.f90', 'signal_dpss.f90', 'signal_windows_extra.f90', 'signal_spectral.f90', 'signal_multitaper.f90')]
DRIVER = ROOT / 'reference' / 'signal_multitaper_channels_reference_driver.f90'


def oracle(samples, nw, kmax, fs, nf, center, eigen):
    """Return frequency, 3-D cross-PSD using independent vectorized FFTs."""
    n, nc = samples.shape
    tapers, ratios = dpss(n, nw, Kmax=kmax, sym=False, return_ratios=True)
    if center:
        samples = samples - np.mean(samples, axis=0, keepdims=True)
    # ft[taper, frequency, channel]
    windowed = tapers[:, :, None] * samples[None, :, :]
    real_input = not np.iscomplexobj(samples)
    ft = (np.fft.rfft(windowed, n=nf, axis=1) if real_input else
          np.fft.fft(windowed, n=nf, axis=1))
    weights = ratios if eigen else np.ones(kmax)
    weights = weights / (fs * np.sum(tapers ** 2, axis=1))
    power = np.einsum('kfi,kfj,k->fij', np.conjugate(ft), ft, weights)
    power /= np.sum(ratios if eigen else np.ones(kmax))
    if real_input:
        frequencies = np.fft.rfftfreq(nf, d=1 / fs)
        if nf % 2 == 0:
            power[1:-1] *= 2
        else:
            power[1:] *= 2
    else:
        frequencies = np.fft.fftfreq(nf, d=1 / fs)
    return frequencies, power


def main():
    compiler = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('FC', 'gfortran')
    rng = np.random.default_rng(20261020)
    cases = []
    request = []
    for i in range(160):
        n = int(rng.integers(10, 65))
        nc = int(rng.integers(1, 5))
        nf = n + int(rng.integers(0, 32))
        nw = float(rng.uniform(0.7, min(5.5, n/2-0.1)))
        kmax = int(rng.integers(1, min(6, n) + 1))
        fs = float(rng.uniform(1.0, 200.0))
        is_complex = bool(i % 2)
        ax = 1 if i % 4 < 2 else 2
        center = bool((i // 2) % 2)
        eigen = bool((i // 4) % 2)
        samples = rng.normal(size=(n,nc))
        if is_complex:
            samples = samples + 1j * rng.normal(size=(n,nc))
        if i % 17 == 0:
            samples[:, -1] = samples[:, 0] * (2-0.7j if is_complex else 1.8)
        cases.append((samples,nw,kmax,fs,nf,center,eigen))
        request.append(f'{n} {nc} {nf} {kmax} {nw:.17g} {fs:.17g} '
                       f'{int(is_complex)} {ax} {int(center)} {int(eigen)}\n')
        for channel in range(nc):
            for t in range(n):
                v = samples[t,channel]
                request.append(f'{v.real:.17g} {v.imag:.17g}\n')

    with tempfile.TemporaryDirectory(prefix='scifort-mt-matrix-') as directory:
        tmp = Path(directory)
        exe = tmp / ('spectra.exe' if os.name == 'nt' else 'spectra')
        subprocess.run([compiler, '-std=f2018', '-O2', '-fcheck=all',
                        '-Werror=implicit-interface', *map(str,SOURCES),
                        str(DRIVER), '-o', str(exe)], cwd=tmp, check=True)
        proc = subprocess.run([str(exe)], input=''.join(request),
                              capture_output=True, text=True, check=True)

    lines = iter(proc.stdout.splitlines())
    max_abs = max_freq = 0.0
    for index, (samples,nw,kmax,fs,nf,center,eigen) in enumerate(cases):
        nc = samples.shape[1]
        freq_ref, spec_ref = oracle(samples,nw,kmax,fs,nf,center,eigen)
        header = next(lines).split()
        assert header == ['CASE',str(len(freq_ref)),str(nc),'0'], (index,header)
        freq = np.fromstring(next(lines),sep=' ')
        spec = np.empty_like(spec_ref)
        for y in range(nc):
            for x in range(nc):
                raw = np.fromstring(next(lines),sep=' ')
                assert len(raw) == 2*len(freq_ref)
                spec[:,x,y] = raw[::2] + 1j*raw[1::2]
        err = float(np.max(np.abs(spec-spec_ref)))
        ferr = float(np.max(np.abs(freq-freq_ref)))
        max_abs = max(err,max_abs)
        max_freq = max(max_freq,ferr)
        if not np.allclose(spec,spec_ref,atol=3e-10,rtol=3e-8):
            raise AssertionError(f'case {index}, diff {err:.3g}, '
                                 f'expected {np.max(abs(spec_ref)):.3g}')
        if ferr > 1e-12:
            raise AssertionError(f'case {index}: frequency mismatch {ferr:.3g}')
        # Additional oracle invariants: Hermitian and positive semidefinite.
        if not np.allclose(spec, np.conjugate(np.swapaxes(spec,1,2)),atol=1e-10):
            raise AssertionError(f'case {index}: matrix not Hermitian')
        if np.min(np.linalg.eigvalsh(spec)) < -1e-9:
            raise AssertionError(f'case {index}: spectral matrix not PSD')
    print(f'Multitaper real/complex multichannel SciPy comparisons: {len(cases)} passed; '
          f'max absolute matrix error {max_abs:.3g}, frequency error {max_freq:.3g}')


if __name__ == '__main__':
    main()
