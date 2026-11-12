"""Compare SciFort multitaper PSD with independent SciPy DPSS + NumPy FFT.

SciPy provides DPSS tapers but not this exact multitaper API.
Usage: python test/compare_signal_multitaper.py [gfortran]
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal.windows import dpss

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / s for s in
           ('kinds.f90', 'signal_dpss.f90', 'signal_windows_extra.f90', 'signal_spectral.f90', 'signal_multitaper.f90')]
DRIVER = ROOT / 'reference' / 'signal_multitaper_reference_driver.f90'


def reference_spectrum(x: np.ndarray, nw: float, k: int,
                       fs: float, nf: int, center: bool, eigen: bool) -> tuple[np.ndarray, np.ndarray]:
    tapers, ratios = dpss(len(x), nw, k, sym=False, return_ratios=True)
    if center:
        x = x - np.mean(x)
    bins = np.fft.rfftfreq(nf, d=1/fs)
    powers = np.abs(np.fft.rfft(tapers * x[None,:], n=nf, axis=-1))**2
    powers /= fs * np.sum(tapers**2, axis=-1)[:,None]
    weights = ratios if eigen else np.ones_like(ratios)
    combined = np.average(powers, axis=0, weights=weights)
    if nf % 2 == 0:
        combined[1:-1] *= 2
    else:
        combined[1:] *= 2
    return bins, combined


def main() -> None:
    compiler = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('FC', 'gfortran')
    rng = np.random.default_rng(20261020)
    cases = []
    inputs = []
    for i in range(240):
        n = int(rng.integers(8, 90))
        nw = float(rng.uniform(0.5, min(6.5, n/2-0.1)))
        k = int(rng.integers(1, min(n, 7)+1))
        nf = int(rng.integers(n, n+60))
        fs = float(rng.uniform(0.1, 200))
        center = bool(i % 2)
        eigen = bool((i // 2) % 2)
        x = rng.normal(size=n)
        if i % 6 == 0:
            x = np.sin(np.arange(n)*2*np.pi*rng.uniform(0.05,0.45)) + .05*x
        if i % 20 == 0:
            x[:] = 1.3
        cases.append((x, nw, k, fs, nf, center, eigen))
        inputs.append(f'{n} {nw:.17g} {k} {fs:.17g} {nf} {int(center)} {int(eigen)}\n')
        inputs.append(' '.join(f'{v:.17g}' for v in x) + '\n')
    with tempfile.TemporaryDirectory(prefix='scifort-multitaper-') as directory:
        tmp = Path(directory)
        exe = tmp / ('multitaper.exe' if os.name == 'nt' else 'multitaper')
        subprocess.run([compiler, '-std=f2018', '-O2', '-fcheck=all',
                        '-Werror=implicit-interface', *map(str, SOURCES),
                        str(DRIVER), '-o', str(exe)], cwd=tmp, check=True)
        proc = subprocess.run([str(exe)], input=''.join(inputs), capture_output=True,
                              text=True, check=True)
    lines = iter(proc.stdout.splitlines())
    max_error = max_freq_error = 0.0
    for i,(x,nw,k,fs,nf,center,eigen) in enumerate(cases):
        header = next(lines).split()
        assert header == ['CASE', str(nf//2+1), '0'], (i,header)
        freq = np.fromstring(next(lines),sep=' ')
        power = np.fromstring(next(lines),sep=' ')
        expected_freq, expected_power = reference_spectrum(x,nw,k,fs,nf,center,eigen)
        err = float(np.max(np.abs(power-expected_power)))
        ferr = float(np.max(np.abs(freq-expected_freq)))
        max_error = max(max_error,err)
        max_freq_error = max(max_freq_error,ferr)
        if not np.allclose(power,expected_power,atol=2e-10,rtol=2e-8):
            raise AssertionError(f'case {i}: max absolute power error {err:.3g}')
        if ferr > 1e-12:
            raise AssertionError(f'case {i}: max frequency error {ferr:.3g}')
    print(f'Multitaper SciPy/NumPy comparisons: {len(cases)} passed, '
          f'max abs PSD error {max_error:.3g}, max frequency error {max_freq_error:.3g}')


if __name__ == '__main__':
    main()
