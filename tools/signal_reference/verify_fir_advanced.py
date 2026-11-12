#!/usr/bin/env python3
"""Cross-check advanced SciFort FIR routines against installed SciPy.

Requires SciPy, NumPy, and gfortran. Execute from any working directory:
    python tools/signal_reference/verify_fir_advanced.py
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261017


def line(values: np.ndarray | list[float]) -> str:
    return ' '.join(f'{v:.17e}' for v in values)


def make_cases() -> tuple[str, list[tuple[str, np.ndarray]]]:
    rng = np.random.default_rng(SEED)
    input_lines: list[str] = []
    expected: list[tuple[str, np.ndarray]] = []

    # Types I-IV, arbitrary gain curves, discontinuities, sample rates,
    # window choices, FFT-friendly and non-power-of-two sampling meshes.
    for _ in range(400):
        n = int(rng.integers(1, 58))
        anti = bool(rng.integers(0, 2))
        fs = float(rng.choice([2, 48, 1000]))
        window = str(rng.choice(['hamming', 'hann', 'blackman', 'boxcar']))
        nedges = int(rng.integers(2, 8))
        inner = np.sort(rng.uniform(.07, .93, size=nedges-2))
        freq = np.r_[0, inner, 1.0]*fs/2
        if nedges > 3 and bool(rng.integers(0, 4) == 0):
            k = int(rng.integers(1, nedges-2))
            freq[k+1] = freq[k]
        gain = rng.uniform(0, 1.7, size=nedges)
        if anti:
            gain[0] = 0.0
            if n%2:
                gain[-1] = 0.0
        elif n%2 == 0:
            gain[-1] = 0.0
        default_grid = 2**int(np.ceil(np.log2(n))) + 1
        nfreq = int(rng.choice([default_grid, max(n+2, default_grid+2), default_grid+4]))
        taps = signal.firwin2(n, freq, gain, nfreqs=nfreq,
                              window=window, antisymmetric=anti, fs=fs)
        expected.append(('firwin2', taps))
        input_lines += ['WIN2',f'{n} {nedges} {nfreq} {int(anti)} {fs:.17e}',
                        window,line(freq),line(gain)]

    # Multi-band weighted least squares, with flat and tilted desired gains.
    accepted_ls = 0
    while accepted_ls < 340:
        n = int(rng.integers(0, 8))*2+1
        fs = float(rng.choice([2, 40, 800]))
        nbands = int(rng.integers(1, 4))
        breaks = np.linspace(0, 1, 2*nbands+1)
        breaks[1:-1] += rng.uniform(-0.035, .035, size=2*nbands-1)
        breaks[0],breaks[-1] = 0, 1
        bands = np.column_stack([breaks[0:-1:2],breaks[1::2]]).ravel()*fs/2
        desired = rng.uniform(0, 1.5, size=2*nbands)
        weights = rng.uniform(.3, 3.5, size=nbands)
        # Filter out severely ill-conditioned normal equations: the two
        # double-precision solvers need not agree on such coefficient vectors.
        m = (n-1)//2
        q = np.zeros(n)
        for j in range(nbands):
            lo, hi = bands[2*j:2*j+2]/(fs/2)
            q[0] += weights[j]*(hi-lo)
            k = np.arange(1,n)
            q[1:] += weights[j]*(np.sin(np.pi*k*hi)-np.sin(np.pi*k*lo))/(np.pi*k)
        indices = np.arange(m+1)
        mat = q[np.abs(indices[:,None]-indices[None,:])] + q[indices[:,None]+indices[None,:]]
        if np.linalg.cond(mat) > 1e6:
            continue
        taps = signal.firls(n, bands, desired, weight=weights, fs=fs)
        accepted_ls += 1
        expected.append(('firls', taps))
        input_lines += ['LS',f'{n} {len(bands)} {fs:.17e}',
                        line(bands),line(desired),line(weights)]

    # Cepstral homomorphic and modified discrete Hilbert conversion.
    for _ in range(330):
        n = int(rng.integers(5, 44))
        cutoff = float(rng.uniform(.12, .81))
        h = signal.firwin(n, cutoff, window='hamming')
        method = str(rng.choice(['homomorphic', 'hilbert']))
        half = method == 'hilbert' or bool(rng.integers(0, 2))
        n_fft = int(rng.choice([512, 1024, 2048, 4096]))
        taps = signal.minimum_phase(h, method=method, half=half, n_fft=n_fft)
        expected.append((f'minimum_phase_{method}', taps))
        input_lines += ['MIN',f'{n} {n_fft} {int(half)}',method,line(h)]

    return '\n'.join(input_lines)+'\n',expected


def main() -> None:
    data, cases = make_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-fir-advanced-') as temp:
        work = Path(temp)
        flags = ['-std=f2018','-O2','-Wimplicit-interface',
                 '-Werror=implicit-interface','-fcheck=all']
        objects: list[str] = []
        for name in ['kinds','signal_dpss', 'signal_windows_extra', 'signal_spectral','signal_fir_advanced']:
            source = ROOT / 'src' / f'{name}.f90'
            out = work / f'{name}.o'
            subprocess.run(['gfortran',*flags,'-I',str(work),'-J',str(work),
                            '-c',str(source),'-o',str(out)],check=True,cwd=work)
            objects.append(str(out))
        exe = work / 'fir_advanced_driver'
        subprocess.run(['gfortran',*flags,'-I',str(work),'-J',str(work),
                        str(ROOT/'tools/signal_reference/fir_advanced_driver.f90'),
                        *objects,'-o',str(exe)],check=True,cwd=work)
        process = subprocess.run([str(exe)],input=data,stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,text=True,check=True,cwd=work)
    tokens = iter(process.stdout.split())
    max_errors: dict[str,float] = {}
    checked: dict[str,int] = {}
    for case_index,(name,expected_taps) in enumerate(cases):
        try:
            status = int(next(tokens))
            n = int(next(tokens))
            actual = np.array([float(next(tokens)) for _ in range(n)])
        except (ValueError,StopIteration) as err:
            raise AssertionError(f'Cannot read result at case {case_index}') from err
        if status != 0 or len(actual) != len(expected_taps):
            raise AssertionError(f'{name} case {case_index}: status {status}, length {n} '
                                 f'(expected {len(expected_taps)})')
        error = float(np.max(np.abs(actual-expected_taps)))
        scale = 1+float(np.max(np.abs(expected_taps)))
        tolerance = {'firwin2':8e-12,'firls':8e-9,'minimum_phase_hilbert':1e-2,'minimum_phase_homomorphic':2e-9}[name]
        if not np.all(np.isfinite(actual)) or error > tolerance*scale:
            raise AssertionError(f'{name} case {case_index}: max error {error:.3e} '
                                 f'over tolerance {tolerance*scale:.3e}\n'
                                 f'actual={actual}\nexpected={expected_taps}')
        checked[name] = checked.get(name,0)+1
        max_errors[name] = max(error,max_errors.get(name,0.0))
    if list(tokens):
        raise AssertionError('Extra output from Fortran reference driver')
    print(f'PASS: {len(cases)} randomized comparisons against SciPy {__import__("scipy").__version__}')
    for name in ['firwin2','firls','minimum_phase_homomorphic','minimum_phase_hilbert']:
        print(f'  {name}: {checked[name]} cases, max absolute error {max_errors[name]:.3e}')


if __name__ == '__main__':
    main()
