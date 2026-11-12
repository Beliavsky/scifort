#!/usr/bin/env python3
"""Validate SciFort independent-axis and elliptical 2D FIR designs.

SciFort's independent per-axis cutoffs and rotated elliptical radial windows
are extensions, not scipy.signal.firwin_2d options. This script constructs
independent references from scipy.signal.firwin and NumPy's interpolation.

Run: python tools/signal_reference/compare_signal_fir2d_anisotropic_scipy.py
Use --optimized for a separate -O3 build. Requires NumPy, SciPy and gfortran.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    'kinds', 'signal_windows_extra', 'signal_dpss', 'signal_spectral',
    'signal_analog_design', 'signal_design', 'signal_fir2d',
)


def separable_reference(mode: int, nr: int, nc: int) -> np.ndarray:
    first = {'window': 'hamming', 'cutoff': .23}
    second = {'window': 'hann', 'cutoff': .37}
    if mode == 1:
        first = {'window': 'hann', 'cutoff': .17}
        second = {'window': 'blackman', 'cutoff': .36, 'fs': 2.8}
    elif mode == 2:
        first = {'window': 'blackman', 'cutoff': .39, 'pass_zero': 'highpass'}
        second = {'window': 'hamming', 'cutoff': .72, 'fs': 3.1,
                  'pass_zero': 'highpass'}
    elif mode == 3:
        first = {'window': 'hamming', 'cutoff': [.12, .29], 'pass_zero': 'bandpass'}
        second = {'window': 'hann', 'cutoff': [.19, .37], 'pass_zero': 'bandstop'}
    elif mode == 4:
        first = {'window': ('kaiser', 6.), 'cutoff': [.13, .35],
                 'pass_zero': 'bandstop'}
        second = {'window': ('kaiser', 3.), 'cutoff': [.20, .48],
                  'pass_zero': 'bandpass', 'scale': False}
    elif mode == 5:
        first = {'window': 'hann', 'cutoff': .16, 'fs': 2.5}
        second = {'window': 'hann', 'cutoff': [.21, .61],
                  'fs': 3.1, 'pass_zero': 'bandpass'}
    elif mode == 6:
        first = {'window': 'hamming', 'cutoff': [.08, .33],
                 'pass_zero': 'bandpass', 'scale': False}
        second = {'window': 'blackman', 'cutoff': .43,
                  'pass_zero': 'highpass', 'scale': True}
    elif mode != 7:
        raise ValueError('invalid separable test mode')
    return np.outer(signal.firwin(nr, **first), signal.firwin(nc, **second))


def elliptical_reference(mode: int, nr: int, nc: int) -> np.ndarray:
    settings = {
        1: ((1.4, .65), 0., 'hamming', .18, 2., 'lowpass'),
        2: ((.65, 1.8), .37, 'blackman', .27, 2., 'lowpass'),
        3: ((.9, 1.2), -.58, 'hamming', [.16, .29], 2., 'bandpass'),
        4: ((1., 1.), 1.2, 'hann', .26, 3., 'lowpass'),
    }
    axes, angle, window, cutoff, fs, pass_zero = settings[mode]
    nrad = 8 * max(nr, nc)
    radial_profile = signal.firwin(nrad, cutoff, window=window, fs=fs,
                                   pass_zero=pass_zero)
    rows, cols = np.meshgrid(np.linspace(-1., 1., nr),
                             np.linspace(-1., 1., nc))
    u = (np.cos(angle) * rows + np.sin(angle) * cols) / axes[0]
    v = (-np.sin(angle) * rows + np.cos(angle) * cols) / axes[1]
    distances = np.hypot(u, v)
    return np.interp(distances, np.linspace(0., 1., nrad), radial_profile)


def compare(output: str) -> tuple[int, int, float, float]:
    count = 0
    values = 0
    max_absolute = 0.
    max_scaled = 0.
    for line in output.splitlines():
        parts = line.split()
        family = parts[0]
        case, nr, nc, status, got_nr, got_nc = map(int, parts[1:7])
        if status != 0:
            raise AssertionError(f'{family} case {case}, ({nr},{nc}): status {status}')
        actual = np.asarray(parts[7:], dtype=np.float64).reshape((got_nr, got_nc), order='F')
        expected = (separable_reference(case, nr, nc) if family == 'S' else
                    elliptical_reference(case, nr, nc))
        if actual.shape != expected.shape:
            raise AssertionError(f'{family} case {case}: shape {actual.shape} != {expected.shape}')
        difference = np.abs(actual - expected)
        absolute = float(np.max(difference, initial=0.))
        scaled = float(np.max(difference / np.maximum(1., np.abs(expected)), initial=0.))
        if not np.isfinite(scaled) or scaled > 2e-11:
            raise AssertionError(f'{family} {case} ({nr},{nc}) difference {absolute:.4e}')
        count += 1
        values += actual.size
        max_absolute = max(max_absolute, absolute)
        max_scaled = max(max_scaled, scaled)
    return count, values, max_absolute, max_scaled


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimized', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort-fir2d-anisotropic-') as folder:
        work = Path(folder)
        flags = ['-std=f2018', '-Werror=implicit-interface',
                 '-O3' if args.optimized else '-O0',
                 '-I', str(work), '-J', str(work)]
        if not args.optimized:
            flags.extend(['-g', '-fcheck=all', '-fbacktrace'])
        objects = []
        for name in SOURCES:
            source = ROOT / 'src' / (name + '.f90')
            obj = work / (name + '.o')
            subprocess.run(['gfortran', *flags, '-c', str(source), '-o', str(obj)],
                           check=True, cwd=work)
            objects.append(str(obj))
        exe = work / 'fir2d_anisotropic_driver'
        subprocess.run(['gfortran', *flags,
                        str(ROOT / 'tools/signal_reference' / 'fir2d_anisotropic_driver.f90'),
                        *objects, '-o', str(exe)], check=True, cwd=work)
        result = subprocess.check_output([str(exe)], text=True)
        ncases, nvalues, absolute, scaled = compare(result)
        print(f'PASS: {ncases} anisotropic and elliptical FIR2D SciPy comparisons, '
              f'{nvalues} coefficients, maximum absolute {absolute:.3e}, '
              f'magnitude-scaled {scaled:.3e}')


if __name__ == '__main__':
    main()
