#!/usr/bin/env python3
"""Check custom Fortran CWT kernels and peak positions against SciPy 1.17.

Run from any directory with Python, NumPy, SciPy and GNU Fortran installed.
Uses the private SciPy _cwt helper solely as a numerical test oracle; production
SciFort has no SciPy dependency.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import subprocess

import numpy as np
from scipy.signal import find_peaks_cwt
from scipy.signal._wavelets import _cwt

ROOT = Path(__file__).resolve().parents[2]


def wavelet(points, width, family):
    u = (np.arange(0, points, dtype=float) - (points-1.0)/2.0)/width
    if family == 1:
        return np.exp(-0.5*u*u)/(np.sqrt(2.0*np.pi)*width)
    if family == 2:
        u = u - 0.35
        return np.exp(-0.5*u*u)/(np.sqrt(2.0*np.pi)*width)
    if family == 3:
        return np.cos(3.3*u)*np.exp(-0.5*u*u)/np.sqrt(width)
    if family == 4:
        return u*np.exp(-0.5*u*u)/np.sqrt(width)
    raise ValueError(family)


def generate_cases(count):
    rng = np.random.default_rng(405892)
    cases = []
    for i in range(count):
        n = int(rng.integers(20, 145))
        nw = int(rng.integers(2, 18))
        family = 1 + i % 4
        widths = rng.uniform(.7, 13.0, size=nw)
        if i % 5 == 0:
            widths = np.linspace(.75, 13.0, nw)
        if i % 9 == 0:
            widths = rng.choice([.75, 1.25, 2.5, 4.25, 8.5], size=nw)
        x = rng.normal(size=n)*rng.uniform(.05, 1.5)
        for _ in range(int(rng.integers(0, 5))):
            center = rng.uniform(0, n)
            scale = rng.uniform(.7, 12)
            x += rng.uniform(.4, 6.)*np.exp(-.5*((np.arange(n)-center)/scale)**2)
        minlen = int(rng.integers(1, nw+1))
        window = int(rng.integers(2, 15))
        gap = float(rng.uniform(0.5, 6.0))
        snr = float(rng.uniform(0.1, 3.0))
        perc = float(rng.uniform(0, 95))
        cases.append((x, widths, family, minlen, window, gap, snr, perc))
    return cases


def run(compute_dir, cases, optimized):
    exe = compute_dir / ('cwt_custom_opt' if optimized else 'cwt_custom_checked')
    flags = ['-std=f2018', '-Werror=implicit-interface', '-O3'] if optimized else [
        '-std=f2018', '-Werror=implicit-interface', '-O0', '-g', '-fcheck=all']
    cmd = ['gfortran', *flags, '-J', str(compute_dir), '-I', str(compute_dir),
           str(ROOT/'src/kinds.f90'), str(ROOT/'src/signal_cwt_peaks.f90'),
           str(ROOT/'tools/signal_reference/cwt_custom_driver.f90'), '-o', str(exe)]
    subprocess.run(cmd, cwd=compute_dir, check=True, capture_output=True, text=True)
    inputs = []
    for x, widths, family, minlen, window, gap, snr, perc in cases:
        inputs.extend((f'{len(x)} {len(widths)} {family} {minlen} {window}',
                       ' '.join(f'{v:.17g}' for v in x),
                       ' '.join(f'{v:.17g}' for v in widths),
                       f'{gap:.17g} {snr:.17g} {perc:.17g}'))
    output = subprocess.run([str(exe)], input='\n'.join(inputs)+'\n', text=True,
                            capture_output=True, check=True, timeout=120).stdout.splitlines()
    index = 0
    max_error = 0.0
    peak_values = 0
    cwt_values = 0
    for case_id, (x, widths, family, minlen, window, gap, snr, perc) in enumerate(cases):
        status, npeak, transform_status = map(int, output[index].split())
        index += 1
        if status != 0 or transform_status != 0:
            raise AssertionError(f'case {case_id}: statuses {status}, {transform_status}')
        peaks = np.fromstring(output[index], dtype=int, sep=' ')
        index += 1
        if len(peaks) != npeak:
            raise AssertionError(f'case {case_id}: wrong number of peaks')
        callback = lambda points, width: wavelet(points, width, family)
        expected = np.asarray(find_peaks_cwt(x, widths, wavelet=callback,
            gap_thresh=gap, min_length=minlen, min_snr=snr,
            noise_perc=perc, window_size=window), dtype=int)+1
        if not np.array_equal(peaks, expected):
            raise AssertionError(
                f'case {case_id}, wavelet {family}: peaks={peaks} vs SciPy={expected}')
        peak_values += len(peaks)
        transform = np.empty((len(widths), len(x)))
        for j in range(len(widths)):
            transform[j] = np.fromstring(output[index], dtype=float, sep=' ')
            index += 1
        reference = _cwt(x, callback, widths)
        error = float(np.max(np.abs(transform-reference)))
        max_error = max(max_error, error)
        if not np.allclose(transform, reference, rtol=5e-12, atol=5e-12):
            raise AssertionError(f'case {case_id}: CWT maximum error {error:.3g}')
        cwt_values += transform.size
    if index != len(output):
        raise AssertionError(f'trailing Fortran output: {len(output)-index} lines')
    return max_error, cwt_values, peak_values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', type=int, default=600)
    args = parser.parse_args()
    cases = generate_cases(args.cases)
    with TemporaryDirectory(prefix='scifort_wavelet_') as tmp:
        work = Path(tmp)
        for opt in (False, True):
            err, coeffs, peaks = run(work, cases, opt)
            print(f"{'optimized' if opt else 'checked'}: {len(cases)}/{len(cases)} "
                  f'custom-wavelet SciPy matches, {coeffs} CWT samples, '
                  f'{peaks} peak indices, maximum absolute error {err:.3e}')


if __name__ == '__main__':
    main()
