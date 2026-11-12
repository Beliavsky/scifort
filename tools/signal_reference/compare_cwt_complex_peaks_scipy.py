#!/usr/bin/env python3
"""Compare SciFort complex CWT magnitude ridges with an independent SciPy reference.

The extension itself is SciFort-specific: we apply SciPy's existing real-valued
ridge finder to the *magnitude* of a complex transform calculated by SciPy's
public convolution routine. No SciPy complex find_peaks_cwt API is claimed.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import subprocess

import numpy as np
from scipy.signal import convolve
from scipy.signal._peak_finding import _identify_ridge_lines, _filter_ridge_lines


def wavelet(points: float, width: float, omega0: float, custom: bool) -> np.ndarray:
    k = np.arange(int(np.ceil(points)), dtype=float)
    r = (k - (points - 1.0) / 2.0) / width
    w = (np.pi ** -0.25 / np.sqrt(width)) * np.exp(-0.5 * r * r) * (
        np.cos(omega0 * r) + 1j * np.sin(omega0 * r)
    )
    if custom:
        w *= (1.0 + 0.08 * r + 0.2j)
    return w


def reference(x: np.ndarray, widths: np.ndarray, omega0: float, custom: bool,
              distances: np.ndarray, gap: float, minlen: int, minsnr: float,
              perc: float, win: int) -> np.ndarray:
    rows = []
    for width in widths:
        points = min(10.0 * width, x.size)
        sampled = wavelet(points, width, 5.0 if custom else omega0, custom)
        rows.append(convolve(x, np.conj(sampled[::-1]), mode='same', method='direct'))
    magnitude = np.abs(np.asarray(rows))
    ridges = _identify_ridge_lines(magnitude, distances, gap)
    good = _filter_ridge_lines(magnitude, ridges, window_size=win, min_length=minlen,
                               min_snr=minsnr, noise_perc=perc)
    return np.sort(np.asarray([r[1][0] + 1 for r in good], dtype=int))


def create_cases():
    rng = np.random.default_rng(4912026)
    cases = []
    # Reproducible seeded signals, varied scale grids, both input types and
    # asymmetric callbacks; avoid artificially flat coefficient maps.
    for k in range(280):
        n = int(rng.integers(35, 106))
        nw = int(rng.integers(3, 11))
        mode = k % 4
        t = np.arange(n, dtype=float)
        x = (rng.normal(0.0, .07, n) + 1j * rng.normal(0.0, .07, n))
        for _ in range(3):
            center = rng.uniform(8, n - 8)
            scale = rng.uniform(2.2, 8.0)
            amp = rng.uniform(.4, 3.2)
            freq = rng.uniform(.4, 2.0)
            phase = rng.uniform(-np.pi, np.pi)
            x += amp * np.exp(-.5 * ((t - center) / scale)**2) * np.exp(1j*(freq*t+phase))
        if mode % 2 == 0:
            x = np.real(x).astype(complex)
        widths = np.sort(rng.uniform(1.0, 10.0, nw))
        if k % 5 == 0:
            widths = np.arange(1, nw+1, dtype=float)
        gap = float(rng.uniform(0.8, 3.6))
        snr = float(rng.uniform(.5, 2.0))
        perc = float(rng.uniform(7.0, 55.0))
        minlen = int(rng.integers(1, nw + 1))
        win = int(rng.choice([5, 7, 9, 11, 13]))
        omega0 = float(rng.uniform(2.0, 6.5))
        use_dist = k % 3 != 0
        distances = rng.uniform(.2, 4.5, nw) if use_dist else widths / 4
        cases.append((x, widths, mode, gap, snr, perc, minlen, win,
                      omega0, use_dist, distances))
    return cases


def run(driver: Path) -> None:
    lines = []
    expected = []
    counts = 0
    for x, widths, mode, gap, snr, perc, minlen, win, omega0, use_dist, distances in create_cases():
        lines.append(f'{x.size} {widths.size} {mode} {minlen} {win} {int(use_dist)}')
        lines.append(' '.join(f'({v.real:.17g},{v.imag:.17g})' for v in x))
        lines.append(' '.join(f'{v:.17g}' for v in widths))
        lines.append(f'{gap:.17g} {snr:.17g} {perc:.17g} {omega0:.17g}')
        if use_dist:
            lines.append(' '.join(f'{v:.17g}' for v in distances))
        expected.append(reference(x, widths, omega0, mode >= 2,
                                  distances, gap, minlen, snr, perc, win))
        counts += int(widths.size * x.size)
    result = subprocess.run([str(driver)], input='\n'.join(lines)+'\n',
                            text=True, capture_output=True, check=True)
    outputs = result.stdout.splitlines()
    if len(outputs) != len(expected):
        raise AssertionError(f'Expected {len(expected)} outputs; got {len(outputs)}: '
                             f'{result.stderr[-500:]}')
    npeaks = 0
    for k, (line, want) in enumerate(zip(outputs, expected)):
        got = np.fromstring(line, dtype=int, sep=' ')
        if len(got) < 2 or got[0] != 0 or got[1] != want.size or not np.array_equal(got[2:], want):
            raise AssertionError(f'Case {k}: Fortran={got}, SciPy magnitude-ridges={want}')
        npeaks += want.size
    print(f'PASS {len(expected)} SciPy magnitude-ridge cases, '
          f'{counts} CWT positions checked, {npeaks} accepted peaks')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver', type=Path)
    args = parser.parse_args()
    run(args.driver.resolve())


if __name__ == '__main__':
    main()
