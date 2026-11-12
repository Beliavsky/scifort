#!/usr/bin/env python3
"""Compare segment-wise Fortran detrending callbacks with SciPy's callable detrend."""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
from scipy import signal


def quadratic_detrend(samples: np.ndarray) -> np.ndarray:
    """Orthogonal quadratic least-squares removal along the last axis."""
    n = samples.shape[-1]
    t = (np.arange(1, n + 1, dtype=float) - 0.5 * (n + 1)) / n
    q = t * t - np.mean(t * t)
    av = np.mean(samples, axis=-1, keepdims=True)
    slope = np.sum(samples * t, axis=-1, keepdims=True) / np.sum(t * t)
    curve = np.sum(samples * q, axis=-1, keepdims=True) / np.sum(q * q)
    return samples - av - slope * t - curve * q


def expected_cases():
    t = np.arange(1, 65, dtype=float)
    x = np.sin(0.19 * t) + 0.25 * np.cos(0.72 * t) + 0.003 * t*t + 0.05*t
    y = np.cos(0.25 * t) + 0.14 * np.sin(0.34 * t) - 0.0004 * t*t + 0.025*t
    z = x + 0.45j*y
    w = y + 1j*(-0.2*x + 0.1*y)
    for ns, no in [(12, 6), (16, 8), (23, 9), (32, 13)]:
        nf = ns + 3
        onesided = ns in (12, 23)
        params = dict(fs=9.0, window='hann', nperseg=ns, noverlap=no,
                      nfft=nf, detrend=quadratic_detrend)
        # Return triples (frequency, time-grid or None, complex values).
        f, v = signal.periodogram(x, fs=9.0, window='hann', nfft=nf,
                                  detrend=quadratic_detrend)
        yield f, None, v
        f, v = signal.welch(x, **params, average='median', return_onesided=onesided)
        yield f, None, v
        f, v = signal.csd(x, y, **params, return_onesided=onesided)
        yield f, None, v
        f, v = signal.coherence(x, y, **params)
        yield f, None, v
        f, time, v = signal.stft(x, **params, boundary=None, padded=False,
                                 return_onesided=onesided)
        yield f, time, v
        f, time, v = signal.spectrogram(x, **params, mode='complex',
                                        return_onesided=onesided)
        yield f, time, v
        f, v = signal.periodogram(z, fs=9.0, window='hann', nfft=nf,
                                  detrend=quadratic_detrend)
        yield f, None, v
        f, v = signal.welch(z, **params, average='median')
        yield f, None, v
        f, v = signal.csd(z, w, **params)
        yield f, None, v
        f, v = signal.coherence(z, w, **params)
        yield f, None, v
        f, time, v = signal.stft(z, **params, boundary=None, padded=False)
        yield f, time, v
        f, time, v = signal.spectrogram(z, **params, mode='complex')
        yield f, time, v


def compare_file(path: Path) -> None:
    lines = path.read_text().splitlines()
    lines = iter(lines)
    max_error = 0.0
    values = 0
    count = 0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FutureWarning)
        warnings.filterwarnings('ignore', message='Input data is complex', category=UserWarning)
        for f, times, expected in expected_cases():
            header = next(lines).split()
            assert header[0] == 'CASE', header
            count += 1
            assert int(header[1]) == count
            rows, cols = map(int, header[2:])
            matrix = np.array([[float(v) for v in next(lines).split()]
                               for _ in range(rows*cols)]).reshape((cols, rows, 4))
            got_frequency = matrix[:, :, 0]
            got_times = matrix[:, :, 1]
            got = (matrix[:, :, 2] + 1j*matrix[:, :, 3]).T
            expected = np.asarray(expected).reshape(rows, cols)
            assert got.shape == expected.shape, (count, got.shape, expected.shape)
            assert np.max(np.abs(got_frequency - f[None, :])) < 2e-12, count
            if times is not None:
                assert np.max(np.abs(got_times - times[:, None])) < 2e-12, count
            elif np.max(np.abs(got_times)) != 0:
                raise AssertionError(f'Unexpected time values for case {count}')
            error = float(np.max(np.abs(got - expected)))
            if error > 5e-11:
                raise AssertionError(f'Case {count}: maximum error {error}')
            max_error = max(max_error, error)
            values += got.size
    assert count == 48
    print(f'PASS: {count} SciPy comparison cases, {values} values, '
          f'max absolute error {max_error:.3e}')


if __name__ == '__main__':
    compare_file(Path(sys.argv[1]))
