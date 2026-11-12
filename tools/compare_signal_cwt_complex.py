#!/usr/bin/env python3
"""Compare complex CWT and Morlet coefficients with SciPy convolution."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.signal import convolve


def morlet(length: float, width: float, omega: float) -> np.ndarray:
    position = (np.arange(int(np.ceil(length))) - (length - 1.0) / 2.0) / width
    return (np.pi ** (-0.25) / np.sqrt(width)) * np.exp(-0.5 * position**2) * np.exp(1j * omega * position)


def asymmetric(length: float, width: float) -> np.ndarray:
    position = (np.arange(int(np.ceil(length))) - (length - 1.0) / 2.0) / width
    return np.exp(-0.5 * position**2) * (1.0 + 0.1 * position + 1j * (0.7 * np.sin(position) + 0.2))


def run(path: Path) -> None:
    lines = iter(path.read_text(encoding="utf-8").splitlines())
    errors = []
    values = 0
    total = 0
    for header in lines:
        parts = header.split()
        if parts[0] == 'TOTAL':
            assert int(parts[1]) == total, (parts, total)
            break
        assert parts[0] == 'CASE', header
        case, kind, mode, n, nw = map(int, parts[1:6])
        omega = float(parts[6])
        assert case == total + 1
        x = np.array([complex(*map(float, next(lines).split())) for _ in range(n)])
        widths = np.array([float(next(lines)) for _ in range(nw)])
        got = np.array([complex(*map(float, next(lines).split())) for _ in range(nw * n)]).reshape(nw, n)
        if kind == 1:
            assert np.array_equal(x.imag, np.zeros(n))
            x = x.real
        expected = np.empty((nw, n), dtype=np.complex128)
        for i, width in enumerate(widths):
            length = min(10 * width, n)
            wave = (morlet(length, width, omega if mode == 1 else 5.0)
                    if mode in (1, 2) else asymmetric(length, width))
            expected[i] = convolve(x, np.conj(wave[::-1]), mode='same', method='direct')
        errors.append(float(np.max(np.abs(got - expected))))
        values += got.size
        total += 1
    assert total == 240, total
    worst = max(errors)
    assert worst < 5e-13, f'largest CWT error {worst:.4g} exceeds tolerance'
    print(f'PASS {total} SciPy complex CWT cases, {values} coefficients, max error {worst:.3g}')


if __name__ == '__main__':
    run(Path(sys.argv[1]))
