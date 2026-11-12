#!/usr/bin/env python3
"""Independently compare complex Wiener filtering against scipy.signal.wiener."""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from scipy.signal import wiener


def inputs() -> tuple[np.ndarray, ...]:
    x1 = np.empty(13, dtype=np.complex128)
    x2 = np.empty((8, 7), dtype=np.complex128)
    x3 = np.empty((5, 4, 3), dtype=np.complex128)
    for i in range(1, 14):
        x1[i - 1] = ((17 * i % 29 - 14) / 7 +
                      1j * ((19 * i + 3) % 23 - 11) / 9)
    for j in range(1, 8):
        for i in range(1, 9):
            x2[i - 1, j - 1] = ((7 * i + 11 * j) % 37 - 18) / 8 + \
                1j * ((13 * i + 26 * j) % 31 - 15) / 10
    for k in range(1, 4):
        for j in range(1, 5):
            for i in range(1, 6):
                x3[i - 1, j - 1, k - 1] = ((11 * i + 7 * j + 5 * k) % 41 - 20) / 11 + \
                    1j * ((17 * i + 13 * j + 19 * k) % 43 - 21) / 13
    return x1, x2, x3


def expected_cases() -> dict[int, np.ndarray]:
    out = {}
    shapes = (
        [None, 2, 5, 3, (4,), (1,), 7, 4, 9],
        [None, 2, (3, 2), (4, 3), (1, 5), (7, 1), (2, 4), 5, (5, 5)],
        [None, 2, (3, 2, 5), (4, 3, 2), (1, 5, 2), (5, 1, 3), (2, 4, 1), 3, (3, 3, 3)],
    )
    noise = (
        [None, None, 0.12, 0.3-0.14j, -0.2+0.13j, None, 0.0, 0.0j, 0.65],
        [None, None, 0.12, 0.3-0.14j, None, -0.2+0.13j, 0.0, 0.0j, 0.65],
        [None, None, 0.12, 0.3-0.14j, None, -0.2+0.13j, 0.0, 0.0j, 0.65],
    )
    for rank, x in enumerate(inputs()):
        for index, (shape, n) in enumerate(zip(shapes[rank], noise[rank]), start=1):
            with np.errstate(invalid='ignore', divide='ignore'):
                result = wiener(x, mysize=shape, noise=n)
            out[100 * rank + index] = np.asarray(result).ravel(order='F')
    return out


def check(path: Path) -> None:
    expected = expected_cases()
    lines = path.read_text().splitlines()
    offset = 0
    total = 0
    max_error = 0.0
    while offset < len(lines):
        tokens = lines[offset].split()
        assert len(tokens) == 3 and tokens[0] == 'CASE', tokens
        key, n = int(tokens[1]), int(tokens[2])
        assert key in expected, key
        actual = np.array([complex(*map(float, line.split()))
                           for line in lines[offset + 1:offset + 1 + n]])
        goal = expected.pop(key)
        assert n == goal.size, (key, n, goal.size)
        assert np.array_equal(np.isnan(actual.real), np.isnan(goal.real)), f'Real NaN mask differs: case {key}'
        assert np.array_equal(np.isnan(actual.imag), np.isnan(goal.imag)), f'Imag NaN mask differs: case {key}'
        error = np.abs(actual - goal)
        finite = np.isfinite(error)
        worst = float(np.max(error[finite])) if np.any(finite) else 0.0
        if worst > 3e-12:
            at = np.where(finite)[0][np.argmax(error[finite])]
            raise AssertionError(f'case {key} sample {at}: {actual[at]} != {goal[at]}')
        max_error = max(max_error, worst)
        total += n
        offset += n + 1
    assert not expected, f'Missing comparisons: {list(expected)}'
    print(f'PASS: 27 complex Wiener SciPy comparisons, {total} output samples; max abs error {max_error:.3e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    check(parser.parse_args().results)
