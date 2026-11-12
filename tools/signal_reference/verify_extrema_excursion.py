#!/usr/bin/env python3
"""Compare SciFort relative extrema against SciPy and excursion detector against Billauer algorithm."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
from scipy import signal


def excursion_reference(x: np.ndarray, delta: float):
    """Independent Billauer-style state machine; final unconfirmed swing is discarded."""
    mx, mn, imx, imn = [], [], [], []
    if not len(x):
        return imx, mx, imn, mn
    largest, smallest = float(x[0]), float(x[0])
    largest_i = smallest_i = 1
    looking_for_maximum = True
    for i, value in enumerate(x, start=1):
        if value > largest:
            largest, largest_i = float(value), i
        if value < smallest:
            smallest, smallest_i = float(value), i
        if looking_for_maximum:
            if value < largest - delta:
                imx.append(largest_i)
                mx.append(largest)
                smallest, smallest_i = float(value), i
                looking_for_maximum = False
        else:
            if value > smallest + delta:
                imn.append(smallest_i)
                mn.append(smallest)
                largest, largest_i = float(value), i
                looking_for_maximum = True
    return imx, mx, imn, mn


def cases():
    rng = np.random.default_rng(20261019)
    for j in range(450):
        n = int(rng.choice([0, 1, 2, 3, 7, 15, 41, 80]))
        x = rng.standard_normal(n)
        if j % 3 == 0:
            x = np.round(x, 1)
        if j % 7 == 0:
            x = np.round(x, 0)
        if j % 11 == 0:
            x = np.sin(np.arange(n) * 0.17) + 0.1 * x
        order = int(rng.choice([1, 2, 3, 4, 7, 20, 81]))
        mode = 'wrap' if j % 2 else 'clip'
        delta = float(rng.uniform(0.04, 2.5))
        yield x, order, mode, delta


def parse_array(line: str, typ):
    parts = line.split()
    count = int(parts[0])
    values = np.array([typ(z) for z in parts[1:]], dtype=typ)
    if len(values) != count:
        raise AssertionError(f'bad stream record: {line}')
    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('driver', type=Path)
    args = parser.parse_args()
    tests = list(cases())
    lines = []
    for x, order, mode, delta in tests:
        lines.append(f'{len(x)} {order} {int(mode == "wrap")} {delta:.17e}')
        if len(x):
            lines.append(' '.join(f'{float(t):.17e}' for t in x))
    run = subprocess.run([str(args.driver.resolve())], input='\n'.join(lines) + '\n',
                         text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    lines = iter(run.stdout.splitlines())
    for case_number, (x, order, mode, delta) in enumerate(tests):
        actual = [parse_array(next(lines), int) for _ in range(3)]
        expected = [np.asarray(signal.argrelmax(x, order=order, mode=mode)[0] + 1),
                    np.asarray(signal.argrelmin(x, order=order, mode=mode)[0] + 1),
                    np.asarray(signal.argrelextrema(x, np.greater_equal,
                                                   order=order, mode=mode)[0] + 1)]
        for label, got, want in zip(('max', 'min', 'greater_equal'), actual, expected):
            if not np.array_equal(got, want):
                raise AssertionError(f'case {case_number}: {label}: {got} != {want}')
        imx = parse_array(next(lines), int)
        mx = parse_array(next(lines), float)
        imn = parse_array(next(lines), int)
        mn = parse_array(next(lines), float)
        exp_imx, exp_mx, exp_imn, exp_mn = excursion_reference(x, delta)
        if not np.array_equal(imx, exp_imx) or not np.array_equal(imn, exp_imn):
            raise AssertionError(f'case {case_number}: excursion indices mismatch')
        if not np.allclose(mx, exp_mx, rtol=0, atol=1e-13) or not np.allclose(
                mn, exp_mn, rtol=0, atol=1e-13):
            raise AssertionError(f'case {case_number}: excursion values mismatch')
    try:
        next(lines)
    except StopIteration:
        pass
    else:
        raise AssertionError('unread stream output')
    print(f'{len(tests)} randomized cases: 1350 SciPy relative-extrema comparisons '
          f'and 450 Billauer-style detector comparisons passed')


if __name__ == '__main__':
    main()
