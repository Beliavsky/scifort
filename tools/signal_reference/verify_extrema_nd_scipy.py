#!/usr/bin/env python3
"""Check SciFort's 1..15-D relative-extrema masks against scipy.signal.

This compiles only scifort_kinds and the relative-extrema kernel in either a
checked or optimized configuration. No SciPy source code is copied.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]


def cases():
    rng = np.random.default_rng(20261027)
    supported_ranks = (1, 2, 3, 4, 5, 15)
    for j in range(720):
        r = supported_ranks[j % len(supported_ranks)]
        shape = [1] * r
        if r <= 5:
            shape = [int(rng.integers(1, 7)) for _ in range(r)]
            if np.prod(shape) > 800:
                shape = [int(rng.integers(1, 4)) for _ in range(r)]
        else:
            for d in (0, 3, 7, 11, 14):
                shape[d] = int(rng.integers(1, 4))
        if j % 41 == 0:
            shape[(j // 41) % r] = 0
        axis0 = int(rng.integers(0, r))
        axis = axis0 + 1 if j % 2 == 0 else axis0 - r
        order = [1, 2, 3, 5, 9][j % 5]
        mode = 'clip' if j % 2 == 0 else 'wrap'
        op = (j // 5) % 3 + 1
        x = rng.standard_normal(size=shape).astype(np.float64)
        if j % 3 == 0:
            x = np.round(x)
        if j % 19 == 0 and x.size:
            x.flat[int(rng.integers(x.size))] = np.nan
        cmp = (np.greater, np.less, np.greater_equal)[op - 1]
        location = signal.argrelextrema(x, cmp, axis=axis0, order=order, mode=mode)
        mask = np.zeros(x.shape, dtype=np.bool_)
        mask[location] = True
        yield x, axis, order, mode, op, mask


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--optimized', action='store_true')
    args = parser.parse_args()
    reference = list(cases())
    data = []
    total = 0
    for x, axis, order, mode, op, _ in reference:
        shape = list(x.shape) + [1] * (15 - x.ndim)
        data.append(' '.join(map(str, [x.ndim, *shape, axis, order, int(mode == 'wrap'), op])))
        if x.size:
            data.append(' '.join(f'{v:.17g}' for v in x.ravel(order='F')))
            total += x.size
    with tempfile.TemporaryDirectory(prefix='scifort-extrema-nd-') as temp:
        temp = Path(temp)
        cmd = ['gfortran', '-std=f2018', '-Werror=implicit-interface',
               '-I', str(temp), '-J', str(temp)]
        cmd += ['-O3'] if args.optimized else ['-O0', '-g', '-fcheck=all']
        cmd += ['-o', str(temp/'test_extrema_nd'),
                str(ROOT/'src/kinds.f90'),
                str(ROOT/'src/signal_relative_extrema.f90'),
                str(ROOT/'tools/signal_reference/extrema_nd_driver.f90')]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        run = subprocess.run([str(temp/'test_extrema_nd')], input='\n'.join(data)+'\n',
                             capture_output=True, text=True, check=True, timeout=90)
        lines = iter(run.stdout.splitlines())
        for i, (_, _, _, _, _, ref) in enumerate(reference):
            code, count = map(int, next(lines).split())
            if code != 0:
                raise AssertionError(f'case {i}: returned error {code}')
            got = np.array([bool(int(v)) for v in next(lines).split()] if ref.size else [],
                           dtype=np.bool_).reshape(ref.shape, order='F')
            if not np.array_equal(got, ref) or count != ref.sum():
                raise AssertionError(f'case {i}: rank={ref.ndim}, ' +
                                     f'expected count={ref.sum()}, got={count}, ' +
                                     f'mismatch count={np.count_nonzero(got!=ref)}')
        if next(lines, None) is not None:
            raise AssertionError('driver emitted additional output')
    print(f'PASS {len(reference)} independent SciPy relative-extrema comparisons; '
          f'{total} positions checked; optimized={args.optimized}')


if __name__ == '__main__':
    main()
