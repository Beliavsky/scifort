#!/usr/bin/env python3
"""Seeded SciPy 1.17 comparison for rank-5/rank-15 lfilter_nd."""
from __future__ import annotations

from pathlib import Path
import argparse
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
SHAPES = {5: (4, 3, 2, 3, 5), 15: (3, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4)}


def encode(z):
    for value in np.ravel(z, order='F'):
        yield f'({value.real:.17g},{value.imag:.17g})'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimized', action='store_true')
    args = parser.parse_args()
    rng = np.random.default_rng(20261018)
    cases = []
    lines = []
    for rank_id, repetitions in ((5, 20), (15, 10)):
        shape = SHAPES[rank_id]
        for combo in range(8):
            for j in range(repetitions):
                axis = (j + combo) % rank_id
                with_state = int(j % 2 == 0)
                b_complex = bool(combo & 1)
                a_complex = bool(combo & 2)
                x_complex = bool(combo & 4)
                raw_b = rng.normal(size=3) * 0.3 + 1j * rng.normal(size=3) * 0.1
                raw_a = np.array([1.3 + rng.random(), -0.25 + 0.1*rng.random()], dtype='complex128')
                raw_a += 1j * np.array([0.07, -0.04])
                raw_x = rng.normal(size=shape) + 0.4j*rng.normal(size=shape)
                state_shape = list(shape)
                state_shape[axis] = 2
                raw_zi = 0.1*rng.normal(size=state_shape) + 0.03j*rng.normal(size=state_shape)
                b = raw_b if b_complex else raw_b.real
                a = raw_a if a_complex else raw_a.real
                x = raw_x if x_complex else raw_x.real
                complex_output = b_complex or a_complex or x_complex
                zi = raw_zi if complex_output else raw_zi.real
                if not with_state:
                    zi = np.zeros(state_shape, dtype='complex128' if complex_output else 'float64')
                y, zf = signal.lfilter(b, a, x, zi=zi, axis=axis)
                cases.append((y, zf))
                lines.append(f'{rank_id} {axis+1} {combo} {with_state}')
                lines.extend(encode(raw_b))
                lines.extend(encode(raw_a))
                lines.extend(encode(raw_x))
                if with_state:
                    lines.extend(encode(raw_zi))

    with tempfile.TemporaryDirectory(prefix='scifort-lfsci-') as td:
        tmp = Path(td)
        flags = ['-std=f2018', '-Wimplicit-interface', '-Werror=implicit-interface']
        flags += ['-O3'] if args.optimized else ['-O0', '-g', '-fcheck=all']
        command = ['gfortran', *flags, f'-J{tmp}', f'-I{tmp}',
                   str(ROOT/'src/kinds.f90'), str(ROOT/'src/signal_lfilter.f90'),
                   str(ROOT/'tools/signal_reference/lfilter_nd_driver.f90'),
                   '-o', str(tmp/'probe')]
        subprocess.run(command, cwd=tmp, check=True, capture_output=True, text=True)
        result = subprocess.run([str(tmp/'probe')], input=str(len(cases))+'\n'+'\n'.join(lines)+'\n',
                                text=True, capture_output=True, check=True, cwd=tmp, timeout=150)
        out = iter(result.stdout.splitlines())
        max_diff = 0.0
        for number, (expected_y, expected_zf) in enumerate(cases):
            code = int(next(out))
            if code:
                raise AssertionError(f'case {number}: Fortran status {code}')
            def read_array(shape):
                data = [complex(*map(float, next(out).split())) for _ in range(int(np.prod(shape)))]
                return np.asarray(data, dtype=np.complex128).reshape(shape, order='F')
            got_y = read_array(expected_y.shape)
            got_zf = read_array(expected_zf.shape)
            delta = max(float(np.max(np.abs(got_y-expected_y))),
                        float(np.max(np.abs(got_zf-expected_zf))))
            max_diff = max(max_diff, delta)
            if not np.allclose(got_y, expected_y, rtol=3e-12, atol=3e-12):
                raise AssertionError(f'case {number}: output differs: {delta}')
            if not np.allclose(got_zf, expected_zf, rtol=3e-12, atol=3e-12):
                raise AssertionError(f'case {number}: states differ: {delta}')
        if next(out, None) is not None:
            raise AssertionError('Unexpected extra output lines')
        print(f'PASS: {len(cases)} rank-5/rank-15 lfilter_nd SciPy cases; '
              f'max absolute difference {max_diff:.3g}; '
              f'{"-O3" if args.optimized else "runtime checked"}')


if __name__ == '__main__':
    main()
