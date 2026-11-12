#!/usr/bin/env python3
"""Compare rank-two/three SciFort SOS and forward-backward filters to SciPy."""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261018)


def encode(array):
    return [f'{v.real:.17g} {v.imag:.17g}' for v in np.ravel(array, order='F')]


def generate():
    cases, lines = [], []
    for op in (1, 2, 3):
        for rank in (2, 3):
            for kind in (0, 1):
                for case in range(40):
                    dims = list(map(int, RNG.integers(12, 21, size=rank)))
                    axis = case % rank
                    full_shape = dims + [1] * (3 - rank)
                    ns = 1 + case % 3
                    haszi = case % 2
                    mode_id = 1 + case % 4
                    padlen = (0, 2, 5, 8)[case % 4]
                    sos = np.zeros((ns, 6), dtype=complex)
                    sos[:, :3] = RNG.normal(scale=0.18, size=(ns, 3))
                    sos[:, 3] = 1
                    sos[:, 4] = RNG.normal(scale=0.14, size=ns)
                    sos[:, 5] = RNG.normal(scale=0.04, size=ns)
                    b = RNG.normal(scale=0.23, size=3).astype(complex)
                    a = np.array([1.0, -0.25, 0.05], dtype=complex)
                    x = RNG.normal(size=dims)
                    if kind:
                        sos[:, :3] += 0.03j * RNG.normal(size=(ns, 3))
                        sos[:, 4:6] += 0.015j * RNG.normal(size=(ns, 2))
                        b += 0.06j * RNG.normal(size=3)
                        a[1:] += 0.02j * RNG.normal(size=2)
                        x = x + 0.1j * RNG.normal(size=dims)
                    if not kind:
                        sos = np.ascontiguousarray(sos.real)
                        b, a = b.real, a.real
                    state_shape = list(dims)
                    state_shape[axis] = 2
                    zi = np.zeros((ns, *state_shape), dtype=np.result_type(sos, x))
                    if haszi:
                        zi = 0.09 * RNG.normal(size=(ns, *state_shape))
                        if kind:
                            zi = zi.astype(complex) + 0.06j * RNG.normal(size=(ns, *state_shape))
                    modes = ('odd', 'even', 'constant', None)
                    if op == 1:
                        expected, zf = signal.sosfilt(sos, x, axis=axis, zi=zi)
                    elif op == 2:
                        expected = signal.filtfilt(b, a, x, axis=axis,
                            padtype=modes[mode_id - 1], padlen=padlen)
                        zf = None
                    else:
                        expected = signal.sosfiltfilt(sos, x, axis=axis,
                            padtype=modes[mode_id - 1], padlen=padlen)
                        zf = None
                    lines.append(' '.join(map(str, (
                        op, kind, rank, *full_shape, axis + 1, ns,
                        haszi, mode_id, padlen))))
                    lines.extend(encode(sos))
                    lines.extend(encode(b))
                    lines.extend(encode(a))
                    lines.extend(encode(x))
                    if haszi:
                        lines.extend(encode(zi))
                    cases.append((expected, zf, rank, full_shape))
    return cases, lines


def main():
    cases, lines = generate()
    with tempfile.TemporaryDirectory(prefix='scifort-filter-axes-') as temp:
        work = Path(temp)
        cmd = ['gfortran', '-std=f2018', '-O2', '-Werror=implicit-interface',
               '-I', str(work), '-J', str(work)]
        cmd += [str(ROOT / 'src' / (name + '.f90')) for name in
                ('kinds', 'signal_lfilter', 'signal_sos', 'signal_filtfilt')]
        cmd += [str(ROOT / 'tools/signal_reference/filter_axes_driver.f90'),
                '-o', str(work / 'driver')]
        subprocess.run(cmd, cwd=work, check=True)
        result = subprocess.run([str(work / 'driver')],
            input=f'{len(cases)}\n' + '\n'.join(lines) + '\n',
            text=True, capture_output=True, timeout=180, check=False)
        if result.returncode:
            raise RuntimeError(f'Fortran driver exit {result.returncode}: {result.stderr[-2500:]}')
        output = result.stdout
    it = iter(output.splitlines())
    worst = 0.0
    for num, (expected, zf, rank, full_shape) in enumerate(cases):
        values = list(map(int, next(it).split()))
        status, *shape = values
        if status != 0 or shape[:3] != full_shape:
            raise AssertionError(f'case {num}: status={status}, dimensions={shape}')
        count = int(np.prod(full_shape))
        actual = np.array([complex(*map(float, next(it).split()))
            for _ in range(count)]).reshape(full_shape, order='F')
        if rank == 2:
            actual = actual[:, :, 0]
        difference = float(np.max(np.abs(actual - expected)))
        worst = max(worst, difference)
        if not np.allclose(actual, expected, rtol=4e-12, atol=4e-12):
            raise AssertionError(f'case {num}: output diff {difference}')
        if zf is None:
            if shape[3:] != [0, 0, 0, 0]:
                raise AssertionError(f'case {num}: unexpected filter state shape {shape[3:]}')
            continue
        state_shape = (zf.shape + (1,)) if rank == 2 else zf.shape
        if shape[3:] != list(state_shape):
            raise AssertionError(f'case {num}: final state shape {shape[3:]} vs {state_shape}')
        count = int(np.prod(state_shape))
        final_state = np.array([complex(*map(float, next(it).split()))
            for _ in range(count)]).reshape(state_shape, order='F')
        if rank == 2:
            final_state = final_state[:, :, :, 0]
        diff_state = float(np.max(np.abs(final_state - zf)))
        worst = max(worst, diff_state)
        if not np.allclose(final_state, zf, atol=4e-12, rtol=4e-12):
            raise AssertionError(f'case {num}: final-state difference {diff_state}')
    print(f'PASS: {len(cases)} multidimensional filter cases; maximum absolute error {worst:.3g}')


if __name__ == '__main__':
    main()
