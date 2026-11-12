#!/usr/bin/env python3
"""Independent SciPy rank-15 comparisons for convolve_nd/correlate_nd."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261025)
CHOICES = (0, 1, 2, 4, 9, 14)


def encode(values):
    return [f"({z.real:.17e},{z.imag:.17e})" for z in np.asarray(values).ravel(order="F")]


def build_and_check(args):
    lines = []
    cases = []
    for index in range(args.cases):
        kind = index % 4  # real/real, complex/complex, real/complex, complex/real
        operation = (index // 4) % 2  # convolve vs correlate
        method = ('direct', 'fft', 'auto')[(index // 8) % 3]
        mode = ('full', 'same', 'valid')[(index // 24) % 3]
        all_axes = index % 13 == 0
        selected = tuple(range(15)) if all_axes else tuple(sorted(
            RNG.choice(CHOICES, size=int(RNG.integers(1, 4)), replace=False)))
        n = [1] * 15
        m = [1] * 15
        reversed_order = index % 2 == 0
        for axis in selected:
            if all_axes and axis not in CHOICES:
                continue
            a, b = map(int, RNG.integers(1, 5, size=2))
            if mode == 'valid':
                a, b = sorted((a, b), reverse=not reversed_order)
            n[axis], m[axis] = a, b
        if not all_axes:
            for axis in RNG.choice([i for i in CHOICES if i not in selected],
                                   size=min(2, len(CHOICES) - len(selected)), replace=False):
                extent = int(RNG.integers(1, 4))
                n[axis], m[axis] = (extent, 1) if index % 3 else (1, extent)
        if np.prod(np.maximum(n, m)) > 450:
            # Keep a few nontrivial dimensions, including high axes, bounded.
            for axis in CHOICES[2:-1]:
                if np.prod(np.maximum(n, m)) <= 450:
                    break
                n[axis] = m[axis] = 1
        # Some ranks have only singleton axes, but the actual Fortran rank is 15.
        a = RNG.normal(size=n)
        b = RNG.normal(size=m)
        if kind in (1, 3):
            a = a + 1j * RNG.normal(size=n)
        if kind in (1, 2):
            b = b + 1j * RNG.normal(size=m)
        axes = None if all_axes else selected
        expected = None
        if all_axes:
            if operation:
                expected = signal.correlate(a, b, mode=mode, method=method)
            else:
                expected = signal.convolve(a, b, mode=mode, method=method)
        else:
            transformed = b if operation == 0 else np.conj(np.flip(b, axis=selected))
            expected = signal.fftconvolve(a, transformed, mode=mode, axes=selected)
        lines.append(f"{kind} {operation} {['direct','fft','auto'].index(method)+1} "
                     f"{['full','same','valid'].index(mode)+1} {len(selected) if axes is not None else -1}")
        lines.append(' '.join(map(str, n)))
        lines.append(' '.join(map(str, m)))
        if axes is not None:
            lines.append(' '.join(str(k+1) for k in axes))
        lines += encode(a)
        lines += encode(b)
        cases.append((index, kind, operation, method, mode, selected, a.shape, b.shape, expected))

    with tempfile.TemporaryDirectory(prefix='signal_convolve_nd_') as tmp:
        build = Path(tmp)
        names = ['kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral',
                 'signal_fftconvolve', 'signal_convolve_nd']
        flags = ['-std=f2018', '-Werror=implicit-interface', '-J'+str(build),
                 '-I'+str(build)]
        if args.checked:
            flags += ['-O0', '-g', '-fcheck=all', '-fbacktrace']
        else:
            flags += ['-O3']
        cmd = ['gfortran', *flags, *(str(ROOT/'src'/f'{name}.f90') for name in names),
               str(ROOT/'tools/signal_reference/convolve_nd_driver.f90'), '-o',
               str(build/'reference_driver')]
        subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=build)
        response = subprocess.run([str(build/'reference_driver')], input='\n'.join(lines)+'\n',
                                  capture_output=True, text=True, check=True, cwd=build,
                                  timeout=120)
        stream = iter(response.stdout.splitlines())
        worst = 0.0
        evaluated = 0
        for index, kind, operation, method, mode, selected, na, nb, expected in cases:
            descriptor = list(map(int, next(stream).split()))
            status, dims = descriptor[0], tuple(descriptor[1:])
            if status != 0 or dims != expected.shape:
                raise AssertionError(f'case {index}: status={status}, output={dims}, '
                                     f'expected={expected.shape}, kind={kind}, '
                                     f'mode={mode}, selected={selected}, shapes={na} {nb}')
            observed = np.array([complex(*map(float, next(stream).split()))
                                 for _ in range(expected.size)]).reshape(dims, order='F')
            diff = np.max(np.abs(observed - expected))
            worst = max(worst, float(diff))
            evaluated += observed.size
            if not np.allclose(observed, expected, atol=2.5e-11, rtol=2e-10):
                raise AssertionError(f'case {index} op={operation} method={method} '
                                     f'mode={mode} axes={selected} n={na} m={nb} '
                                     f'max error={diff}')
        print(f'PASS {len(cases)} SciPy convolve/correlate cases; '
              f'{evaluated} outputs; max abs error {worst:.3e}; '
              f'{"checked" if args.checked else "optimized"} build')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', type=int, default=1008)
    parser.add_argument('--checked', action='store_true')
    build_and_check(parser.parse_args())
