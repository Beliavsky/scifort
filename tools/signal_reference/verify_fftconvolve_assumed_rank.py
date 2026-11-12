#!/usr/bin/env python3
"""Seeded SciPy 1.17 comparisons for the rank-15 caller-output FFT API."""

from __future__ import annotations

import pathlib
import subprocess
import tempfile
import numpy as np
from scipy import signal

ROOT = pathlib.Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261019)


def encode_array(values: np.ndarray) -> list[str]:
    values = values.ravel(order='F')
    return [f'{float(v.real):.17e} {float(v.imag):.17e}' for v in values]


def main() -> None:
    lines: list[str] = []
    cases = []
    choices = [0, 1, 2, 3, 4, 14]
    for case in range(576):
        kind = case % 4
        mode = ('full', 'same', 'valid')[(case // 4) % 3]
        n = [1] * 15
        m = [1] * 15
        nconv = 1 + (case % 3)
        axes = tuple(sorted(RNG.choice(choices, size=nconv, replace=False).tolist()))
        # Pick a common length ordering for valid-mode comparisons.
        reverse = bool(case % 2)
        for axis in axes:
            length_a, length_b = map(int, RNG.integers(1, 4, size=2))
            if mode == 'valid':
                length_a, length_b = sorted((length_a,length_b),reverse=not reverse)
            n[axis], m[axis] = length_a, length_b
        for axis in choices:
            if axis in axes:
                continue
            length = int(RNG.integers(1, 4))
            if case % 2:
                n[axis], m[axis] = length, 1
            else:
                n[axis], m[axis] = 1, length
        # Bound the total number of array values, but keep strided axes.
        if np.prod(np.maximum(n,m)) > 400:
            for axis in choices[2:]:
                if axis not in axes:
                    n[axis],m[axis] = 1,1
        a = RNG.normal(size=n)
        b = RNG.normal(size=m)
        if kind in (1,3):
            a = a + 1j * RNG.normal(size=n)
        if kind in (1,2):
            b = b + 1j * RNG.normal(size=m)
        mode_index = ('full', 'same', 'valid').index(mode) + 1
        mask = sum(1 << ax for ax in axes)
        lines.append(f'{kind} {mode_index} {mask}')
        lines.append(' '.join(map(str,n)))
        lines.append(' '.join(map(str,m)))
        lines.extend(encode_array(a))
        lines.extend(encode_array(b))
        expected = signal.fftconvolve(a,b,mode=mode,axes=axes)
        cases.append((case,kind,mode,axes,n,m,expected))
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        cmd = ['gfortran','-std=f2018','-O2','-Werror=implicit-interface',
               f'-J{work}',f'-I{work}']
        cmd += [str(ROOT / 'src' / f'{name}.f90') for name in
                ('kinds','signal_dpss', 'signal_windows_extra', 'signal_spectral','signal_convolve','signal_fftconvolve')]
        cmd += [str(ROOT/'tools/signal_reference/fftconvolve_assumed_rank_driver.f90'),
                '-o',str(work/'driver')]
        subprocess.run(cmd,cwd=work,check=True)
        out = subprocess.run([str(work/'driver')],input='\n'.join(lines)+'\n',
                             text=True,capture_output=True,check=True,timeout=120)
        stream = iter(out.stdout.splitlines())
        worst = 0.0
        for case,kind,mode,axes,n,m,expected in cases:
            info = list(map(int,next(stream).split()))
            status, shape = info[0], tuple(info[1:])
            if status != 0 or shape != expected.shape:
                raise AssertionError(f'case {case} status={status} shape={shape} '
                                     f'expected={expected.shape}')
            obtained = np.array([complex(*map(float,next(stream).split()))
                                 for _ in range(expected.size)]).reshape(shape,order='F')
            diff = np.max(np.abs(obtained-expected))
            worst = max(worst,diff)
            if not np.allclose(obtained,expected,atol=3e-11,rtol=1e-10):
                raise AssertionError(f'case {case}: kind={kind} mode={mode} '
                                     f'axes={axes} n={n} m={m} error={diff}')
        print(f'PASS: {len(cases)} rank-15 SciPy fftconvolve comparisons; '
              f'maximum absolute error {worst:.3e}')


if __name__ == '__main__':
    main()
