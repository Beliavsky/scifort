#!/usr/bin/env python3
"""Compare SciFort chirp/sweep_poly/gausspulse with installed SciPy.

Run from any directory: python tools/signal_reference/verify_waveforms.py
Needs NumPy, SciPy and gfortran. Does not introduce a package dependency.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261016)


def numbers(x: np.ndarray) -> str:
    return " ".join(f"{float(v):.17e}" for v in np.ravel(x))


def main() -> None:
    lines: list[str] = []
    reference: list[np.ndarray | float] = []
    mode_list: list[int] = []
    for idx in range(620):
        method_id = idx % 4 + 1
        method = ['linear', 'quadratic', 'logarithmic', 'hyperbolic'][method_id-1]
        f0 = float(RNG.uniform(.15, 17)) * (-1 if RNG.random() < .3 else 1)
        f1 = float(RNG.uniform(.15, 17)) * (-1 if RNG.random() < .2 else 1)
        if method == 'logarithmic':
            f1 = np.copysign(abs(f1),f0)
        if RNG.random() < .08:
            f1 = f0
        t1 = float(RNG.uniform(.4, 5))
        phi = float(RNG.uniform(-200, 200))
        vertex_zero = bool(RNG.integers(0, 2))
        n = int(RNG.integers(0, 48))
        t = RNG.uniform(-.3*t1, 1.2*t1, n)
        if n > 0 and RNG.random() < .2:
            t[0] = 0.0
        mode = 1 if idx < 360 else 2
        expected = signal.chirp(t, f0, t1, f1, method=method,
                                phi=phi, vertex_zero=vertex_zero, complex=mode == 2)
        lines.extend([str(mode),
                      f"{n} {f0:.17e} {t1:.17e} {f1:.17e} {phi:.17e} "
                      f"{int(vertex_zero)} {method_id}"])
        if n:
            lines.append(numbers(t))
        reference.append(expected)
        mode_list.append(mode)

    for idx in range(250):
        n = int(RNG.integers(0, 48))
        degree = int(RNG.integers(0, 9))
        phi = float(RNG.uniform(-180, 180))
        t = RNG.uniform(-1.5, 2.0, size=n)
        poly = RNG.uniform(-.6, .6, size=degree)
        expected = signal.sweep_poly(t, poly, phi=phi)
        lines.extend(['3',f'{n} {degree} {phi:.17e}'])
        if degree:
            lines.append(numbers(poly))
        if n:
            lines.append(numbers(t))
        reference.append(expected)
        mode_list.append(3)

    for idx in range(240):
        n = int(RNG.integers(0, 50))
        fc = float(RNG.uniform(.1, 30))
        if idx % 30 == 0:
            fc = 0.0
        bw = float(RNG.uniform(.05, 2))
        bwr = float(RNG.uniform(-40, -.1))
        t = RNG.uniform(-1.0, 1.0, size=n)
        expected = np.array(signal.gausspulse(t, fc=fc,bw=bw,bwr=bwr,
                                               retquad=True,retenv=True)).T
        lines.extend(['4',f'{n} {fc:.17e} {bw:.17e} {bwr:.17e}'])
        if n:
            lines.append(numbers(t))
        reference.append(expected)
        mode_list.append(4)

    for idx in range(100):
        fc = float(RNG.uniform(.05, 30))
        bw = float(RNG.uniform(.05, 2))
        bwr = float(RNG.uniform(-40, -.1))
        tpr = float(RNG.uniform(-130, -5))
        expected = float(signal.gausspulse('cutoff', fc=fc,bw=bw,bwr=bwr,tpr=tpr))
        lines.extend(['5',f'{fc:.17e} {bw:.17e} {bwr:.17e} {tpr:.17e}'])
        reference.append(expected)
        mode_list.append(5)

    with tempfile.TemporaryDirectory(prefix='scifort-waveforms-') as temp:
        work = Path(temp)
        executable = work / 'waveforms_driver'
        cmd = ['gfortran','-std=f2018','-O2','-fcheck=all',
               '-Werror=implicit-interface','-J',str(work),'-I',str(work),
               str(ROOT / 'src/kinds.f90'),
               str(ROOT / 'src/signal_waveforms.f90'),
               str(ROOT / 'tools/signal_reference/waveforms_driver.f90'),
               '-o',str(executable)]
        subprocess.run(cmd,check=True,cwd=work,capture_output=True,text=True)
        proc = subprocess.run([str(executable)],input='\n'.join(lines)+'\n',
                              cwd=work,check=True,capture_output=True,text=True)
    outputs = iter(proc.stdout.splitlines())
    largest = 0.0
    for case, (mode, expected) in enumerate(zip(mode_list,reference)):
        if mode == 5:
            state, value = next(outputs).split()
            if int(state) != 0:
                raise AssertionError(f'case {case}: cutoff status {state}')
            error = abs(float(value)-expected)
            largest = max(largest,error)
            if not np.isclose(float(value), expected, atol=4e-12,rtol=2e-12):
                raise AssertionError(f'case {case}: cutoff error {error:.3e}')
        else:
            state, ntext = map(int,next(outputs).split())
            if state != 0 or ntext != len(expected):
                raise AssertionError(f'case {case}: state={state}, n={ntext}, expected={len(expected)}')
            rows = [list(map(float,next(outputs).split())) for _ in range(ntext)]
            if mode == 2:
                result = np.array([complex(v[0],v[1]) for v in rows])
            elif mode == 4:
                result = np.array(rows).reshape((-1,3))
            else:
                result = np.array([v[0] for v in rows])
            error = float(np.max(np.abs(result-expected))) if result.size else 0.0
            largest = max(largest,error)
            if not np.allclose(result,expected,atol=6e-12,rtol=3e-12):
                raise AssertionError(f'case {case}: mode {mode} error {error:.3e};\n'
                                     f'got={result}\nref={expected}')
    try:
        extra = next(outputs)
        raise AssertionError(f'unexpected driver output: {extra}')
    except StopIteration:
        pass
    print(f'Passed {len(reference)} randomized SciPy waveform comparisons')
    print(f'Maximum absolute error: {largest:.3e}')


if __name__ == '__main__':
    main()
