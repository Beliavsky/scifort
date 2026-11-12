#!/usr/bin/env python3
"""Compare SciFort explicit FIR/SOS decimation against SciPy 1.17 dlti filters."""
from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy import signal


def compare(executable: str) -> None:
    rng = np.random.default_rng(20261028)
    cases = []
    firs = [
        np.array([0.15, 0.7, 0.15]),
        np.array([0.04, 0.12, 0.24, 0.28, 0.20, 0.12]),
        np.array([0.12, -0.05, 0.38, 0.31, 0.24]),
        np.array([1.0]),
    ]
    sections = [signal.butter(2, 0.35, output='sos'),
                signal.cheby1(4, 0.3, 0.36, output='sos'),
                signal.butter(6, 0.3, output='sos')]
    for kind in (0, 1):
        for phase in (0, 1):
            for q in (2, 3, 4):
                for n in (49, 64, 91, 137):
                    for ftype in range(1, 8):
                        x = rng.normal(size=n) + 0.2*np.cos(np.arange(n)*0.17)
                        if kind:
                            x = x + 1j*(rng.normal(size=n) + 0.3*np.sin(np.arange(n)*0.3))
                        taps = firs[ftype-1] if ftype <= len(firs) else np.empty(0)
                        sos = sections[ftype-len(firs)-1] if ftype > len(firs) else np.empty((0, 6))
                        cases.append((kind, ftype <= len(firs), phase, q, x, taps, sos))
    lines = [str(len(cases))]
    for kind, is_fir, phase, q, x, h, sos in cases:
        lines.append(f'{kind} {1 if is_fir else 2} {phase} {len(x)} {q} {len(h)} {len(sos)}')
        if len(h):
            lines.append(' '.join(f'{v:.17e}' for v in h))
        for row in sos:
            lines.append(' '.join(f'{v:.17e}' for v in row))
        lines.append(' '.join(f'{v:.17e}' for v in x.real))
        lines.append(' '.join(f'{v:.17e}' for v in x.imag))
    run = subprocess.run([executable], input='\n'.join(lines)+'\n', text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         timeout=90, check=True)
    results = iter(run.stdout.splitlines())
    compared = 0
    maximum = 0.0
    for i, (kind, is_fir, phase, q, x, h, sos) in enumerate(cases):
        if is_fir:
            filt = signal.dlti(h, [1.0])
        else:
            filt = signal.dlti(*signal.sos2tf(sos))
        reference = signal.decimate(x, q, ftype=filt, zero_phase=bool(phase))
        status, nout = map(int, next(results).split())
        if status != 0 or nout != len(reference):
            raise AssertionError(f'case {i}: wrong status/length {(status,nout)} != {len(reference)}')
        actual = np.array([complex(*map(float, next(results).split())) for _ in range(nout)])
        diff = abs(actual - reference)
        maximum = max(maximum, float(np.max(diff)))
        compared += nout
        if np.max(diff) > 1.0e-10:
            k = int(np.argmax(diff))
            raise AssertionError(f'case {i}: FIR={is_fir} phase={phase} q={q} n={len(x)} '
                                 f'element={k} actual={actual[k]} expected={reference[k]} '
                                 f'error={diff[k]:.4g}')
    print(f'PASS: {len(cases)} custom-filter SciPy decimation cases / '
          f'{compared} samples; maximum absolute error {maximum:.3e}', flush=True)


if __name__ == '__main__':
    compare(sys.argv[1])
