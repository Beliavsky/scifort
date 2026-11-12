#!/usr/bin/env python3
"""Numerically compare SciFort b/a decimation against SciPy custom dlti filters."""
from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy import signal


def main(executable: str) -> None:
    rng = np.random.default_rng(28672026)
    designs = [
        # (numerator, denominator), deliberately including nonunit a0.
        ([0.3, 0.15, 0.08], [1.0, -0.54, 0.16]),
        ([0.25, 0.15, 0.25], [1.7, -0.78, 0.25]),
        ([0.12, 0.08, 0.05, 0.03], [1.0, -0.4]),
        ([0.08, -0.12, 0.25], [1.0, -0.3, 0.1]),
        ([0.06, 0.19, 0.50, 0.19, 0.06], [1.5]),
        ([0.06, 0.22, 0.47, 0.22, 0.06], [1.0]),
        ([0.1, -0.03, 0.36, 0.22], [2.0]),
        ([0.75], [1.2]),
    ]
    cases = []
    for mode in (0, 1, 2):
        for phase in (0, 1):
            for q in (2, 3, 5):
                for n in (41, 64, 91, 127):
                    for j, (b0, a0) in enumerate(designs):
                        b = np.array(b0, dtype=np.complex128)
                        a = np.array(a0, dtype=np.complex128)
                        if mode != 0:
                            # Some systems are complex FIR, others complex IIR.
                            b = b + 1j*np.array([0.045*np.cos(k+j) for k in range(len(b))])
                            a = a + 1j*np.array([0.05*np.sin(k+j) if k else 0.03
                                                  for k in range(len(a))])
                        x = rng.normal(size=n) + 0.25*np.cos(np.arange(n)*0.15)
                        if mode == 1:
                            x = x + 1j*(rng.normal(size=n) + 0.1*np.sin(np.arange(n)*0.27))
                        cases.append((mode, phase, q, x, b, a))
    lines = [str(len(cases))]
    for mode, phase, q, x, b, a in cases:
        lines.append(f'{mode} {phase} {len(x)} {q} {len(b)} {len(a)}')
        for vals in (b.real,b.imag,a.real,a.imag,x.real,np.imag(x)):
            lines.append(' '.join(f'{v:.17e}' for v in vals))
    proc = subprocess.run([executable], input='\n'.join(lines)+'\n', text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=180, check=True)
    it = iter(proc.stdout.splitlines())
    max_error = 0.0
    compared = 0
    for idx, (mode, phase, q, x, b, a) in enumerate(cases):
        tf = signal.dlti(b, a)
        ref = signal.decimate(x, q, ftype=tf, zero_phase=bool(phase))
        status, count = map(int, next(it).split())
        if status != 0 or count != len(ref):
            raise AssertionError(f'case {idx}: status/length={status,count}, expected {len(ref)}')
        out = np.array([complex(*map(float, next(it).split())) for _ in range(count)])
        err = abs(out-ref)
        delta = float(np.max(err))
        compared += count
        max_error = max(max_error, delta)
        if delta > 1.0e-10:
            k = int(np.argmax(err))
            raise AssertionError(f'case {idx}: mode={mode} phase={phase} q={q} '
                                 f'np={len(b)}, na={len(a)}, n={len(x)}, '
                                 f'index={k}, got={out[k]}, reference={ref[k]}, err={delta:.4g}')
    print(f'PASS: {len(cases)} SciPy transfer-function decimation cases / '
          f'{compared} output samples; maximum absolute error {max_error:.3e}',flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
