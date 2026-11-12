#!/usr/bin/env python3
"""Check complex SOS decimation against SciPy's sosfilt/sosfiltfilt."""
from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy import signal


def main(executable: str) -> None:
    rng = np.random.default_rng(20261028)
    cases = []
    for mode in (0, 1, 2):  # real x/complex SOS, complex x/complex SOS, complex x/real SOS
        for phase in (0, 1):
            for q in (2, 3, 4, 6):
                for n in (41, 65, 92, 128):
                    for sections in (1, 2, 3):
                        sos = np.zeros((sections, 6), dtype=np.complex128)
                        for j in range(sections):
                            sos[j] = [
                                0.19 + 0.03*j + 0.025*1j,
                                0.12 - 0.02*j - 0.018*1j,
                                0.03 + 0.006*j + 0.015*1j,
                                1.0,
                                -0.39 + 0.06*j + 0.042*1j,
                                0.09 - 0.016*j - 0.028*1j,
                            ]
                        # Exercise SOS padding correction for first-order sections.
                        if sections == 3:
                            sos[0, 2] = sos[0, 5] = 0.0
                        if mode == 2:
                            sos = sos.real.astype(np.complex128)
                        x = (rng.standard_normal(n) + 0.22*np.cos(np.arange(n)*0.17))
                        if mode != 0:
                            x = x + 1j*(rng.standard_normal(n) + 0.13*np.sin(np.arange(n)*0.14))
                        cases.append((mode, phase, q, x, sos))
    lines = [str(len(cases))]
    for mode, phase, q, x, sos in cases:
        lines.append(f'{mode} {phase} {len(x)} {q} {len(sos)}')
        for section in sos:
            for values in (section.real, section.imag):
                lines.append(' '.join(f'{v:.17e}' for v in values))
        for values in (np.real(x),np.imag(x)):
            lines.append(' '.join(f'{v:.17e}' for v in values))
    proc = subprocess.run([executable], input='\n'.join(lines)+'\n', text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          timeout=180, check=True)
    it = iter(proc.stdout.splitlines())
    compared, max_error = 0, 0.0
    for index, (mode, phase, q, x, sos) in enumerate(cases):
        ref = (signal.sosfiltfilt(sos, x) if phase else signal.sosfilt(sos, x))[::q]
        status, count = map(int, next(it).split())
        if status != 0 or count != len(ref):
            raise AssertionError(f'case {index}: status/count {status}/{count}; wanted {len(ref)}')
        observed = np.array([complex(*map(float, next(it).split())) for _ in range(count)])
        err = np.abs(observed-ref)
        delta = float(np.max(err))
        compared += count
        max_error = max(max_error, delta)
        if delta > 1.0e-10:
            k = int(np.argmax(err))
            raise AssertionError(f'case {index}: mode={mode} phase={phase} q={q} '
                                 f'sections={len(sos)} n={len(x)} index={k} '
                                 f'actual={observed[k]} expected={ref[k]} error={delta:.3g}')
    print(f'PASS: {len(cases)} SciPy complex SOS decimation cases / '
          f'{compared} samples, maximum absolute error {max_error:.3e}', flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
