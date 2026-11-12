#!/usr/bin/env python3
"""Randomized SciPy-reference tests for weighted Lomb-Scargle and vector strength.

Requires SciPy, NumPy, and gfortran. No runtime dependency is added to SciFort.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.signal import lombscargle, vectorstrength

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261016
COUNT = 900


def main() -> None:
    rng = np.random.default_rng(SEED)
    chunks = []
    refs = []
    for run in range(COUNT):
        n = int(rng.integers(12, 80))
        m = int(rng.integers(1, 25))
        x = np.sort(rng.uniform(-10.0, 15.0, n))
        f = rng.uniform(0.15, 18.0, m)
        # Mix multiple harmonics, noise, offsets, zero weights, and unequal gaps.
        y = (rng.uniform(-3, 3) + rng.uniform(0.25, 3) *
             np.sin(f[0] * x + rng.uniform(-2, 2)) + rng.normal(0, 0.4, n))
        use_weights = run % 4 != 0
        precenter = run % 3 == 0
        floating = run % 2 == 0
        weights = rng.uniform(0.1, 3.0, n)
        if run % 7 == 0:
            weights[:n//5] = 0.0
        if run % 13 == 0:
            weights *= 1000.0
        kw = dict(precenter=precenter, floating_mean=floating)
        if use_weights:
            kw['weights'] = weights
        p = lombscargle(x, y, f, normalize='power', **kw)
        pn = lombscargle(x, y, f, normalize='normalize', **kw)
        a = lombscargle(x, y, f, normalize='amplitude', **kw)
        strength, phase = vectorstrength(x, 2*np.pi/f)
        refs.append(np.column_stack((p, pn, a.real, a.imag, strength, phase)))
        chunks.append(f'{n} {m} {int(precenter)} {int(floating)} {int(use_weights)}\n')
        for values in (x, y, f, weights):
            chunks.append(' '.join(f'{t:.17e}' for t in values) + '\n')

    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp) / 'irregular_test'
        cmd = ['gfortran', '-std=f2018', '-O2', '-fcheck=all',
               '-Wimplicit-interface', '-Werror=implicit-interface',
               '-J', tmp, '-I', tmp, '-o', str(exe),
               str(ROOT / 'src/kinds.f90'),
               str(ROOT / 'src/signal_irregular.f90'),
               str(ROOT / 'tools/signal_reference/irregular_driver.f90')]
        subprocess.run(cmd, check=True, cwd=tmp)
        run = subprocess.run([str(exe)], input=''.join(chunks), capture_output=True,
                             check=True, text=True, cwd=tmp)

    lines = iter(run.stdout.splitlines())
    max_abs = 0.0
    for idx, expected in enumerate(refs):
        count = int(next(lines))
        assert count == len(expected), (idx, count, len(expected))
        actual = np.array([[float(s) for s in next(lines).split()]
                           for _ in range(count)])
        err = np.abs(actual - expected)
        max_abs = max(max_abs, float(np.max(err)))
        # Phase is in radians and sensitive only when vector strength vanishes;
        # no test inputs here have truly zero resultant strength.
        np.testing.assert_allclose(actual, expected, rtol=2e-9, atol=2e-10,
                                   err_msg=f'case {idx}')
    print(f'PASS: {COUNT} randomized cases, {sum(len(v) for v in refs)} frequency bins, '
          f'max absolute error {max_abs:.3e}')


if __name__ == '__main__':
    main()
