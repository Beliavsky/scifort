#!/usr/bin/env python3
"""Check general_cosine and get_window DPSS against SciPy 1.17."""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.signal import get_window, windows

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261021)
CASES = []

for n in (0, 1, 2, 3, 4, 5, 8, 9, 16, 17, 32, 33, 64, 65, 128):
    for symmetric in (0, 1):
        for coeff in ([], [0], [1], [0.5, 0.5], [0.42, 0.5, 0.08],
                      [0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368],
                      [1.1, -0.4, 0.2, -0.1]):
            for family in (1, 2):
                CASES.append((family, n, symmetric, list(coeff), 0.0))
        for nw in (0.3, 0.5, 1.0, 1.5, 2.5, 3.5, 4.0, 6.0):
            if n > 1 and nw >= n/2:
                continue
            for family in (3, 4):
                CASES.append((family, n, symmetric, [], nw))

for _ in range(280):
    n = int(RNG.integers(1, 115))
    symmetric = int(RNG.integers(2))
    if RNG.random() < 0.62:
        family = int(RNG.choice([1, 2]))
        k = int(RNG.integers(0, 11))
        coeff = RNG.normal(0, 0.65, k).tolist()
        nw = 0.0
    else:
        family = int(RNG.choice([3, 4]))
        coeff = []
        nw = float(RNG.uniform(0.15, max(0.16, n/2 - 0.01))) if n > 1 else 2.0
        if n > 1 and nw >= n/2:
            nw = n/4
    CASES.append((family, n, symmetric, coeff, nw))


def scipy_reference(case):
    family, n, symmetric, coeff, nw = case
    if family == 1:
        return windows.general_cosine(n, coeff, sym=bool(symmetric))
    if family == 2:
        # SciFort's established get_window convention accepts an empty output.
        if n == 0:
            return windows.general_cosine(n, coeff, sym=bool(symmetric))
        return get_window(('general_cosine', coeff), n, fftbins=not symmetric)
    if n < 1:
        return None
    return get_window(('dpss', nw), n, fftbins=not symmetric)


def run():
    with tempfile.TemporaryDirectory(prefix='scifort-cosine-dpss-') as td:
        directory = Path(td)
        sources = ['kinds.f90', 'signal_dpss.f90',
                   'signal_windows_extra.f90', 'signal_spectral.f90']
        cmd = ['gfortran', '-std=f2018', '-O2', '-fcheck=bounds',
               '-J', str(directory), '-I', str(directory),
               *[str(ROOT/'src'/x) for x in sources],
               str(ROOT/'reference/signal_cosine_dpss_reference_driver.f90'),
               '-o', str(directory/'test')]
        subprocess.run(cmd, check=True, cwd=directory)
        inp = []
        for kind, n, symmetric, coeff, nw in CASES:
            inp.append(f'{kind} {n} {symmetric} {len(coeff)} {nw:.17g}')
            if coeff:
                inp.append(' '.join(f'{a:.17g}' for a in coeff))
        output = subprocess.run([str(directory/'test')], input='\n'.join(inp)+'\n',
            text=True, capture_output=True, check=True, cwd=directory).stdout.splitlines()

    assert len(output) == len(CASES), (len(output), len(CASES))
    count, max_error = 0, 0.0
    for case, line in zip(CASES, output):
        toks = line.split()
        status, n = int(toks[0]), int(toks[1])
        result = np.array([float(x.replace('D', 'E')) for x in toks[2:]])
        expected = scipy_reference(case)
        if expected is None:
            assert status == 1 and n == 0, (case, status, n)
        else:
            assert status == 0 and n == len(result) == len(expected), (case, status, n)
            err = np.max(abs(result-expected), initial=0)
            max_error = max(max_error, float(err))
            if not np.allclose(result, expected, atol=8e-10, rtol=8e-10):
                j = int(np.argmax(abs(result-expected)))
                raise AssertionError((case, err, j, result[j], expected[j]))
        count += 1
    print(f'PASS: {count}/{len(CASES)} SciPy 1.17 cosine/DPSS window cases; max absolute error {max_error:.3e}')


if __name__ == '__main__':
    run()
