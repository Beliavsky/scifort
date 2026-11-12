"""Compare reusable Lomb-Scargle plans against independent SciPy results.

Run: python test/compare_signal_irregular_plan_scipy.py [checked|optimized]
Builds the driver and its dependencies in an isolated temporary directory.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import lombscargle

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "reference" / "signal_irregular_plan_reference_driver.f90"


def cases():
    return [(n, nf, nc, variant) for n in (9, 24, 49, 103)
            for nf in (1, 7, 21) for nc in (1, 4) for variant in range(8)]


def make_data(n, nf, nc, variant):
    index = np.arange(1, n + 1)
    x = 0.17*(index-1)+0.009*index*index+0.07*np.sin(0.21*index)
    w = 0.17+0.14*((index*5)%11)
    fidx = np.arange(1, nf + 1)
    freq = 0.15+0.19*fidx+0.003*fidx*fidx
    y = np.empty((2, n, nc), dtype=float)
    for channel in range(1, nc + 1):
        for k in (1, 2):
            y[k-1, :, channel-1] = (0.27*np.cos(0.7*x+0.14*k*channel)
                + np.sin(1.8*x-0.3*channel)+0.05*index/n+0.27*k)
    if variant == 6:
        w[index%4 == 0] = 0.
    return x, freq, (None if variant in (3, 7) else w), y


def reference(n, nf, nc, variant):
    x, freq, weights, y = make_data(n, nf, nc, variant)
    mode = ('amplitude' if variant in (2, 5, 7) else
            'normalize' if variant in (1, 6) else 'power')
    output = np.empty((2, nf, nc), dtype=complex if mode == 'amplitude' else float)
    for channel in range(nc):
        for k in range(2):
            output[k, :, channel] = lombscargle(x, y[k, :, channel], freq,
                weights=weights, precenter=variant in (4, 5, 7),
                floating_mean=variant in (1, 2, 6), normalize=mode)
    return output


def main():
    build = sys.argv[1] if len(sys.argv) > 1 else 'checked'
    if build not in ('checked', 'optimized'):
        raise SystemExit('expected checked or optimized')
    flags = ['-std=f2018', '-Werror=implicit-interface']
    flags += (['-O0', '-g', '-fcheck=all', '-fbacktrace'] if build == 'checked'
              else ['-O3'])
    work = cases()
    with tempfile.TemporaryDirectory(prefix='scifort-lomb-plan-') as temp:
        tmp = Path(temp)
        sources = [ROOT/'src'/name for name in
                   ('kinds.f90', 'signal_irregular.f90', 'signal_irregular_nd.f90')]
        command = [os.getenv('FC', 'gfortran'), *flags, *map(str, sources),
                   str(DRIVER), '-o', str(tmp/'reference')]
        subprocess.run(command, cwd=tmp, capture_output=True, text=True, check=True)
        proc = subprocess.run([str(tmp/'reference')], cwd=tmp,
            input=''.join(' '.join(map(str, c))+'\n' for c in work),
            capture_output=True, text=True, check=True)
    lines = iter(proc.stdout.splitlines())
    max_abs, max_scaled, total = 0., 0., 0
    for n, nf, nc, variant in work:
        header = next(lines).split()
        assert header == ['CASE', str(n), str(nf), str(nc), str(variant), '0'], header
        is_complex = variant in (2, 5, 7)
        data = []
        for _ in range(2*nf*nc):
            row = next(lines)
            data.append(complex(*map(float, row.split())) if is_complex else float(row))
        actual = np.asarray(data).reshape((2, nf, nc), order='F')
        expected = reference(n, nf, nc, variant)
        error = np.max(np.abs(actual-expected))
        scaled = error/max(1., np.max(np.abs(expected)))
        if not np.isfinite(scaled) or scaled > 1e-10:
            raise AssertionError((n, nf, nc, variant, error, scaled))
        max_abs = max(max_abs, error)
        max_scaled = max(max_scaled, scaled)
        total += actual.size
    assert next(lines, None) is None
    print(f'PASS {build}: {len(work)} independent SciPy cases, '
          f'{total} values, max abs {max_abs:.3g}, max scaled {max_scaled:.3g}')


if __name__ == '__main__':
    main()
