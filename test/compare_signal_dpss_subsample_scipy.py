#!/usr/bin/env python3
"""Independent SciPy 1.17 checks for DPSS subsample normalization.

Usage: python test/compare_signal_dpss_subsample_scipy.py [checked|optimized]
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal.windows import dpss

ROOT = Path(__file__).resolve().parents[1]
MODE = sys.argv[1] if len(sys.argv) > 1 else 'checked'
if MODE not in ('checked', 'optimized'):
    raise SystemExit('expected checked or optimized')
FLAGS = ['-std=f2018', '-Werror=implicit-interface']
FLAGS += ['-O0', '-g', '-fcheck=all'] if MODE == 'checked' else ['-O3']


def cases():
    rng = np.random.default_rng(20261027)
    out = []
    for m in (2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 24, 25, 31, 32,
              33, 48, 63, 64, 65, 97, 128, 129):
        for nw in (0.5, 1.1, 2.5, 3.5):
            if nw >= m / 2:
                continue
            for sym in (True, False):
                out.append((m, nw, min(m, 4, max(1, int(2*nw))), sym))
    for _ in range(72):
        m = int(rng.integers(3, 115))
        nw = float(rng.uniform(0.15, 0.49*m))
        k = int(rng.integers(1, min(5, m, max(1, int(2*nw)))+1))
        out.append((m, nw, k, bool(rng.integers(2))))
    return out


def main():
    selected = cases()
    with tempfile.TemporaryDirectory(prefix='scifort-dpss-subsample-') as d:
        directory = Path(d)
        driver = ROOT / 'reference' / 'signal_dpss_reference_driver.f90'
        exe = directory / 'run'
        subprocess.run(['gfortran', *FLAGS, '-J', str(directory), '-I', str(directory),
                        str(ROOT/'src/kinds.f90'), str(ROOT/'src/signal_dpss.f90'),
                        str(driver), '-o', str(exe)], check=True, cwd=directory)
        request = ''.join(f'{m} {nw:.17g} {k} {int(sym)} 2\n'
                          for m, nw, k, sym in selected)
        output = subprocess.run([str(exe)], input=request, text=True,
                                capture_output=True, check=True, cwd=directory).stdout
    it = iter(output.splitlines())
    max_taper = 0.0
    max_ratio = 0.0
    num_values = 0
    for idx, (m, nw, k, sym) in enumerate(selected):
        head = next(it).split()
        assert head == ['CASE', str(m), str(k), '0'], (idx, head)
        ratios = np.fromstring(next(it), sep=' ')
        got = np.array([np.fromstring(next(it), sep=' ') for _ in range(k)])
        expect, eratio = dpss(m, nw, Kmax=k, sym=sym, norm='subsample',
                              return_ratios=True)
        taper_error = np.max(np.abs(got - expect))
        ratio_error = np.max(np.abs(ratios - eratio))
        num_values += got.size
        max_taper = max(max_taper, float(taper_error))
        max_ratio = max(max_ratio, float(ratio_error))
        if taper_error > 3e-8 or ratio_error > 1e-8:
            raise AssertionError(f'case {idx}, {m,nw,k,sym}: '
                                 f'taper {taper_error:.4g}, ratio {ratio_error:.4g}')
    try:
        print('extra output:',next(it))
        raise AssertionError('Unexpected extra Fortran output')
    except StopIteration:
        pass
    print(f'DPSS subsample {MODE}: {len(selected)} SciPy cases, '
          f'{num_values} taper values, max absolute error {max_taper:.3g}, '
          f'max concentration-ratio error {max_ratio:.3g}')


if __name__ == '__main__':
    main()
