#!/usr/bin/env python3
"""Numerical comparison for SciFort tensor-product B-spline evaluation.

Run with: python tools/signal_reference/compare_signal_spline2d_eval_scipy.py
Requires GNU Fortran, NumPy, and SciPy. Coordinates are shifted by 1 to
compare Fortran one-based indices with SciPy zero-based map_coordinates.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261021)


def main():
    cases = []
    for degree in (2, 3):
        for complex_input in (False, True):
            for m, n in ((1, 1), (1, 13), (14, 1), (2, 3), (3, 2),
                         (7, 9), (23, 27), (35, 16)):
                for trial in range(8):
                    coeff = RNG.normal(size=(m, n))
                    if complex_input:
                        coeff = coeff + 1j * RNG.normal(size=(m, n))
                    # Independent scalars, in-bounds positions, mirrored negatives,
                    # and coordinates beyond multiple reflection periods.
                    q = 70
                    coords = np.vstack((RNG.uniform(-3*m, 4*m, size=q),
                                        RNG.uniform(-3*n, 4*n, size=q)))
                    coords[:, 0:4] = [[0, m-1, -1, 2*m], [0, n-1, -2, 2*n]]
                    expected = map_coordinates(coeff, coords, order=degree,
                                               mode='mirror', prefilter=False)
                    cases.append((degree, coeff, coords, expected))
    sources = [ROOT / 'src' / (n + '.f90') for n in
               ('kinds', 'signal_spline', 'signal_spline2d_eval')]
    driver = ROOT / 'tools/signal_reference/spline2d_eval_reference_driver.f90'
    with tempfile.TemporaryDirectory(prefix='scifort-spline2d-eval-') as temp_dir:
        temp = Path(temp_dir)
        objs = []
        for src in sources:
            obj = temp / (src.stem + '.o')
            subprocess.run(['gfortran', '-std=f2018', '-O2', '-fcheck=bounds',
                            '-I', str(temp), '-J', str(temp), '-c', str(src),
                            '-o', str(obj)], check=True, cwd=temp)
            objs.append(str(obj))
        exe = temp / 'spline2d_eval_driver'
        subprocess.run(['gfortran', '-std=f2018', '-O2', '-I', str(temp), '-J', str(temp),
                        str(driver), *objs, '-o', str(exe)], check=True, cwd=temp)
        lines = [str(len(cases))]
        for degree, coeff, coords, _ in cases:
            m, n = coeff.shape
            lines.append(f'{m} {n} {degree} {int(np.iscomplexobj(coeff))} {coords.shape[1]}')
            lines.extend(f'{v.real:.17e} {v.imag:.17e}' for v in coeff.T.flat)
            lines.extend(f'{y+1:.17e} {x+1:.17e}' for y, x in coords.T)
        proc = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                              text=True, check=True, capture_output=True, cwd=temp)
        out = iter(proc.stdout.splitlines())
        worst = 0.0
        for degree, coeff, coords, expected in cases:
            status, size = map(int, next(out).split())
            if status != 0 or size != expected.size:
                raise AssertionError(f'bad status/size: {status=}, {size=}')
            got = np.array([complex(*map(float, next(out).split()))
                            for _ in range(size)])
            err = float(np.max(np.abs(got - expected)))
            if err > 2e-12:
                loc = np.argmax(np.abs(got-expected))
                raise AssertionError(f'{degree=}, {coeff.shape=}, '
                                     f'complex={np.iscomplexobj(coeff)}, '
                                     f'error={err} at {coords[:,loc]}: '
                                     f'{got[loc]} vs {expected[loc]}')
            worst = max(worst, err)
        print(f'PASS: {len(cases)}/{len(cases)} cases, '
              f'{sum(len(c[3]) for c in cases)} point evaluations; '
              f'max absolute error {worst:.3g}')


if __name__ == '__main__':
    main()
