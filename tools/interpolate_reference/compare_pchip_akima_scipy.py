#!/usr/bin/env python3
"""Independent numerical SciPy 1-D cubic interpolation regression.

Compiles native Fortran sources and tests PCHIP / Akima / Makima / cubic Hermite against
scipy.interpolate's PchipInterpolator and Akima1DInterpolator. No library
Fortran sources or compiled artifacts are generated inside the project.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import Akima1DInterpolator, CubicHermiteSpline, PchipInterpolator

ROOT = Path(__file__).resolve().parents[2]


def cases(seed: int):
    rng = np.random.default_rng(seed)
    for flavor in (1, 2, 3, 4):
        for n in (2, 3, 4, 5, 8, 13, 31, 65):
            for pattern in range(15):
                # Test strongly nonuniform spacing, flat intervals, sign changes,
                # and ties as well as smooth and random data.
                h = rng.uniform(0.1, 2.0, size=n - 1)
                if pattern % 4 == 0:
                    h = 10.0 ** rng.uniform(-1.0, 1.0, size=n - 1)
                x = np.r_[0.0, np.cumsum(h)]
                x -= rng.uniform(0.0, 5.0)
                if pattern == 0:
                    y = np.full(n, rng.uniform(-3.0, 3.0))
                elif pattern == 1:
                    y = rng.uniform(-2.0, 2.0) + rng.uniform(-3.0, 3.0)*x
                elif pattern == 2:
                    y = (x - np.mean(x))**2
                elif pattern == 3:
                    y = np.sin(x * 1.4)
                elif pattern == 4:
                    y = np.cumsum(rng.uniform(0.0, 1.0, n))
                elif pattern == 5:
                    y = np.cumsum(-rng.uniform(0.0, 1.0, n))
                elif pattern == 6:
                    y = rng.integers(-2, 3, n).astype(float)
                elif pattern == 7:
                    y = np.cumsum(rng.normal(0.0, 0.001, n))
                elif pattern == 8:
                    y = np.full(n, 0.0)
                    y[n // 2 :] = 2.0
                elif pattern == 9:
                    y = 0.1*x**3 - 0.2*x**2
                else:
                    y = rng.normal(size=n)
                q = np.r_[x, (x[:-1]+x[1:])/2.0,
                          x[0] - 0.25*h[0], x[-1] + 0.2*h[-1],
                          rng.uniform(x[0]-h[0], x[-1]+h[-1], 16)]
                ext = (-1, 0, 1)[pattern % 3]
                a = float(rng.uniform(x[0]-h[0], x[-1]+h[-1]))
                b = float(rng.uniform(x[0]-h[0], x[-1]+h[-1]))
                dydx = rng.normal(size=n)
                if pattern == 0:
                    dydx.fill(0.0)
                if pattern == 1:
                    dydx.fill((y[-1]-y[0])/(x[-1]-x[0]))
                yield (flavor, n, x, y, q, ext, a, b, dydx)


def run(build_dir: Path, flags: list[str], *, seed: int) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    compiler = os.environ.get('FC', 'gfortran')
    sources = (ROOT / 'src/kinds.f90', ROOT / 'src/interpolate_bspline.f90', ROOT / 'src/interpolate_regular_grid.f90', ROOT/'src/interpolate_piecewise.f90', ROOT / 'src/interpolate_barycentric.f90', ROOT/'src/interpolate_global.f90', ROOT / 'src/interpolate_kdtree.f90', ROOT / 'src/interpolate_rbf.f90', ROOT / 'src/interpolate_legacy.f90', ROOT / 'src/interpolate_scattered.f90', ROOT / 'src/interpolate_clough_tocher.f90', ROOT / 'src/interpolate_fitpack_bivariate.f90', ROOT / 'src/interpolate_fitpack_compat.f90', ROOT / 'src/interpolate.f90',
               ROOT / 'tools/interpolate_reference/pchip_akima_driver.f90')
    proc = subprocess.run([compiler, '-std=f2018', *flags,
                           '-J' + str(build_dir), '-I' + str(build_dir),
                           *(str(src) for src in sources), '-o', str(build_dir / 'reference')],
                          cwd=build_dir, capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(f'Fortran compile failure:\n{proc.stderr}')
    testcases = list(cases(seed))
    payload = [str(len(testcases))]
    for flavor, n, x, y, q, ext, a, b, dydx in testcases:
        payload.extend([f'{n} {len(q)} {flavor} {ext}',
                        ' '.join(format(v, '.17g') for v in x),
                        ' '.join(format(v, '.17g') for v in y)])
        if flavor == 4:
            payload.append(' '.join(format(v, '.17g') for v in dydx))
        payload.extend([' '.join(format(v, '.17g') for v in q),
                        f'{a:.17g} {b:.17g}'])
    proc = subprocess.run([str(build_dir / 'reference')],
                          input='\n'.join(payload)+'\n', capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError('Driver failed:\n' + proc.stderr)
    lines = iter(proc.stdout.splitlines())
    maximum_absolute = 0.0
    comparisons = 0
    failures = []
    for index, (flavor, n, x, y, q, ext, a, b, dydx) in enumerate(testcases, 1):
        header = next(lines).split()
        if header != ['CASE', str(index), '0']:
            raise RuntimeError(f'Unexpected case header {header!r} for {index}')
        got = np.array([[float(z) for z in next(lines).split()] for _ in q])
        integral = float(next(lines))
        coefficients = np.array([[float(z) for z in next(lines).split()]
                                 for _ in range(n - 1)]).T
        extr = None if ext < 0 else bool(ext)
        if flavor == 1:
            ref = PchipInterpolator(x, y, extrapolate=extr)
        elif flavor in (2, 3):
            ref = Akima1DInterpolator(x, y, method='akima' if flavor == 2 else 'makima',
                                      extrapolate=extr)
        else:
            ref = CubicHermiteSpline(x, y, dydx, extrapolate=extr)
        expected = np.array([ref(q, nu=nu) for nu in range(5)]).T
        values = [
            ('derivatives', got, expected, 5.0e-9, 1.5e-10),
            ('integral', np.array([integral]), np.array([ref.integrate(a, b)]),
             1.0e-8, 2.0e-10),
            ('coefficients', coefficients, ref.c[::-1], 5.0e-9, 1.5e-10),
        ]
        for label, lhs, rhs, atol, rtol in values:
            if lhs.shape != rhs.shape:
                raise AssertionError(f'{label} shape {lhs.shape} != {rhs.shape}')
            equal = np.isclose(lhs, rhs, atol=atol, rtol=rtol, equal_nan=True)
            if not np.all(equal):
                where = np.argwhere(~equal)[0]
                bad = tuple(where)
                failures.append((index, flavor, n, label, bad, lhs[bad], rhs[bad]))
            finite = np.isfinite(lhs) & np.isfinite(rhs)
            if np.any(finite):
                maximum_absolute = max(maximum_absolute,
                                       float(np.max(np.abs(lhs[finite]-rhs[finite]))))
            comparisons += lhs.size
    print(f'{len(testcases)} SciPy reference cases, {comparisons} compared values, '
          f'largest absolute difference {maximum_absolute:.3g}')
    if failures:
        print('FIRST FAILURES:', *failures[:12], sep='\n  ')
        raise AssertionError(f'{len(failures)} reference comparisons failed')
    print('PASS PCHIP/Akima/Makima/Hermite Fortran-SciPy regression')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=20261031)
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_interpolate_') as directory:
        root = Path(directory)
        if args.mode in ('checked', 'both'):
            print('CHECKED:', flush=True)
            run(root/'checked', ['-O0', '-g', '-Wall', '-Wextra', '-fcheck=all',
                                 '-ffpe-trap=invalid,zero,overflow'], seed=args.seed)
        if args.mode in ('optimized', 'both'):
            print('OPTIMIZED:', flush=True)
            run(root/'optimized', ['-O3', '-Wall', '-Wextra'], seed=args.seed)


if __name__ == '__main__':
    main()
