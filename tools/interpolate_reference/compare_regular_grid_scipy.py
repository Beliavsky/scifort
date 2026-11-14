#!/usr/bin/env python3
"""Check SciFort rectilinear-grid interpolation against SciPy 1.17.

Usage:
    python tools/interpolate_reference/compare_regular_grid_scipy.py
    python tools/interpolate_reference/compare_regular_grid_scipy.py --seed 123

Only NumPy and SciPy are used in this test harness, not in the Fortran library.
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[2]


def make_cases(seed: int):
    rng = np.random.default_rng(seed)
    for ndim in range(1, 7):
        for case in range(90):
            axes = []
            for d in range(ndim):
                n = int(rng.integers(1, 6 if ndim <= 3 else 4))
                # One or more singleton axes, varied spacing and orientation.
                if case % 17 == 0 and d == 0:
                    n = 1
                points = rng.uniform(-2, 2) + np.r_[0., np.cumsum(10.0 ** rng.uniform(-1.5, 1.1, n-1))]
                if (d + case) % 3 == 0:
                    points = points[::-1].copy()
                axes.append(points)
            nchan = 1 + case % 3
            iscomplex = (case // 2) % 2 == 1
            method = ('linear', 'nearest')[(case // 3) % 2]
            fill = (case // 5) % 3  # NaN, constant, or extrapolate
            bounds = (case % 9 == 0)
            shape = tuple(len(v) for v in axes) + (nchan,)
            vals = rng.normal(size=shape)
            if iscomplex:
                vals = vals + 1j * rng.normal(size=shape)
            # Include constant/affine functions as exactness checks.
            if case % 11 == 0:
                vals = np.zeros(shape, dtype=complex if iscomplex else float)
                for d, coords in enumerate(axes):
                    reshape = (1,) * d + (len(coords),) + (1,) * (ndim-d-1) + (1,)
                    vals += rng.uniform(-2,2) * coords.reshape(reshape)
            nq = 17
            x = np.empty((nq, ndim))
            for d, coords in enumerate(axes):
                a, b = sorted((coords[0],coords[-1]))
                span = max(0.1,b-a)
                x[:,d] = rng.uniform(a-0.8*span,b+0.8*span,size=nq)
                x[0,d] = coords[0]
                x[1,d] = coords[-1]
                x[2,d] = coords[len(coords)//2]
                x[3,d] = 0.5 * (coords[0]+coords[-1])
                x[4,d] = a-0.8*span
                x[5,d] = b+0.8*span
                if len(coords) > 1:
                    x[6,d] = 0.5 * (coords[0]+coords[1])
            if case % 7 == 0:
                x[-1,0] = np.nan
            yield axes, vals, x, method, bounds, fill, iscomplex


def format_values(values):
    return ' '.join(f'{float(v):.17g}' for v in np.asarray(values).ravel(order='F'))


def make_input(cases):
    lines = [str(len(cases))]
    for axes, values, x, method, bounds, fill, iscomplex in cases:
        ndim = len(axes)
        nchan = values.shape[-1]
        lines.append(f'{ndim} {nchan} {len(x)} {int(method == "nearest")} {int(bounds)} {fill} {int(iscomplex)}')
        for coords in axes:
            lines.append(str(len(coords)))
            lines.append(format_values(coords))
        lines.append(format_values(values.real))
        lines.append(format_values(values.imag if iscomplex else np.zeros_like(values)))
        for row in x:
            lines.append(format_values(row))
    return '\n'.join(lines) + '\n'


def run_test(flags, seed):
    cases = list(make_cases(seed))
    with tempfile.TemporaryDirectory(prefix='scifort_grid_') as td:
        tmp = Path(td)
        sources = [ROOT/'src/kinds.f90', ROOT/'src/interpolate_bspline.f90',ROOT/'src/interpolate_regular_grid.f90',
                   ROOT/'src/interpolate_piecewise.f90', ROOT/'src/interpolate_barycentric.f90', ROOT/'src/interpolate_global.f90', ROOT/'src/interpolate_kdtree.f90', ROOT/'src/interpolate_rbf.f90', ROOT/'src/interpolate_legacy.f90', ROOT/'src/interpolate_scattered.f90', ROOT/'src/interpolate_clough_tocher.f90', ROOT/'src/interpolate_fitpack_bivariate.f90', ROOT/'src/interpolate_curfit_core.f90', ROOT/'src/interpolate_fitpack_compat.f90', ROOT/'src/interpolate.f90',
                   ROOT/'tools/interpolate_reference/regular_grid_driver.f90']
        compiler = os.environ.get('FC', 'gfortran')
        cmd = [compiler, '-std=f2018', *flags, '-J'+str(tmp), '-I'+str(tmp),
               *map(str,sources), '-o', str(tmp/'driver')]
        result = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(f'Compilation failed:\n{result.stderr}')
        result = subprocess.run([str(tmp/'driver')], input=make_input(cases),
                                capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(f'Fortran failed:\n{result.stderr[-4500:]}\n{result.stdout[-1000:]}')
        lines = iter(result.stdout.splitlines())
        nvalues = 0
        maxerr = 0.0
        for case_id,(axes, values, x, method, bounds, fill, iscomplex) in enumerate(cases,1):
            fields = next(lines).split()
            if fields[:2] != ['CASE',str(case_id)] or fields[2] != '0':
                raise AssertionError(f'case {case_id} setup failed {fields}')
            gotstatus = int(fields[3])
            expected = np.empty((len(x),values.shape[-1]),dtype=complex)
            for j in range(len(x)):
                raw = np.fromstring(next(lines),sep=' ')
                if raw.size != 2*values.shape[-1]:
                    raise AssertionError((case_id,j,raw.size))
                expected[j,:] = raw[::2] + 1j*raw[1::2]
            fill_value = (None if fill == 2 else (np.nan if fill == 0 else
                          (-17.375 + 3.25j if iscomplex else -17.375)))
            ref = RegularGridInterpolator(tuple(axes),values,method=method,
                                          bounds_error=bounds,fill_value=fill_value)
            try:
                target = ref(x)
                target_status = 0
            except ValueError as err:
                if 'out of bounds' not in str(err):
                    raise
                target_status = 2
            if gotstatus != target_status:
                raise AssertionError(f'case {case_id} status {gotstatus} vs SciPy {target_status}')
            if target_status != 0:
                continue
            if not iscomplex:
                target = target.astype(complex)
            # Preserve nan+0j, which NumPy gives for nan query coordinates.
            difference = np.abs(expected-target)
            if not np.allclose(expected,target,rtol=5e-12,atol=5e-12,equal_nan=True):
                where = np.unravel_index(np.nanargmax(difference),difference.shape)
                raise AssertionError(f'case {case_id}: index {where}: got {expected[where]}, expected {target[where]}, '
                                     f'axes {axes}, method={method}, bounds={bounds}, fill={fill}')
            maxerr = max(maxerr,float(np.nanmax(difference,initial=0)))
            nvalues += target.size
        rest=list(lines)
        if rest: raise AssertionError(f'unexpected trailing lines: {rest[:3]}')
        print(f'{len(cases)} scenarios; {nvalues} result values; maximum abs error {maxerr:.3g}; '
              f'flags={" ".join(flags)}; seed={seed}')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=20261031)
    args=parser.parse_args()
    for flags in (['-O0','-g','-fcheck=all','-fbacktrace'], ['-O3']):
        run_test(flags,args.seed)


if __name__ == '__main__':
    main()
