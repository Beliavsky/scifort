#!/usr/bin/env python3
"""Compare SciFort multidimensional recursive PCHIP with scipy.interpolate.

Includes ascending/descending axes, nonuniform spacing, nonlinear/flat data,
multichannel outputs, 1D--6D, extrapolation, fill, bounds checks, and NaNs.

Run: python tools/interpolate_reference/compare_tensor_pchip_scipy.py
     python tools/interpolate_reference/compare_tensor_pchip_scipy.py --seed 8675309
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[2]


def make_cases(seed):
    rng = np.random.default_rng(seed)
    for ndim, number in ((1, 80), (2, 70), (3, 44), (4, 24), (5, 10), (6, 4)):
        for case in range(number):
            axes = []
            for d in range(ndim):
                n = 4 + int(rng.integers(0, 4 if ndim <= 3 else 2))
                spacing = rng.uniform(0.12, 1.5, n - 1)
                if case % 7 == 0:
                    spacing[:] = 0.15 + 0.1 * d
                if case % 7 == 1:
                    spacing[0] *= 0.01
                axis = rng.uniform(-3.0, 3.0) + np.r_[0.0, np.cumsum(spacing)]
                if (case + d) % 3 == 0:
                    axis = axis[::-1].copy()
                axes.append(axis)
            nchan = 1 + case % 3
            shape = tuple(len(a) for a in axes) + (nchan,)
            vals = rng.normal(size=shape)
            if case % 5 == 0:
                vals[...] = 3.25  # constant plateaus
            elif case % 5 == 1:
                vals[...] = 0.
                for d,a in enumerate(axes):
                    vals += (a*a).reshape((1,)*d + (len(a),) + (1,)*(ndim-d)) / (d+1)
            elif case % 5 == 2:
                vals[::2, ...] = vals[1::2, ...][-1:, ...] if shape[0] % 2 else 0.0
            bounds = case % 7 == 0
            fill = case % 3
            q = np.empty((9, ndim))
            for d, a in enumerate(axes):
                lo, hi = min(a), max(a)
                width = hi - lo
                q[:, d] = rng.uniform(lo-0.3*width, hi+0.3*width, len(q))
                q[0,d], q[1,d] = a[0], a[-1]
                q[2,d], q[3,d] = a[1], a[-2]
                q[4,d] = (a[0]+a[1])*0.5
                q[5,d] = lo-0.1*width
                q[6,d] = hi+0.1*width
            if case % 11 == 2:
                q[-1,0] = np.nan
            yield axes, vals, q, bounds, fill


def serialize(samples):
    lines = [str(len(samples))]
    def fstring(arr):
        return ' '.join(f'{float(x):.17g}' for x in np.asarray(arr).ravel(order='F'))
    for axes, vals, q, bounds, fill in samples:
        nd, nch = len(axes), vals.shape[-1]
        # method code 5 is PCHIP in the shared driver.
        lines.append(f'{nd} {nch} {len(q)} 5 {int(bounds)} {fill} 0')
        for a in axes:
            lines.append(str(len(a)))
            lines.append(fstring(a))
        lines.append(' '.join(['-1']*nd))
        lines.append(fstring(vals))
        lines.append(fstring(np.zeros_like(vals)))
        for row in q:
            lines.append(fstring(row))
    return '\n'.join(lines) + '\n'


def check(flags, seed):
    samples = list(make_cases(seed))
    with tempfile.TemporaryDirectory(prefix='scifort_pchip_') as td:
        work = Path(td)
        sources = [ROOT / 'src/kinds.f90', ROOT / 'src/interpolate_bspline.f90',
                   ROOT / 'src/interpolate_regular_grid.f90', ROOT/'src/interpolate_piecewise.f90', ROOT / 'src/interpolate_barycentric.f90', ROOT/'src/interpolate_global.f90', ROOT / 'src/interpolate_kdtree.f90', ROOT / 'src/interpolate_rbf.f90', ROOT / 'src/interpolate_legacy.f90', ROOT / 'src/interpolate_scattered.f90', ROOT / 'src/interpolate_clough_tocher.f90', ROOT / 'src/interpolate_fitpack_bivariate.f90', ROOT / 'src/interpolate_curfit_core.f90', ROOT / 'src/interpolate_fitpack_compat.f90', ROOT / 'src/interpolate.f90',
                   ROOT / 'tools/interpolate_reference/tensor_grid_driver.f90']
        cmd = [os.environ.get('FC', 'gfortran'), '-std=f2018', *flags,
               '-J'+str(work),'-I'+str(work),*map(str,sources),'-o',str(work/'test')]
        cp = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
        if cp.returncode:
            raise RuntimeError(f'Build failed: {cp.stderr[-5000:]}')
        cp = subprocess.run([str(work / 'test')], input=serialize(samples), cwd=work,
                            capture_output=True, text=True)
        if cp.returncode:
            raise RuntimeError(f'Native driver failed: {cp.stderr[-5000:]}\n{cp.stdout[-1000:]}')
        lines = iter(cp.stdout.splitlines())
        count = 0
        max_abs = max_scaled = 0.0
        for cid, (axes, vals, q, bounds, fill) in enumerate(samples,1):
            header = next(lines).split()
            if header[:3] != ['CASE',str(cid),'0']:
                raise AssertionError(f'case {cid} setup: {header}')
            code = int(header[3])
            actual = np.empty((len(q), vals.shape[-1]), dtype=complex)
            for j in range(len(q)):
                numbers = np.fromstring(next(lines), sep=' ')
                if len(numbers) != 2*vals.shape[-1]:
                    raise AssertionError(f'case {cid} malformed output: {numbers}')
                actual[j] = numbers[::2] + 1j*numbers[1::2]
            fill_value = None if fill == 2 else (np.nan if fill == 0 else -17.375)
            with warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                ref = RegularGridInterpolator(tuple(axes), vals, method='pchip',
                                              bounds_error=bounds, fill_value=fill_value)
                try:
                    expected = ref(q)
                    expected_code = 0
                except ValueError as e:
                    if 'out of bounds' not in str(e):
                        raise
                    expected_code = 2
            if code != expected_code:
                raise AssertionError(f'case {cid} wrong status {code} vs {expected_code}')
            if code:
                continue
            diff = np.abs(actual-expected)
            tolerance = 5e-10+5e-10*np.abs(expected)
            valid = (diff <= tolerance) | (np.isnan(actual.real) & np.isnan(expected.real))
            if not np.all(valid):
                i,j = np.argwhere(~valid)[0]
                raise AssertionError(f'case {cid}: {len(axes)}D, channel {j}, point {i}, '
                                     f'got {actual[i,j]}, expected {expected[i,j]}, '
                                     f'axes={axes}, query={q[i]}')
            count += actual.size
            max_abs = max(max_abs, float(np.nanmax(diff, initial=0.)))
            max_scaled = max(max_scaled, float(np.nanmax(diff/(1+np.abs(expected)), initial=0.)))
        remaining = list(lines)
        if remaining:
            raise AssertionError(f'extra output: {remaining[:3]}')
    print(f'{len(samples)} scenarios, {count} scalar values, max abs {max_abs:.3g}, '
          f'max scaled {max_scaled:.3g}; seed={seed}; flags={" ".join(flags)}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=20261101)
    args = parser.parse_args()
    for flags in (['-O0','-g','-fcheck=all','-fbacktrace'],['-O3']):
        check(flags, args.seed)

if __name__ == '__main__':
    main()
