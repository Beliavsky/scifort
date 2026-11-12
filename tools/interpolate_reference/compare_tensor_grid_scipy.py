#!/usr/bin/env python3
"""Randomized comparisons of SciFort rectilinear tensor-product splines with SciPy.

Cubic and quintic interpolation uses independent tensor-product LU solves in
SciFort. SciPy's sparse iterative reference is configured with tighter-than-
default tolerances for fair comparisons with a direct solve.

Run: python tools/interpolate_reference/compare_tensor_grid_scipy.py
     python tools/interpolate_reference/compare_tensor_grid_scipy.py --seed 8675309
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
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[2]
METHODS = ('slinear', 'cubic', 'quintic')


def make_cases(seed):
    rng = np.random.default_rng(seed)
    for ndim in range(1, 5):
        for method in METHODS:
            degree = {'slinear': 1, 'cubic': 3, 'quintic': 5}[method]
            number = 22 if ndim <= 2 else (12 if ndim == 3 else 7)
            for case in range(number):
                axes = []
                for d in range(ndim):
                    n = degree + 1 + int(rng.integers(0, 3 if ndim <= 2 else 2))
                    # Regular/strongly irregular spacing; some negative axes.
                    spacing = 10.0 ** rng.uniform(-0.5, 0.5, n-1)
                    if case % 5 == 0:
                        spacing[:] = 0.5 + 0.3 * d
                    x = rng.uniform(-3, 1) + np.r_[0., np.cumsum(spacing)]
                    if (d+case) % 3 == 0:
                        x = x[::-1].copy()
                    axes.append(x)
                nchan = 1 + case % 3
                complex_data = case % 2 == 1
                shape = tuple(len(a) for a in axes) + (nchan,)
                vals = rng.standard_normal(size=shape)
                if complex_data:
                    vals = vals + 1j*rng.standard_normal(size=shape)
                if case % 4 == 0:
                    # Low-degree polynomials should interpolate exactly,
                    # even with descending and irregular grid axes.
                    vals[...] = 1.0
                    for d,a in enumerate(axes):
                        shp = (1,)*d + (len(a),) + (1,)*(ndim-d)  # channels dimension
                        vals += 0.2*a.reshape(shp)
                    if complex_data:
                        vals = vals.astype(complex) + 0.3j*vals
                fill = case % 3  # 0: nan, 1: constant, 2: extrapolate
                bounds = case % 6 == 0
                nq = 9
                q = np.empty((nq, ndim))
                for d,a in enumerate(axes):
                    lo,hi = min(a),max(a)
                    w = hi-lo
                    q[:,d] = rng.uniform(lo-0.12*w,hi+0.12*w,nq)
                    q[0,d] = a[0]
                    q[1,d] = a[-1]
                    q[2,d] = a[len(a)//2]
                    q[3,d] = (a[0]+a[1])/2
                    q[4,d] = lo-(0.2 if ndim <= 2 else 0.05)*w
                    q[5,d] = hi+(0.3 if ndim <= 2 else 0.07)*w
                if case % 7 == 1:
                    q[-1,0] = np.nan
                nu = None
                if case % 3 != 0:
                    nu = np.zeros(ndim, dtype=int)
                    if case % 5 == 2:
                        nu[case % ndim] = degree+1  # higher-than-polynomial derivative is zero
                    elif case % 5 == 3 and ndim > 1:
                        nu[0] = 1
                        nu[-1] = 2 if degree > 1 else 1
                    else:
                        nu[case % ndim] = min(degree,1+(case//3) % 3)
                yield axes, vals, q, method, bounds, fill, complex_data, nu


def serialize(samples):
    lines = [str(len(samples))]
    def fstring(a):
        return ' '.join(f'{float(v):.17g}' for v in np.asarray(a).ravel(order='F'))
    for axes,vals,q,method,bounds,fill,is_complex,nu in samples:
        nd,nchan = len(axes),vals.shape[-1]
        lines.append(f'{nd} {nchan} {len(q)} {2+METHODS.index(method)} {int(bounds)} {fill} {int(is_complex)}')
        for a in axes:
            lines.append(str(len(a)))
            lines.append(fstring(a))
        lines.append(' '.join(map(str, nu)) if nu is not None else ' '.join(['-1']*nd))
        lines.append(fstring(vals.real))
        lines.append(fstring(vals.imag if is_complex else np.zeros_like(vals)))
        for row in q:
            lines.append(fstring(row))
    return '\n'.join(lines)+'\n'


def check(flags, seed):
    samples = list(make_cases(seed))
    with tempfile.TemporaryDirectory(prefix='scifort_tensor_grid_') as td:
        work = Path(td)
        sources = [ROOT/'src/kinds.f90', ROOT/'src/interpolate_bspline.f90',
                   ROOT/'src/interpolate_regular_grid.f90', ROOT/'src/interpolate_piecewise.f90', ROOT/'src/interpolate_barycentric.f90', ROOT/'src/interpolate_global.f90', ROOT/'src/interpolate_kdtree.f90', ROOT/'src/interpolate_rbf.f90', ROOT/'src/interpolate_legacy.f90', ROOT/'src/interpolate_scattered.f90', ROOT/'src/interpolate_clough_tocher.f90', ROOT/'src/interpolate_fitpack_bivariate.f90', ROOT/'src/interpolate_fitpack_compat.f90', ROOT/'src/interpolate.f90',
                   ROOT/'tools/interpolate_reference/tensor_grid_driver.f90']
        cmd = [os.environ.get('FC','gfortran'), '-std=f2018', *flags,
               '-J'+str(work),'-I'+str(work),*map(str,sources),'-o',str(work/'test')]
        cp = subprocess.run(cmd,cwd=work,text=True,capture_output=True)
        if cp.returncode:
            raise RuntimeError('Failed to compile: '+cp.stderr[-5000:])
        cp = subprocess.run([str(work/'test')], input=serialize(samples),
                            text=True,capture_output=True,cwd=work)
        if cp.returncode:
            raise RuntimeError('Native test aborted: '+cp.stderr[-4500:]+'\n'+cp.stdout[-500:])
        lines = iter(cp.stdout.splitlines())
        count = 0
        maxdiff = 0.
        maxscaled = 0.
        for cid,(axes, vals, q, method, bounds, fill, is_complex, nu) in enumerate(samples,1):
            fields = next(lines).split()
            if fields[:2] != ['CASE',str(cid)]:
                raise AssertionError(f'malformed case header: {fields}')
            if fields[2] != '0':
                raise AssertionError(f'case {cid} setup failed: {fields}; method={method} axes={axes}')
            gotstatus = int(fields[3])
            got = np.empty((len(q),vals.shape[-1]),dtype=complex)
            for j in range(len(q)):
                numbers = np.fromstring(next(lines),sep=' ')
                if numbers.size != 2*vals.shape[-1]:
                    raise AssertionError(f'case {cid} row {j}: invalid {numbers}')
                got[j,:] = numbers[::2]+1j*numbers[1::2]
            fill_value = None if fill == 2 else (
                np.nan if fill == 0 else (-17.375+3.25j if is_complex else -17.375))
            with warnings.catch_warnings():
                warnings.simplefilter('error',RuntimeWarning)
                ref = RegularGridInterpolator(tuple(axes),vals,method=method,
                                              bounds_error=bounds,fill_value=fill_value,
                                              solver=spsolve)
                try:
                    expected = ref(q, nu=nu).astype(complex)
                    expstatus = 0
                except ValueError as err:
                    if 'out of bounds' not in str(err):
                        raise
                    expstatus = 2
            if gotstatus != expstatus:
                raise AssertionError(f'case {cid} status {gotstatus} != {expstatus}')
            if expstatus != 0:
                continue
            # With NaNs SciPy gives NaN+0j even for complex input.
            difference = np.abs(got-expected)
            tolerance = 2e-9+2e-9*np.abs(expected)
            if not np.all((difference <= tolerance)|(
                    np.isnan(got.real)&np.isnan(expected.real))):
                w = np.nanargmax(difference/tolerance)
                i,j = np.unravel_index(w,got.shape)
                raise AssertionError(f'case {cid} method={method} ndim={len(axes)} index={(i,j)}: '
                                     f'got={got[i,j]} expected={expected[i,j]} diff={difference[i,j]} '
                                     f'grid_shape={vals.shape}, axes={axes}, q={q[i,:]}, nu={nu}, fill={fill}')
            if difference.size:
                maxdiff = max(maxdiff,float(np.nanmax(difference,initial=0.)))
                maxscaled = max(maxscaled,float(np.nanmax(difference/(1+np.abs(expected)),initial=0.)))
            count += got.size
        remaining=list(lines)
        if remaining:
            raise AssertionError(f'unexpected additional output {remaining[:2]}')
    print(f'{len(samples)} cases, {count} comparisons, max abs error {maxdiff:.3g}, '
          f'max scaled {maxscaled:.3g} [seed={seed}, flags={" ".join(flags)}]')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=20261031)
    args=parser.parse_args()
    for flags in (['-O0','-g','-fcheck=all','-fbacktrace'],['-O3']):
        check(flags,args.seed)

if __name__=='__main__':
    main()
