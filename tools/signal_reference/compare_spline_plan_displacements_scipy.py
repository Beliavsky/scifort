#!/usr/bin/env python3
"""Verify cached real/complex displacement-field warps against SciPy."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
MODES = ('mirror', 'constant', 'grid-constant', 'nearest',
         'reflect', 'grid-mirror', 'wrap', 'grid-wrap')


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-warp-ref-') as folder:
        tmp = Path(folder)
        flags = ['gfortran', '-std=f2018', '-O0', '-fcheck=all',
                 '-ffpe-trap=invalid,zero,overflow', '-J', str(tmp), '-I', str(tmp)]
        objects = []
        for filename in ('kinds.f90', 'signal_spline.f90', 'signal_spline_nd.f90'):
            obj = tmp / (Path(filename).stem + '.o')
            subprocess.run([*flags, '-c', str(ROOT / 'src' / filename), '-o', str(obj)],
                           check=True, cwd=tmp)
            objects.append(str(obj))
        exe = tmp / 'displacement_reference'
        subprocess.run([*flags, *objects,
                        str(ROOT / 'tools/signal_reference/spline_plan_displacements_driver.f90'),
                        '-o', str(exe)], check=True, cwd=tmp)
        stdout = subprocess.check_output([str(exe)], text=True, cwd=tmp)

    actual = {}
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) != 7:
            raise ValueError(f'Unexpected Fortran output: {line!r}')
        mode, degree, kind, i, j = map(int, fields[:5])
        key = mode, degree, kind
        if key not in actual:
            actual[key] = np.empty((27,) if mode == 0 else (14, 12), dtype=complex)
        if mode == 0:
            actual[key][i-1] = complex(float(fields[5]), float(fields[6]))
        else:
            actual[key][i-1, j-1] = complex(float(fields[5]), float(fields[6]))
    assert len(actual) == 8 * 6 * 2 + 1, len(actual)

    ii = np.arange(1, 17)[:, None]
    jj = np.arange(1, 14)[None, :]
    x = np.sin(0.17*ii + 0.29*jj) + 0.002*ii*jj
    z = x + 1j*(0.32*x+0.1)
    oi = np.arange(1, 15)[:, None]
    oj = np.arange(1, 13)[None, :]
    coordinates = np.empty((2, 14, 12))
    coordinates[0] = oi + 0.21*np.sin(0.4*(oi+oj)) - 1
    coordinates[1] = oj - 0.31*np.cos(0.27*(2*oi-oj)) - 1

    count, maximum, reflect_max = 0, 0.0, 0.0
    for (mi, order, kind), got in actual.items():
        if mi == 0:
            affine_coords = np.vstack((
                1+0.1+0.45*np.arange(27)+0.1*np.sin(0.31*np.arange(1,28)),
                1-0.3+0.27*np.arange(27)+0.15*np.cos(0.23*np.arange(1,28)),
            )) - 1
            ref = ndimage.map_coordinates(x, affine_coords, order=3,
                                          mode='grid-wrap', prefilter=True)
            tol = 3e-8
        else:
            mode = MODES[mi-1]
            data = x if kind == 0 else z
            fill = -0.4 if kind == 0 else (-0.4+0.25j)
            ref = ndimage.map_coordinates(data, coordinates, order=order, mode=mode,
                                          cval=fill, prefilter=True)
            tol = 3e-8
            if mode in ('reflect', 'grid-mirror') and order > 1:
                tol = 1e-5  # known SciPy recursion vs exact finite-line solve
        err = float(np.max(np.abs(got-ref)))
        if err > tol:
            raise AssertionError(f'{mi=} {order=} {kind=} {err=} {tol=}')
        if tol == 1e-5:
            reflect_max = max(reflect_max, err)
        else:
            maximum = max(maximum, err)
        count += got.size * (2 if kind else 1)
    print(f'PASS: {count} real/complex displacement-field SciPy scalar comparisons')
    print(f'Max normal-case error: {maximum:.3g}; max short reflect discrepancy: {reflect_max:.3g}')


if __name__ == '__main__':
    main()
