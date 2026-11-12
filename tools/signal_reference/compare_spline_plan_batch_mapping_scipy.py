#!/usr/bin/env python3
"""Compare batched Fortran spline-plan nonlinear warps with SciPy 1.17+."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
MODES = ['mirror', 'constant', 'grid-constant', 'nearest',
         'reflect', 'grid-mirror', 'wrap', 'grid-wrap']


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-spline-map-batch-') as directory:
        tmp = Path(directory)
        flags = ['gfortran', '-std=f2018', '-O0', '-fcheck=all',
                 '-ffpe-trap=invalid,zero,overflow', '-J', str(tmp), '-I', str(tmp)]
        objects = []
        for source in ('kinds.f90', 'signal_spline.f90', 'signal_spline_nd.f90'):
            obj = tmp / (Path(source).stem + '.o')
            subprocess.run([*flags, '-c', str(ROOT / 'src' / source), '-o', str(obj)],
                           check=True, cwd=tmp)
            objects.append(str(obj))
        executable = tmp / 'compare_batch'
        subprocess.run([*flags, *objects,
            str(ROOT / 'tools/signal_reference/spline_plan_batch_mapping_driver.f90'),
            '-o', str(executable)], check=True, cwd=tmp)
        lines = subprocess.check_output([str(executable)], text=True, cwd=tmp).splitlines()

    ii = np.arange(1, 20)[:, None]
    jj = np.arange(1, 18)[None, :]
    x = np.sin(0.17 * ii + 0.31 * jj) + 0.004 * ii * jj
    z = x + 1j * (-0.37 * x + 0.29)
    oi = np.arange(1, 13)[:, None]
    oj = np.arange(1, 12)[None, :]
    coords = np.empty((2, 12, 11))
    coords[0] = 0.52 + 0.86 * oi - 0.13 * oj + 0.11 * np.sin(0.13 * oj) - 1
    coords[1] = 0.72 + 0.22 * oi + 0.79 * oj - 0.17 * np.cos(0.21 * oi) - 1

    actual = {}
    for line in lines:
        fields = line.split()
        if len(fields) != 7:
            raise ValueError(f'Unexpected Fortran output: {line!r}')
        mi, order, kind, i, j = map(int, fields[:5])
        key = (mi, order, kind)
        if key not in actual:
            actual[key] = np.empty((12, 11), dtype=np.complex128)
        actual[key][i-1, j-1] = float(fields[5]) + 1j * float(fields[6])
    assert len(actual) == 8 * 6 * 2, len(actual)

    count, maximum, reflect = 0, 0.0, []
    for (mi, order, kind), got in actual.items():
        mode = MODES[mi-1]
        data = x if kind == 0 else z
        fill = 1.125 if kind == 0 else 1.125-0.35j
        ref = ndimage.map_coordinates(data, coords, mode=mode, order=order,
                                      cval=fill, prefilter=True)
        err = np.max(np.abs(got - ref))
        tolerance = 3e-8
        if mode in ('reflect', 'grid-mirror') and order > 1:
            tolerance = 1e-5  # documented SciPy recursion vs exact collocation
            reflect.append(err)
        if err > tolerance:
            raise AssertionError(f'{mode} degree={order} complex={kind}: {err:g} > {tolerance:g}')
        if tolerance == 3e-8:
            maximum = max(maximum, float(err))
        count += got.size * (2 if kind else 1)
    print(f'PASS: {count} batched nonlinear spline scalar SciPy comparisons')
    print(f'Max error outside short reflect special cases: {maximum:.3g}')
    if reflect:
        print(f'Max reflect/grid-mirror high-order discrepancy: {max(reflect):.3g}')


if __name__ == '__main__':
    main()
