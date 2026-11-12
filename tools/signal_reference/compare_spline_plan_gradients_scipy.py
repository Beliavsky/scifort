#!/usr/bin/env python3
"""Compare analytic cached-spline gradients against SciPy central differences."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
MODES = ('mirror', 'constant', 'grid-constant', 'nearest',
         'reflect', 'grid-mirror', 'wrap', 'grid-wrap')
POINTS = np.array([
    [1.21, 1.37], [3.42, 2.71], [6.68, 4.21], [-0.31, 1.83],
    [8.34, 2.34], [2.14, -1.23], [3.29, 7.41], [-3.27, 8.31],
    [5.32, 4.89], [4.0, 3.0]], dtype=np.float64).T - 1.0


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-grad-ref-') as td:
        wd = Path(td)
        flags = ['gfortran', '-std=f2018', '-O0', '-fcheck=all',
                 '-ffpe-trap=invalid,zero,overflow', '-J', str(wd), '-I', str(wd)]
        objects = []
        for filename in ('kinds.f90', 'signal_spline.f90', 'signal_spline_nd.f90'):
            obj = wd / (Path(filename).stem + '.o')
            subprocess.run([*flags, '-c', str(ROOT / 'src' / filename), '-o', str(obj)],
                           check=True, cwd=wd)
            objects.append(str(obj))
        executable = wd / 'gradient_reference'
        subprocess.run([*flags, *objects,
                        str(ROOT / 'tools/signal_reference/spline_plan_gradient_driver.f90'),
                        '-o', str(executable)], check=True, cwd=wd)
        stdout = subprocess.check_output([str(executable)], text=True, cwd=wd)

    ii = np.arange(1, 8)[:, None]
    jj = np.arange(1, 6)[None, :]
    x = np.sin(0.21*ii + 0.37*jj) + 0.01*ii*jj
    data = [x, x+1j*(0.3*x+0.04), x[:1]]
    count = 0
    maximum = 0.0
    reflect_max = 0.0
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) != 11:
            raise ValueError(f'Unexpected driver result: {line!r}')
        mi, degree, kind, filtered, k = (int(v) for v in fields[:5])
        nums = [float(v) for v in fields[5:]]
        got = np.array([nums[0]+1j*nums[1], nums[2]+1j*nums[3],
                        nums[4]+1j*nums[5]])
        source = data[kind]
        mode = MODES[mi-1]
        cval = (0.25-0.13j) if kind==1 else 0.25
        coord = POINTS[:, k-1:k]
        def interpolated(at):
            return ndimage.map_coordinates(source, at, order=degree,
                mode=mode, cval=cval, prefilter=bool(filtered))[0]
        step = 1e-5
        expected = np.empty(3, dtype=complex)
        expected[0] = interpolated(coord)
        for axis in range(2):
            positive = coord.copy()
            negative = coord.copy()
            positive[axis] += step
            negative[axis] -= step
            expected[axis+1]=(interpolated(positive)-interpolated(negative))/(2*step)
        error = float(np.max(np.abs(got-expected)))
        tolerance = 4e-5
        if bool(filtered) and degree>1 and mode in ('reflect','grid-mirror'):
            tolerance = 5e-3  # documented short-line recursion discrepancy
        if error > tolerance:
            raise AssertionError((mi, degree, kind, filtered, k, error,
                                  got, expected, tolerance))
        if tolerance > 1e-3:
            reflect_max = max(reflect_max, error)
        else:
            maximum = max(maximum, error)
        count += 3 * (2 if kind==1 else 1)
    assert count == 8 * 6 * 3 * 2 * 10 * 3 + 8 * 6 * 2 * 10 * 3, count
    print(f'PASS: {count} real/imag gradient and value SciPy comparisons')
    print(f'Maximum regular-case error {maximum:.3g}; short reflect {reflect_max:.3g}')


if __name__ == '__main__':
    main()
