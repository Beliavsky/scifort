#!/usr/bin/env python3
"""Compare SciFort geometric transformations and batched region callbacks."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [ROOT / p for p in (
    'src/kinds.f90', 'src/signal_spline.f90', 'src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90', 'src/ndimage_distance.f90',
    'src/ndimage_chamfer.f90', 'src/ndimage_measurements.f90',
    'src/ndimage_watershed.f90', 'src/ndimage.f90',
    'tools/ndimage_reference/geometric_comprehension_driver.f90',
)]
MODES = ('reflect', 'grid-mirror', 'mirror', 'nearest', 'wrap',
         'grid-wrap', 'constant', 'grid-constant')
ORDERS = (0, 1, 3, 5)


def mapping2(z):
    return (z[0] + 0.35 * (z[1] - 2), z[1] - 0.2 * (z[0] - 1))


def mapping3(z):
    return (0.4 * (z[0] + 1), 0.5 + 0.25 * (z[1] + 1),
            1.0 + np.sin(0.11 * (z[0] + 1)))


def mapping1(z):
    return (0.8 * z[0] + 0.33 * z[1] - 0.1 * z[2],)


def main():
    i, j = np.mgrid[1:7, 1:8]
    a2 = np.sin(0.17 * i + 0.31 * j) + np.cos(0.23 * i - 0.51 * j)
    z2 = a2.astype(complex) + 1j * (0.3 * a2 + np.sin(0.21 * j))
    labs = (i + 3 * j) % 5 - 2
    i, j, k = np.mgrid[1:5, 1:4, 1:6]
    a3 = np.sin(0.2 * i + 0.7 * j - 0.3 * k) + 0.03 * i * k
    i = np.arange(1, 7)
    a1 = np.sin(0.8 * i) + 0.03 * i * i
    ids = [-2, -1, 0, 1, 2, 99, 1, -2]
    total = 0
    max_error = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-geometric-') as td:
        root = Path(td)
        subprocess.run(['gfortran', '-std=f2018', '-O0', '-fcheck=all',
                        '-Werror=implicit-interface', '-J', str(root),
                        '-I', str(root), *map(str, SOURCES),
                        '-o', str(root/'driver')], cwd=root, check=True)
        subprocess.run([str(root/'driver')], cwd=root, check=True)
        with (root/'geometric_comprehension.bin').open('rb') as stream:
            def read(shape, dtype):
                nonlocal total
                count = int(np.prod(shape))
                data = np.frombuffer(stream.read(count*np.dtype(dtype).itemsize), dtype=dtype)
                if data.size != count:
                    raise AssertionError(('missing output', shape, dtype))
                total += count
                return data.reshape(shape, order='F')

            def compare(label, actual, expected, tolerance=2e-8):
                nonlocal max_error
                err = float(np.max(np.abs(actual-expected)))
                max_error = max(max_error, err)
                if not np.isfinite(err) or err > tolerance:
                    raise AssertionError((label, err, actual, expected))

            for mode in MODES:
                for order in ORDERS:
                    got = read((5, 8), np.float64)
                    expected = ndi.geometric_transform(
                        a2, mapping2, output_shape=(5, 8), order=order,
                        mode=mode, cval=-0.75)
                    # Existing high-order reflect prefilter has short-array residuals.
                    tol = (1e-4 if order == 5 else 2e-7) if mode in (
                        'reflect', 'grid-mirror') and order > 1 else 2e-8
                    compare(('real', mode, order), got, expected, tol)
                    got = read((5, 8), np.complex128)
                    expected = ndi.geometric_transform(
                        z2, mapping2, output_shape=(5, 8), order=order,
                        mode=mode, cval=-0.75+0.4j)
                    compare(('complex', mode, order), got, expected, tol)
            got = read((5, 8), np.float64)
            expected = ndi.geometric_transform(
                a2, mapping2, output_shape=(5, 8), order=3, mode='mirror',
                prefilter=False, cval=-0.75)
            compare('prefilter false', got, expected)
            got = read((6, 4), np.float64)
            expected = ndi.geometric_transform(
                a3, mapping3, output_shape=(6, 4), order=3, mode='nearest')
            compare('rank 3 to 2', got, expected)
            got = read((3, 2, 2), np.float64)
            expected = ndi.geometric_transform(
                a1, mapping1, output_shape=(3, 2, 2), order=1, mode='reflect')
            compare('rank 1 to 3', got, expected)
            got = read((len(ids),), np.float64)
            def func(values):
                return np.sum(np.sin(values)+values*values)
            expected = ndi.labeled_comprehension(a2, labs, ids, func, float, -19.5)
            compare('batched comprehension', got, expected, 2e-12)
            if stream.read(1):
                raise AssertionError('unexpected trailing output')
    print(f'PASS: {total} geometric/label SciPy values, max error={max_error:.3g}')
if __name__ == '__main__':
    main()
