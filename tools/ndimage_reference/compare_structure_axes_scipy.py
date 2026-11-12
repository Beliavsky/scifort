#!/usr/bin/env python3
"""Compare binary stencils, numeric component labels and derivative axes with SciPy."""
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
    'tools/ndimage_reference/structure_axes_driver.f90')]


def read_values(stream, shape, dtype):
    n = int(np.prod(shape))
    arr = np.frombuffer(stream.read(n * np.dtype(dtype).itemsize), dtype=dtype)
    if arr.size != n:
        raise AssertionError('truncated binary reference output')
    return arr.reshape(shape, order='F')


def first_derivative(x, axis, output, mode, cval):
    return ndi.correlate1d(x, [-.5, 0, .5], axis=axis,
                           output=output, mode=mode, cval=cval)


def second_derivative(x, axis, output, mode, cval):
    return ndi.correlate1d(x, [1., -2., 1.], axis=axis,
                           output=output, mode=mode, cval=cval)


def check_labels(f, x, connectivity):
    n = int(read_values(f, (1,), np.int32)[0])
    labels = read_values(f, x.shape, np.int32)
    structure = ndi.generate_binary_structure(x.ndim, connectivity)
    expected, n_expected = ndi.label(x, structure=structure)
    if n != n_expected:
        raise AssertionError(('number of labels', n, n_expected, x.shape))
    # SciFort enumerates components in column-major order; SciPy numbers
    # components in row-major order. Match partitions, not their numeric IDs.
    got_flat = labels.ravel(order='F')
    exp_flat = expected.ravel(order='F')
    if not np.array_equal(got_flat == 0, exp_flat == 0):
        raise AssertionError('foreground mask differs')
    if not np.array_equal(got_flat[:, None] == got_flat[None, :],
                          exp_flat[:, None] == exp_flat[None, :]):
        raise AssertionError('connected-component partitions differ')
    return labels.size


def check_real(f, x, derivative, axes, mode, cval=0.0):
    got = read_values(f, x.shape, np.float64)
    fun = ndi.generic_gradient_magnitude if derivative == 'gradient' else ndi.generic_laplace
    callback = first_derivative if derivative == 'gradient' else second_derivative
    expected = fun(x, callback, axes=axes, mode=mode, cval=cval)
    err = np.max(np.abs(got - expected))
    if err > 2e-12:
        raise AssertionError((derivative, x.shape, axes, mode, err))
    return got.size, err


def main():
    count = 0
    largest_error = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-structure-axes-') as td:
        td = Path(td)
        cmd = ['gfortran', '-std=f2018', '-O2', '-fcheck=all',
               '-ffpe-trap=invalid,zero,overflow', '-Werror=implicit-interface',
               '-J', str(td), '-I', str(td), *(str(src) for src in SOURCES),
               '-o', str(td / 'reference')]
        subprocess.run(cmd, cwd=td, check=True)
        subprocess.run([str(td / 'reference')], cwd=td, check=True)
        with (td / 'structure_axes_results.bin').open('rb') as f:
            for nd, conn in ((1, 1), (2, 1), (2, 2), (3, 2), (4, 3)):
                got = read_values(f, (3,) * nd, np.int32).astype(bool)
                expected = ndi.generate_binary_structure(nd, conn)
                if not np.array_equal(got, expected):
                    raise AssertionError(('generated structure', nd, conn))
                count += got.size
            cross = ndi.generate_binary_structure(2, 1)
            span = np.array([True, False, True])
            even = np.array([[True, True], [False, True]])
            asym = np.fromfunction(lambda i, j, k:
                                   ((2*(i+1)+3*(j+1)+5*(k+1)+(i+1)*(j+1)) % 7) > 3,
                                   (2,3,3), dtype=int)
            for original, niter in ((cross, 2), (cross, 3),
                                    (span, 2), (span, 3), (even, 3), (asym, 2)):
                expected = ndi.iterate_structure(original, niter)
                got = read_values(f, expected.shape, np.int32).astype(bool)
                if not np.array_equal(got, expected):
                    raise AssertionError(('iterated structure', original.shape, niter,
                                          got.astype(int), expected.astype(int)))
                count += got.size
            ii = np.fromfunction(lambda i, j:
                                 ((i+1)**2+3*(j+1)+2*(i+1)*(j+1))%7-4,
                                 (6,5), dtype=int)
            count += check_labels(f, ii, 1)
            count += check_labels(f, ii, 2)
            rr = np.fromfunction(lambda i,j,k:
                                 ((2*(i+1)+3*(j+1)+5*(k+1)+(i+1)*(j+1))%5-2)/3.0,
                                 (4,3,2), dtype=int)
            count += check_labels(f, rr, 1)
            count += check_labels(f, rr, 3)
            a2 = np.fromfunction(lambda i,j:
                                 np.sin(.32*(i+1)+.41*(j+1))+(i+1)*(j+1)/13.0,
                                 (5,6))
            a3 = np.fromfunction(lambda i,j,k:
                                 np.cos(.23*(i+1)+.14*(j+1)+.77*(k+1))+
                                 (i+1)*(j+1)*(k+1)/19.0, (4,3,2))
            cases = [
                (a2, 'gradient', (0,), 'nearest', 0.0),
                (a2, 'laplace', (1,), 'reflect', 0.0),
                (a2, 'gradient', (1, 0), 'wrap', 0.0),
                (a3, 'gradient', (0, 2), ('nearest','constant'), -0.7),
                (a3, 'laplace', (2,), ('grid-wrap',), 0.25),
                (a3, 'laplace', (1, 0), ('constant','wrap'), 1.3),
                (a3, 'gradient', (1, 2, 0), ('grid-wrap','nearest','reflect'), 0.5),
            ]
            for x, derivative, axes, mode, cval in cases:
                n, err = check_real(f, x, derivative, axes, mode, cval)
                count += n
                largest_error = max(largest_error, err)
            if f.read(1):
                raise AssertionError('extra bytes in reference output')
    print(f'PASS: {count} structure/label/derivative comparisons, '
          f'max derivative error {largest_error:.3g}')


if __name__ == '__main__':
    main()
