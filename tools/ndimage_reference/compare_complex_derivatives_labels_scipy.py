#!/usr/bin/env python3
"""Compare complex generic derivatives and indexed region statistics to SciPy."""
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [ROOT / p for p in (
    'src/kinds.f90', 'src/signal_spline.f90', 'src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90', 'src/ndimage_distance.f90',
    'src/ndimage_chamfer.f90', 'src/ndimage_measurements.f90',
    'src/ndimage_watershed.f90', 'src/ndimage.f90',
    'tools/ndimage_reference/complex_derivatives_labels_driver.f90',
)]
MODES = ('reflect', 'grid-mirror', 'mirror', 'nearest',
         'wrap', 'grid-wrap', 'constant', 'grid-constant')
FILL = .37 - .6j


def first_derivative(x, axis, output, mode, cval):
    return ndi.correlate1d(x, [-.5, 0, .5], axis=axis, output=output,
                           mode=mode, cval=cval)


def second_derivative(x, axis, output, mode, cval):
    return ndi.correlate1d(x, [1, -2, 1], axis=axis, output=output,
                           mode=mode, cval=cval)


def assert_close(got, expected, context, atol=2e-12):
    np.testing.assert_allclose(got, expected, atol=atol, rtol=2e-12,
                               equal_nan=True, err_msg=str(context))
    return np.asarray(got).size


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-complex-deriv-label-') as td:
        tmp = Path(td)
        args = ['gfortran', '-std=f2018', '-O0', '-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow', '-Werror=implicit-interface',
                '-J', str(tmp), '-I', str(tmp)]
        objs = []
        for source in SOURCES[:-1]:
            obj = tmp / (source.stem + '.o')
            subprocess.run([*args, '-c', str(source), '-o', str(obj)],
                           cwd=tmp, check=True, capture_output=True, text=True)
            objs.append(str(obj))
        subprocess.run([*args, *objs, str(SOURCES[-1]), '-o', str(tmp/'driver')],
                       cwd=tmp, check=True, capture_output=True, text=True)
        subprocess.run([str(tmp/'driver')], cwd=tmp, check=True,
                       capture_output=True, text=True)
        payload = (tmp/'complex_derivatives_labels.bin').read_bytes()
    offset = 0
    total = 0
    def take(shape, dtype):
        nonlocal offset
        shape = tuple(shape)
        n = int(np.prod(shape))
        values = np.frombuffer(payload, dtype=dtype, count=n, offset=offset)
        assert values.size == n, ('truncated stream', shape)
        offset += n * np.dtype(dtype).itemsize
        return values.reshape(shape, order='F')

    i, j = np.mgrid[1:7, 1:6]
    a2 = np.sin(.25*i+.37*j) + 1j*np.cos(.43*i-.16*j)
    i, j, k = np.mgrid[1:5, 1:4, 1:3]
    a3 = np.sin(.29*i+.41*j+.14*k) + 1j*np.cos(.11*i-.17*j+.31*k)
    for mode in MODES:
        total += assert_close(take(a2.shape, np.complex128),
                              ndi.generic_gradient_magnitude(a2, first_derivative,
                                                             mode=mode, cval=FILL),
                              ('complex gradient', mode))
        total += assert_close(take(a2.shape, np.complex128),
                              ndi.generic_laplace(a2, second_derivative,
                                                  mode=mode, cval=FILL),
                              ('complex Laplace', mode))
    for p, mode in enumerate(MODES):
        selected_modes = (mode, MODES[(p+3) % 8])
        kwargs = dict(axes=(2, 0), mode=selected_modes, cval=FILL)
        total += assert_close(take(a3.shape, np.complex128),
                              ndi.generic_gradient_magnitude(a3, first_derivative, **kwargs),
                              ('selected complex gradient', selected_modes))
        total += assert_close(take(a3.shape, np.complex128),
                              ndi.generic_laplace(a3, second_derivative, **kwargs),
                              ('selected complex Laplace', selected_modes))
    i, j = np.mgrid[1:6, 1:5]
    x = (3*i-j)/5.0 + np.cos(i*j)
    labels = np.full(x.shape, 3, dtype=np.int32)
    labels[(2*i+3*j)%6==0] = -7
    labels[(2*i+3*j)%6==1] = 0
    labels[(2*i+3*j)%6==2] = 500000001
    index = np.array([3, -7, 0, 3, 500000001, -999, -7, 2**31-1])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for name, method in (
            ('sum', ndi.sum_labels), ('mean', ndi.mean),
            ('variance', ndi.variance), ('standard deviation', ndi.standard_deviation),
        ):
            ref = method(x, labels=labels, index=index)
            # Label absent -> NaN for mean/variance/std. The sums are zero.
            total += assert_close(take((8,), np.float64), ref, name)
        expected_centers = np.array(ndi.center_of_mass(x, labels, index)) + 1
        total += assert_close(take((2,8), np.float64).T, expected_centers, 'center of mass')
    ref_counts = [np.count_nonzero(labels==key) for key in index]
    np.testing.assert_array_equal(take((8,), np.int32), ref_counts)
    total += len(index)
    hist = take((7,8), np.int32)
    for n, label in enumerate(index):
        expected = ndi.histogram(x, -5., 5., 7, labels=labels, index=int(label))
        if expected is None:
            expected = np.zeros(7, dtype=int)
        np.testing.assert_array_equal(hist[:, n], expected)
        total += 7
    mins = take((8,), np.float64)
    maxs = take((8,), np.float64)
    minpos = take((2,8), np.int32)
    maxpos = take((2,8), np.int32)
    low = take((2,8), np.int32)
    high = take((2,8), np.int32)
    # The SciFort positions and bounding boxes use one-based coordinates.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        minimums, maximums, pmins, pmaxs = ndi.extrema(x, labels, index)
    total += assert_close(mins, minimums, 'region minima')
    total += assert_close(maxs, maximums, 'region maxima')
    for n, label in enumerate(index):
        total += assert_close(minpos[:,n], np.array(pmins[n])+1,
                              ('minimum position', label)) if np.any(labels==label) else 0
        total += assert_close(maxpos[:,n], np.array(pmaxs[n])+1,
                              ('maximum position', label)) if np.any(labels==label) else 0
        if np.any(labels==label):
            points = np.array(np.nonzero(labels==label))
            np.testing.assert_array_equal(low[:,n], points.min(axis=1)+1)
            np.testing.assert_array_equal(high[:,n], points.max(axis=1)+1)
        else:
            np.testing.assert_array_equal(low[:,n], 0)
            np.testing.assert_array_equal(high[:,n], 0)
        total += 4
    assert offset == len(payload), (offset, len(payload))
    print(f'PASS: {total} SciPy complex derivative and labeled-region values')


if __name__ == '__main__':
    main()
