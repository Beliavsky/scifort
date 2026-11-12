#!/usr/bin/env python3
"""Compare SciFort vectorized neighborhood batches with scipy.ndimage."""
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
MODES = ('reflect', 'nearest', 'mirror', 'wrap', 'constant')


def main():
    i, j, k = np.indices((5, 4, 3))
    i, j, k = i + 1, j + 1, k + 1
    x = np.asfortranarray(np.sin(.31 * i + .23 * j) + .17 * k + .09 * i * j)
    z = x + 1j * (.37 * x + .21)
    sparse = np.array([[1, 0, 1], [0, 1, 1]], bool)
    cases = []
    for mode in MODES:
        cases.append((f'real {mode}', ndi.vectorized_filter(x, np.mean, size=(2, 3, 2),
                      mode=mode, **({'cval': -1.75} if mode == 'constant' else {}))))
        cases.append((f'complex {mode}', ndi.vectorized_filter(z, np.mean, size=(2, 3, 2),
                      mode=mode, **({'cval': -1.75 + .6j} if mode == 'constant' else {}))))

    def weighted(w, *, axis):
        return np.sum(w * np.arange(1, 5).reshape((1,) * (w.ndim - 1) + (4,)), axis=axis)

    cases.append(('sparse reordered footprint', ndi.vectorized_filter(
        x, weighted, footprint=sparse, axes=(2, 0), origin=(0, 1), mode='nearest')))
    cases.append(('valid real', ndi.vectorized_filter(x, np.mean, size=(2, 3, 2), mode='valid')))
    cases.append(('valid complex', ndi.vectorized_filter(z, np.mean, size=(2, 3, 2), mode='valid')))
    cases.append(('sparse selected valid', ndi.vectorized_filter(
        x, weighted, footprint=sparse, axes=(2, 0), mode='valid')))
    cases.append(('shifted origins', ndi.vectorized_filter(
        x, np.mean, size=(2, 3, 2), mode='mirror', origin=(-1, 1, 0))))
    nan_x = x.copy(order='F')
    nan_x[1, 2, 0] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        cases.append(('nanmean', ndi.vectorized_filter(
            nan_x, np.nanmean, size=(2, 2, 3), mode='constant', cval=np.nan)))
    cases.append(('no selected axes', x))

    with tempfile.TemporaryDirectory(prefix='scifort-vectorized-') as tmp_str:
        tmp = Path(tmp_str)
        sources = ['src/kinds.f90', 'src/signal_spline.f90',
                   'src/signal_spline_nd.f90', 'src/ndimage_morphology.f90',
                   'src/ndimage_distance.f90', 'src/ndimage_chamfer.f90',
                   'src/ndimage_measurements.f90', 'src/ndimage_watershed.f90',
                   'src/ndimage.f90',
                   'tools/ndimage_reference/vectorized_filter_driver.f90']
        cmd = ['gfortran', '-std=f2018', '-O0', '-fcheck=all',
               '-Werror=implicit-interface', '-J', str(tmp), '-I', str(tmp)]
        cmd += [str(ROOT / s) for s in sources] + ['-o', str(tmp / 'driver')]
        subprocess.run(cmd, check=True, cwd=tmp)
        subprocess.run([str(tmp / 'driver')], check=True, cwd=tmp,
                       stdout=subprocess.PIPE, text=True)
        with (tmp / 'vectorized_reference.bin').open('rb') as f:
            total = 0
            errors = {}
            for name, expected in cases:
                dtype = np.complex128 if np.iscomplexobj(expected) else np.float64
                data = f.read(expected.size * np.dtype(dtype).itemsize)
                if len(data) != expected.size * np.dtype(dtype).itemsize:
                    raise AssertionError(f'{name}: truncated binary output')
                observed = np.frombuffer(data, dtype).reshape(expected.shape, order='F')
                if not np.array_equal(np.isnan(observed), np.isnan(expected)):
                    raise AssertionError(f'{name}: NaN positions do not match')
                err = float(np.nanmax(np.abs(observed - expected)))
                if err > 1e-11:
                    index = np.unravel_index(np.nanargmax(np.abs(observed - expected)), expected.shape)
                    raise AssertionError(f'{name}: error={err}, index={index}, '
                                         f'got={observed[index]}, expected={expected[index]}')
                errors[name] = err
                total += expected.size
            if f.read(1):
                raise AssertionError('Trailing bytes in reference data')
    print(f'PASS: {total} SciPy vectorized_filter elements, {len(cases)} cases')
    print(f'Maximum absolute error: {max(errors.values()):.3g}')


if __name__ == '__main__':
    main()
