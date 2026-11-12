#!/usr/bin/env python3
"""Compare SciFort multichannel vectorized filters with SciPy 1.17+."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]


def features(w, *, axis):
    return np.stack((np.mean(w, axis=axis), np.mean(w * w, axis=axis)), axis=-1)


def main():
    i, j, k = np.indices((5, 4, 3))
    i, j, k = i + 1, j + 1, k + 1
    x = np.asfortranarray(np.sin(.31 * i + .23 * j) + .17 * k + .09 * i * j)
    z = x + 1j * (.37 * x + .21)
    sparse = np.array([[1, 0, 1], [0, 1, 1]], bool)
    cases = []
    for mode in ('reflect', 'nearest', 'mirror', 'wrap', 'constant'):
        cases.append((f'real {mode}', ndi.vectorized_filter(
            x, features, size=(2, 3, 2), mode=mode,
            **({'cval': -1.75} if mode == 'constant' else {}))))
        cases.append((f'complex {mode}', ndi.vectorized_filter(
            z, features, size=(2, 3, 2), mode=mode,
            **({'cval': -1.75 + .6j} if mode == 'constant' else {}))))
    # SciPy 1.17's multi-valued return with nondefault axes permutes its
    # extra channel dimension unexpectedly. Two scalar-result callbacks
    # provide the unambiguous (original input shape, channels) reference.
    def two_scalar_features(data, **kwargs):
        mean = ndi.vectorized_filter(data, np.mean, **kwargs)
        square = ndi.vectorized_filter(
            data, lambda w, *, axis: np.mean(w * w, axis=axis), **kwargs)
        return np.stack((mean, square), axis=-1)

    cases.append(('selected reordered sparse', two_scalar_features(
        x, footprint=sparse, axes=(2, 0), origin=(0, 1), mode='nearest')))
    cases.append(('valid real', ndi.vectorized_filter(x, features,
                                                      size=(2, 3, 2), mode='valid')))
    cases.append(('valid complex', ndi.vectorized_filter(z, features,
                                                         size=(2, 3, 2), mode='valid')))
    cases.append(('valid selected footprint', two_scalar_features(
        x, footprint=sparse, axes=(2, 0), mode='valid')))
    cases.append(('shifted origin', ndi.vectorized_filter(
        x, features, size=(2, 3, 2), mode='mirror', origin=(-1, 1, 0))))
    nearest = ndi.vectorized_filter(x, features, size=(2, 3, 2), mode='nearest')
    cases.append(('flat channels', nearest.reshape((x.size, 2), order='F')))
    cases.append(('strided destination', nearest))
    cases.append(('scalar reference, fast path', ndi.vectorized_filter(
        x, np.mean, size=(2, 3, 2), mode='nearest')))

    with tempfile.TemporaryDirectory(prefix='scifort-vectorized-multi-') as tmp_str:
        tmp = Path(tmp_str)
        sources = ['src/kinds.f90', 'src/signal_spline.f90',
                   'src/signal_spline_nd.f90', 'src/ndimage_morphology.f90',
                   'src/ndimage_distance.f90', 'src/ndimage_chamfer.f90',
                   'src/ndimage_measurements.f90', 'src/ndimage_watershed.f90',
                   'src/ndimage.f90',
                   'tools/ndimage_reference/vectorized_multi_driver.f90']
        cmd = ['gfortran', '-std=f2018', '-O0', '-fcheck=all',
               '-Werror=implicit-interface', '-J', str(tmp), '-I', str(tmp)]
        cmd += [str(ROOT / name) for name in sources] + ['-o', str(tmp / 'driver')]
        subprocess.run(cmd, check=True, cwd=tmp)
        subprocess.run([str(tmp / 'driver')], check=True, cwd=tmp,
                       stdout=subprocess.PIPE, text=True)
        total = 0
        largest = 0.0
        with (tmp / 'vectorized_multi_reference.bin').open('rb') as f:
            for name, expected in cases:
                dtype = np.complex128 if np.iscomplexobj(expected) else np.float64
                length = expected.size * np.dtype(dtype).itemsize
                data = f.read(length)
                if len(data) != length:
                    raise AssertionError(f'{name}: truncated binary output')
                actual = np.frombuffer(data, dtype).reshape(expected.shape, order='F')
                difference = np.abs(actual - expected)
                error = float(np.max(difference))
                if not np.isfinite(error) or error > 1e-11:
                    index = np.unravel_index(np.argmax(difference), expected.shape)
                    raise AssertionError(f'{name}: at {index}, got {actual[index]}, '
                                         f'expected {expected[index]}, error {error}')
                largest = max(largest, error)
                total += expected.size
            if f.read(1):
                raise AssertionError('trailing reference data')
    print(f'PASS: {total} multi-channel SciPy values across {len(cases)} cases')
    print(f'Maximum absolute error: {largest:.3g}')


if __name__ == '__main__':
    main()
