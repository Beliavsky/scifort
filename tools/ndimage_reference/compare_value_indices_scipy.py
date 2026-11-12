#!/usr/bin/env python3
"""Compare SciFort one-based CSR value_indices with scipy.ndimage.value_indices."""
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
    'tools/ndimage_reference/value_indices_driver.f90',
)]


def main():
    i = np.mgrid[1:12]
    a1 = (7 * i) % 11 - 5
    i, j = np.mgrid[1:10, 1:8]
    a2 = (13*i + 17*j + 3*i*j) % 31 - 13
    i, j, k = np.mgrid[1:5, 1:6, 1:4]
    a3 = (17*i + 31*j + 7*k + i*j*k) % 13 - 6
    i, j = np.mgrid[1:65, 1:49]
    a4 = (39*i + 17*j + 3*i*j) % 25 - 12
    empty_labels = np.ones((2, 2), dtype=int)
    high = np.full((2,) + (1,)*14, 17, dtype=np.int32)
    high[(0,) + (0,)*14] = -4
    tests = [(a1, None), (a2, 0), (a3, -1), (a3, 55),
             (a4, 0), (empty_labels, 1), (high, 17)]

    with tempfile.TemporaryDirectory(prefix='scifort-value-indices-') as temp:
        temp = Path(temp)
        common = ['gfortran', '-std=f2018', '-O0', '-Werror=implicit-interface',
                  '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                  '-J', str(temp), '-I', str(temp)]
        objects = []
        for source in SOURCES[:-1]:
            obj = temp / (source.stem + '.o')
            subprocess.run([*common, '-c', str(source), '-o', str(obj)],
                           cwd=temp, check=True, capture_output=True, text=True)
            objects.append(str(obj))
        exe = temp / 'value_indices_driver'
        subprocess.run([*common, *objects, str(SOURCES[-1]), '-o', str(exe)],
                       cwd=temp, check=True, capture_output=True, text=True)
        subprocess.run([str(exe)], cwd=temp, check=True, capture_output=True, text=True)
        payload = (temp / 'value_indices_results.bin').read_bytes()

    offset = 0
    total_indices = 0

    def take(count):
        nonlocal offset
        result = np.frombuffer(payload, dtype=np.int32, count=count,
                               offset=offset).copy()
        offset += count * 4
        return result

    for arr, ignore in tests:
        rank, nkeys, npoints = take(3)
        assert rank == arr.ndim
        keys = take(nkeys)
        starts = take(nkeys+1)
        coords = take(rank * npoints).reshape((rank, npoints), order='F')
        python_indices = ndi.value_indices(arr.astype(np.int32), ignore_value=ignore)
        reference_keys = np.array(sorted(python_indices), dtype=np.int32)
        np.testing.assert_array_equal(keys, reference_keys)
        assert starts[0] == 1 and starts[-1] == npoints+1
        assert len(starts) == len(keys)+1
        for j, key in enumerate(keys):
            python_coords = np.array(python_indices[int(key)], dtype=np.int32)+1
            np.testing.assert_array_equal(coords[:, starts[j]-1:starts[j+1]-1],
                                          python_coords, err_msg=f'value {key}')
        total_indices += npoints
    assert offset == len(payload), (offset, len(payload))
    print(f'PASS: {total_indices} SciPy value_indices coordinates, '
          f'{len(tests)} grouped-array cases')


if __name__ == '__main__':
    main()
