#!/usr/bin/env python3
"""Check complex generic_filter1d extension and selected binary axes against SciPy 1.17.

SciPy generic_filter1d disallows complex input, so the linear complex callback
is decomposed into four SciPy real-valued calls as an independent reference.
"""
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
    'tools/ndimage_reference/complex_line_binary_axes_driver.f90',
)]
MODES = ['reflect', 'grid-mirror', 'mirror', 'nearest', 'wrap',
         'grid-wrap', 'constant', 'grid-constant']


def complex_scipy_reference(x, width, axis, mode, cval, origin):
    k = np.arange(1, width + 1, dtype=float)
    real_weights = .2 * k
    imag_weights = .03 * np.where(k % 2, -1., 1.) * k

    def make_reducer(weights):
        def reducer(padded, output):
            for j in range(len(output)):
                output[j] = np.dot(padded[j:j+width], weights)
        return reducer

    def real_pass(xpart, weights, exterior):
        return ndi.generic_filter1d(xpart, make_reducer(weights), width,
                                    axis=axis, mode=mode, cval=exterior,
                                    origin=origin)

    a = real_pass(x.real, real_weights, cval.real)
    b = real_pass(x.imag, real_weights, cval.imag)
    c = real_pass(x.real, imag_weights, cval.real)
    d = real_pass(x.imag, imag_weights, cval.imag)
    return (a-d) + 1j*(b+c)


def main():
    i, j = np.mgrid[1:6, 1:7]
    a2 = np.sin(.37*i+.13*j) + 1j*np.cos(.41*j-.19*i)
    i, j, k = np.mgrid[1:4, 1:5, 1:3]
    a3 = np.sin(.14*i*j+.7*k) + 1j*np.cos(.21*i-.31*j*k)
    i, j, k = np.mgrid[1:5, 1:6, 1:4]
    m3 = (3*i+5*j+7*k) % 11 < 4
    mask3 = (i+2*j+3*k) % 5 > 0
    i, j, k, l = np.mgrid[1:3, 1:4, 1:5, 1:3]
    m4 = (i+2*j+3*k+5*l) % 7 < 3
    mask4 = (2*i+j+k+3*l) % 4 != 0
    s23 = np.array([1, 0, 1, 1, 1, 0], dtype=bool).reshape((2, 3), order='F')
    s33 = np.array([1, 1, 0, 0, 1, 1, 1, 0, 1], dtype=bool).reshape((3, 3), order="F")
    s4 = np.array([1, 1, 0, 1], dtype=bool)
    s323 = (np.arange(1, 19) % 4 != 0).reshape((3, 2, 3), order='F')
    morphology = [ndi.binary_dilation, ndi.binary_erosion,
                  ndi.binary_opening, ndi.binary_closing]
    expected = []
    for mode in MODES:
        for x, width, axis, origin in [
                (a2, 4, 0, 1), (a2, 3, 1, -1),
                (a3, 2, 2, 0), (a3, 5, 1, 1)]:
            expected.append(complex_scipy_reference(x, width, axis, mode,
                                                    1.5-.75j, origin))

    options = [
        dict(axes=(0, 2), structure=s23, origin=(-1, 1),
             iterations=2, mask=mask3),
        dict(axes=(2, 0), structure=s33, origin=(1, -1),
             iterations=1, border_value=1),
        dict(axes=(1,), structure=s4, origin=(1,),
             iterations=1, mask=mask3),
        dict(axes=(0, 1, 2), structure=s323, origin=(1, 0, -1),
             iterations=1, border_value=1),
        dict(axes=()),
        dict(axes=(1, 2), structure=ndi.generate_binary_structure(2, 1),
             origin=(0, 1), iterations=2),
        dict(axes=(0, 1), structure=ndi.generate_binary_structure(2, 2),
             iterations=0),
        dict(axes=(0, 1, 2), iterations=0, mask=mask3),
        dict(axes=(0, 2), origin=(-1, 0)),
    ]
    for op in morphology:
        for kwargs in options:
            expected.append(op(m3, **kwargs))
        expected.append(op(m4, structure=s23, axes=(0, 3),
                           origin=(0, 0), iterations=2, mask=mask4))

    with tempfile.TemporaryDirectory(prefix='scifort-callback-binary-') as temp:
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
        exe = temp / 'complex_line_binary_axes_driver'
        subprocess.run([*common, *objects, str(SOURCES[-1]), '-o', str(exe)],
                       cwd=temp, check=True, capture_output=True, text=True)
        subprocess.run([str(exe)], cwd=temp, check=True, capture_output=True, text=True)
        payload = (temp / 'complex_line_binary_axes_results.bin').read_bytes()

    offset = 0
    total_values = 0
    for case_id, value in enumerate(expected):
        dtype = np.complex128 if np.iscomplexobj(value) else np.int32
        count = value.size
        length = count * np.dtype(dtype).itemsize
        actual = np.frombuffer(payload, dtype=dtype, count=count, offset=offset)
        offset += length
        actual = actual.reshape(value.shape, order='F')
        if dtype == np.complex128:
            np.testing.assert_allclose(actual, value, rtol=2e-13, atol=2e-13,
                                       err_msg=f'complex filter case {case_id}')
        else:
            np.testing.assert_array_equal(actual != 0, value,
                                          err_msg=f'axis binary case {case_id}')
        total_values += count
    assert offset == len(payload), (offset, len(payload))
    print(f'PASS: {total_values} complex line-filter and binary-axes reference values')


if __name__ == '__main__':
    main()
