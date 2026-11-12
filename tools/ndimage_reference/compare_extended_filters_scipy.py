#!/usr/bin/env python3
"""Compare asymmetric ndimage kernels and neighborhood percentiles to SciPy."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
MODES = ('reflect', 'mirror', 'nearest', 'wrap', 'constant',
         'grid-constant', 'grid-mirror', 'grid-wrap')


def main():
    x = np.empty((4, 5), order='F')
    z = np.empty((4, 5), dtype=np.complex128, order='F')
    w = np.empty((2, 3), order='F')
    cw = np.empty((2, 3), dtype=np.complex128, order='F')
    v = np.empty((4, 2), order='F')
    a = np.empty((4, 3, 2), order='F')
    w3 = np.empty((2, 3, 2), order='F')
    for j in range(1, 6):
        for i in range(1, 5):
            x[i-1, j-1] = np.sin(.3*i+.7*j)+(17*i+13*j)%9
            z[i-1, j-1] = x[i-1, j-1]+1j*(-.37*x[i-1, j-1]+.11*i)
    for j in range(1, 4):
        for i in range(1, 3):
            w[i-1, j-1] = np.sin(.6*i-.29*j)+.1*i*j
            cw[i-1, j-1] = w[i-1, j-1]+1j*(.16*j-.23*i)
    for j in range(1, 3):
        for i in range(1, 5):
            v[i-1, j-1] = np.cos(.21*i+.44*j)+.13*i*j
    for k in range(1, 3):
        for j in range(1, 4):
            for i in range(1, 5):
                a[i-1, j-1, k-1] = np.sin(.52*i-.19*j+.11*k)+(13*i+7*j+11*k)%7
    for k in range(1, 3):
        for j in range(1, 4):
            for i in range(1, 3):
                w3[i-1, j-1, k-1] = np.sin(.22*i+.17*j-.28*k)
    count = 0
    largest = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-extended-filter-') as tmp:
        tmp = Path(tmp)
        subprocess.run(['gfortran', '-std=f2018', '-O2', '-fcheck=all',
                        '-Werror=implicit-interface', '-J', str(tmp), '-I', str(tmp),
                        *(str(ROOT/f) for f in ('src/kinds.f90', 'src/signal_spline.f90',
                         'src/signal_spline_nd.f90', 'src/ndimage_morphology.f90',
                         'src/ndimage_distance.f90','src/ndimage_chamfer.f90', 'src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
                         'tools/ndimage_reference/extended_filter_driver.f90')),
                        '-o', str(tmp/'driver')], cwd=tmp, check=True)
        subprocess.run([str(tmp/'driver')], cwd=tmp, check=True)
        with (tmp/'extended_filter_results.bin').open('rb') as f:
            for mode in MODES:
                kw = dict(mode=mode, cval=-2.3)
                cases = [
                    ('correlate', ndi.correlate(x, w, **kw)),
                    ('convolve', ndi.convolve(x, w, **kw)),
                    ('correlate_even', ndi.correlate(x, v, **kw)),
                    ('convolve_even', ndi.convolve(x, v, **kw)),
                    ('median_vector', ndi.median_filter(x, size=[2,3], **kw)),
                    ('rank0', ndi.rank_filter(x, 0, size=[3,2], **kw)),
                    ('rank_negative', ndi.rank_filter(x, -2, size=[2,3], **kw)),
                    ('percentile25', ndi.percentile_filter(x, 25, size=[2,3], **kw)),
                    ('percentile75', ndi.percentile_filter(x, 75, size=[3,2], **kw)),
                    ('median_scalar', ndi.median_filter(x, size=3, **kw)),
                    ('prewitt2', ndi.prewitt(x,axis=1,**kw)),
                    ('complex_correlate', ndi.correlate(z, cw, mode=mode,cval=-2.3+.7j)),
                    ('complex_convolve', ndi.convolve(z, cw, mode=mode,cval=-2.3+.7j)),
                    ('correlate3', ndi.correlate(a, w3, **kw)),
                    ('convolve3', ndi.convolve(a, w3, **kw)),
                    ('rank3', ndi.rank_filter(a, 4, size=[3,2,2], **kw)),
                    ('percentile3', ndi.percentile_filter(a,95,size=[2,3,2],**kw)),
                    ('median3', ndi.median_filter(a,size=[3,1,2],**kw)),
                    ('prewitt3', ndi.prewitt(a,axis=2,**kw)),
                ]
                for name, expected in cases:
                    dtype = np.complex128 if np.iscomplexobj(expected) else np.float64
                    result = np.frombuffer(f.read(expected.size*np.dtype(dtype).itemsize),
                                           dtype=dtype).reshape(expected.shape, order='F')
                    error = float(np.max(np.abs(result-expected)))
                    largest = max(largest, error)
                    if error > 5e-11:
                        ix = np.unravel_index(np.argmax(np.abs(result-expected)),expected.shape)
                        raise AssertionError((mode,name,error,ix,result[ix],expected[ix]))
                    count += 1
            if f.read(1):
                raise AssertionError('trailing bytes')
    print(f'PASS: {count} full-kernel/order-statistic SciPy cases, maximum error={largest:.3g}')


if __name__ == '__main__':
    main()
