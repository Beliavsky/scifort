#!/usr/bin/env python3
"""Check SciPy-compatible axis-specific 1D uniform and extrema filters."""
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
    'tools/ndimage_reference/axis_extrema_driver.f90')]
MODES = ('reflect', 'grid-mirror', 'mirror', 'nearest', 'wrap',
         'grid-wrap', 'constant', 'grid-constant')
WIDTHS = (1, 2, 3, 4, 6)


def read_values(stream, shape, dtype):
    n = int(np.prod(shape))
    data = np.frombuffer(stream.read(n*np.dtype(dtype).itemsize), dtype=dtype)
    if data.size != n:
        raise AssertionError('incomplete Fortran results')
    return data.reshape(shape, order='F')


def main():
    a = np.fromfunction(lambda i,j:
                        np.sin(.4*((i+1)*2+(j+1)))+(i+1)*(j+1)/17.,(5,7))
    z = a + 1j*np.fromfunction(lambda i,j:
                                np.cos(.6*(i+1)+.17*(j+1)),(5,7))
    count = 0
    max_err = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-axis-extrema-') as td:
        td = Path(td)
        cmd = ['gfortran','-std=f2018','-O2','-fcheck=all',
               '-ffpe-trap=invalid,zero,overflow','-Werror=implicit-interface',
               '-J',str(td),'-I',str(td),*(str(src) for src in SOURCES),
               '-o',str(td/'reference')]
        subprocess.run(cmd,cwd=td,check=True)
        subprocess.run([str(td/'reference')],cwd=td,check=True)
        with (td/'axis_extrema_results.bin').open('rb') as stream:
            for mode in MODES:
                for width in WIDTHS:
                    for origin in (-width//2 if False else -(width//2), 0, (width-1)//2):
                        for axis in range(2):
                            for routine in (ndi.uniform_filter1d,
                                            ndi.minimum_filter1d,ndi.maximum_filter1d):
                                got = read_values(stream,a.shape,np.float64)
                                expected = routine(a,width,axis=axis,mode=mode,
                                                   origin=origin,cval=-1.75)
                                err = np.max(np.abs(got-expected))
                                max_err = max(max_err,float(err))
                                if err>3e-13:
                                    raise AssertionError((routine.__name__,mode,width,
                                                          origin,axis,err,got,expected))
                                count += got.size
                            got = read_values(stream,z.shape,np.complex128)
                            expected = ndi.uniform_filter1d(z,width,axis=axis,mode=mode,
                                                            origin=origin,cval=-1.75+.6j)
                            err = np.max(np.abs(got-expected))
                            max_err = max(max_err,float(err))
                            if err>3e-13:
                                raise AssertionError(('complex uniform',mode,width,origin,
                                                      axis,err,got,expected))
                            count += 2*got.size
            if stream.read(1):
                raise AssertionError('extra bytes in Fortran output')
    print(f'PASS: {count} 1D axis-filter real/complex scalar comparisons, '
          f'max error {max_err:.3g}')


if __name__ == '__main__':
    main()
