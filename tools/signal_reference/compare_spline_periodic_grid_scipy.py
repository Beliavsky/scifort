#!/usr/bin/env python3
"""Compare linear-time periodic spline filtering and shaped coordinates to SciPy."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
LENGTHS = (2, 3, 4, 5, 7, 11, 17, 63, 128, 512, 1024)
MODES = ('grid-wrap', 'nearest', 'grid-constant', 'constant', 'mirror')
SOURCES = tuple(ROOT / name for name in (
    'src/kinds.f90', 'src/signal_spline.f90', 'src/signal_spline_nd.f90',
    'tools/signal_reference/spline_periodic_grid_driver.f90',
))


def main():
    count = 0
    maximum = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-spline-periodic-') as folder:
        folder = Path(folder)
        subprocess.run(['gfortran', '-std=f2018', '-O0', '-fcheck=all',
                        '-Werror=implicit-interface', '-I', str(folder),
                        '-J', str(folder), *map(str, SOURCES), '-o',
                        str(folder/'driver')], cwd=folder, check=True)
        subprocess.run([str(folder/'driver')], cwd=folder, check=True)

        def check(label, actual, expected, tol=4e-10):
            nonlocal count, maximum
            if actual.shape != expected.shape:
                raise AssertionError((label, actual.shape, expected.shape))
            err = float(np.max(np.abs(actual-expected)))
            if not np.isfinite(err) or err > tol:
                raise AssertionError((label, err, actual, expected))
            count += actual.size
            maximum = max(maximum, err)

        with (folder/'spline_periodic_grid.bin').open('rb') as f:
            def read(shape, dtype):
                n = int(np.prod(shape))
                b = f.read(n*np.dtype(dtype).itemsize)
                if len(b) != n*np.dtype(dtype).itemsize:
                    raise AssertionError(('short binary stream', shape, dtype))
                return np.frombuffer(b,dtype).reshape(shape,order='F')

            for n in LENGTHS:
                x = np.sin(.37*np.arange(1,n+1))+.013*np.arange(1,n+1)**2
                for order in range(2,6):
                    check(('periodic coefficients',n,order),read((n,),np.float64),
                          ndi.spline_filter1d(x,order=order,mode='grid-wrap'))
                    z=x.astype(complex)+1j*(-.11*x+.07)
                    check(('complex periodic coefficients',n,order),
                          read((n,),np.complex128),
                          ndi.spline_filter1d(z,order=order,mode='grid-wrap',
                                              output=np.complex128))
            i, j = np.mgrid[1:8,1:10]
            sample = np.sin(.21*i+.33*j)+np.cos(.07*i*j)
            z = sample.astype(complex)+1j*(.7*sample-.13)
            c = np.empty((2,5,4),dtype=float)
            i, j = np.mgrid[1:6,1:5]
            c[0] = -.2+.91*i+.43*j
            c[1] = 1.3+.55*i+1.2*j
            for mode in MODES:
                for order in range(6):
                    for image, fill, dtype in ((sample,-1.25,np.float64),
                                               (z,-1.25+.42j,np.complex128)):
                        expected=ndi.map_coordinates(image,c-1,mode=mode,
                            order=order,cval=fill,prefilter=True)
                        check(('grid output',mode,order,dtype),read((5,4),dtype),expected)
            check('prefilter false',read((5,4),np.float64),
                  ndi.map_coordinates(sample,c-1,mode='grid-wrap',order=5,
                                      prefilter=False))
            if f.read(1):
                raise AssertionError('trailing output')
    print(f'PASS: {count} spline periodic/grid SciPy values, max error {maximum:.3g}')


if __name__ == '__main__':
    main()
