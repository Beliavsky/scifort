#!/usr/bin/env python3
"""Numerically compare shifted, non-flat morphological stencils with SciPy."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[2]
MODES = ('reflect','mirror','nearest','wrap','constant','grid-constant',
         'grid-mirror','grid-wrap')
METHODS = (ndi.grey_dilation, ndi.grey_erosion, ndi.grey_opening,
           ndi.grey_closing, ndi.morphological_gradient, ndi.morphological_laplace,
           ndi.white_tophat, ndi.black_tophat)


def main():
    i,j=np.mgrid[1:6,1:7]
    image=np.sin(i*0.31)+(i*i+3*j+i*j)%11
    heights=np.array([0.2,-0.7,1.3,0.4,-1.6,2.2]).reshape((2,3),order='F')
    footprint=np.array([1,0,1,1,0,1],dtype=bool).reshape((2,3),order='F')
    sources = [ROOT/x for x in ('src/kinds.f90','src/signal_spline.f90',
         'src/signal_spline_nd.f90','src/ndimage_morphology.f90',
         'src/ndimage_distance.f90','src/ndimage_chamfer.f90','src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
         'tools/ndimage_reference/advanced_gray_driver.f90')]
    with tempfile.TemporaryDirectory(prefix='scifort-advanced-gray-') as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                        *map(str,sources),'-o',str(tmp/'driver')], check=True, cwd=tmp)
        subprocess.run([str(tmp/'driver')], check=True, cwd=tmp)
        raw=np.fromfile(tmp/'advanced_gray_results.bin',dtype=np.float64)
        count=0
        maxerr=0.
        for mode in MODES:
            for origin in [(0,0),(-1,1),(0,-1),(-1,-1)]:
                for method in METHODS:
                    expected=method(image,footprint=footprint,structure=heights,
                                    origin=origin,mode=mode,cval=-2.5)
                    actual=raw[count*30:(count+1)*30].reshape((5,6),order='F')
                    err=float(np.max(np.abs(actual-expected)))
                    if err>2.e-13:
                        raise AssertionError((method.__name__, mode,origin,err,actual,expected))
                    maxerr=max(maxerr,err)
                    count+=1
        assert count*30==raw.size, (count*30,raw.size)
        print(f'PASS: advanced grayscale {count} SciPy comparisons, max error {maxerr:.3g}')


if __name__=='__main__':
    main()
