#!/usr/bin/env python3
"""SciPy comparisons for rank-2 and rank-3 labeled image statistics."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]


def main():
    i,j=np.mgrid[1:6,1:7]
    image2=np.sin(i+2*j)+(i*i+3*j)%17
    labels2=(i+2*j)%5-1
    i,j,k=np.mgrid[1:4,1:3,1:5]
    image3=np.cos(i-j*3+k)+(i*j+5*k)%13
    labels3=(i+3*j+2*k)%5-1
    ids=np.array([-1,0,1,2,3,99,1])
    sources=[ROOT/z for z in ('src/kinds.f90','src/signal_spline.f90',
        'src/signal_spline_nd.f90','src/ndimage_morphology.f90',
        'src/ndimage_distance.f90','src/ndimage_chamfer.f90','src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
        'tools/ndimage_reference/measurements_driver.f90')]
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O0','-fcheck=all',
                '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                *map(str,sources),'-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        values=np.fromfile(tmp/'measurements_results.bin',np.float64)
        count=0
        errors=[]
        cursor=0
        for x,labels in ((image2,labels2),(image3,labels3)):
            with np.errstate(divide='ignore',invalid='ignore'):
                for function in (ndi.sum_labels,ndi.mean,ndi.variance,ndi.center_of_mass):
                    expected=np.array(function(x,labels,ids),dtype=float)
                    if function==ndi.center_of_mass:
                        expected=expected.T+1.0 # Fortran centers are one-based (rank, selected label)
                    flat=np.asarray(expected).ravel(order='F')
                    got=values[cursor:cursor+flat.size]
                    assert got.size==flat.size
                    for a,b in zip(flat,got):
                        if np.isnan(a):
                            assert np.isnan(b),(function.__name__,a,b)
                        else:
                            assert abs(a-b)<1e-12,(function.__name__,a,b)
                    cursor+=flat.size
                    count+=flat.size
        assert cursor==values.size,(cursor,values.size)
        print(f'PASS: labeled measurements {count} SciPy scalar comparisons')


if __name__=='__main__': main()
