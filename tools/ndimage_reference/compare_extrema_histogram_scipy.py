#!/usr/bin/env python3
"""Compare footprint extrema and region counts/histograms/standard deviations."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]
MODES=['reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']
SOURCES=[ROOT/p for p in (
    'src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90','src/ndimage_distance.f90','src/ndimage_chamfer.f90',
    'src/ndimage_measurements.f90','src/ndimage_watershed.f90','src/ndimage.f90',
    'tools/ndimage_reference/extrema_histogram_driver.f90')]

def main():
    i,j=np.mgrid[1:6,1:7]
    image2=np.sin(i*2+j)+(i*i+3*j)%13
    labels2=(i+2*j)%6-2
    i,j,k=np.mgrid[1:5,1:4,1:3]
    image3=np.cos(2*i-j+3*k)+(i*i+2*j+4*k)%13
    labels3=(i+2*j+3*k)%6-2
    ids=np.array([-2,0,1,2,9,1])
    feet=[np.array([1,0,1,1,1,0],bool).reshape((2,3),order='F'),
          np.array([1,0,1,0,1,1,1,1,0],bool).reshape((3,3),order='F')]
    foot3=np.array([1,0,1,0,1,1,0,1,1,0,1,1],bool).reshape((3,2,2),order='F')
    with tempfile.TemporaryDirectory() as work:
        work=Path(work)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all','-Werror=implicit-interface',
            '-J',str(work),'-I',str(work),*[str(s) for s in SOURCES],'-o',str(work/'driver')],
            check=True,cwd=work)
        subprocess.run([str(work/'driver')],check=True,cwd=work)
        data=np.fromfile(work/'extrema_histogram_results.bin',dtype=np.float64)
    offset=0
    count=0
    def check(x,expected,tag):
        nonlocal offset,count
        expected=np.asarray(expected,dtype=np.float64).ravel(order='F')
        actual=data[offset:offset+len(expected)]
        if len(actual)!=len(expected):
            raise AssertionError((tag,'truncated output'))
        if not np.allclose(actual,expected,rtol=0,atol=2e-13,equal_nan=True):
            raise AssertionError((tag,np.max(np.abs(actual-expected)),actual,expected))
        offset+=len(expected)
        count+=len(expected)
    for foot in feet:
        for mode in MODES:
            for origin in [(0,0),(-1,-1),(0,1)]:
                for f in (ndi.minimum_filter,ndi.maximum_filter):
                    check(image2,f(image2,footprint=foot,origin=origin,cval=-3.25,mode=mode),
                          (f.__name__,'2D',mode,origin,foot.shape))
    for mode in MODES:
        for origin in [(0,0,0),(-1,-1,0),(1,0,-1)]:
            for f in (ndi.minimum_filter,ndi.maximum_filter):
                check(image3,f(image3,footprint=foot3,origin=origin,cval=-3.25,mode=mode),
                      (f.__name__,'3D',mode,origin))
    for mode in MODES:
        for f in (ndi.minimum_filter,ndi.maximum_filter):
            check(image2,f(image2,size=(2,3),origin=(-1,1),mode=mode,cval=-3.25),
                  (f.__name__,'rectangle origin',mode))
    for values,labels in [(image2,labels2),(image3,labels3)]:
        check(values,[np.count_nonzero(labels==z) for z in ids],('counts',values.shape))
        histograms=[(np.histogram(values[labels==z],bins=7,range=(-2,14))[0]
                     if np.any(labels==z) else np.zeros(7)) for z in ids]
        check(values,np.stack(histograms,axis=1),('histograms',values.shape))
        with np.errstate(divide='ignore',invalid='ignore'):
            check(values,ndi.standard_deviation(values,labels,ids),('std',values.shape))
    assert offset==len(data),(offset,len(data))
    print(f'PASS: {count} SciPy extrema-filter and region histogram scalar comparisons')

if __name__=='__main__':
    main()
