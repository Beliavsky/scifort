#!/usr/bin/env python3
"""Cross-check Fourier-domain filters with SciPy complex FFT conventions."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]
SOURCES=[ROOT/p for p in (
    'src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90','src/ndimage_distance.f90','src/ndimage_chamfer.f90',
    'src/ndimage_measurements.f90','src/ndimage_watershed.f90','src/ndimage.f90',
    'tools/ndimage_reference/fourier_driver.f90')]


def main():
    i,j=np.mgrid[1:6,1:7]
    a=np.sin(.4*i+.23*j)+1j*np.cos(.27*i-.49*j)
    i,j,k=np.mgrid[1:5,1:4,1:3]
    c=np.sin(.4*i+.7*j-.2*k)+1j*np.sin(.31*i*j*k)
    i,j=np.mgrid[1:6,1:5]
    r=np.sin(.3*i+j)+1j*np.cos(.44*i-j)
    i,j,k=np.mgrid[1:3,1:4,1:6]
    t=np.sin(.5*i+.6*j+.3*k)+1j*np.cos(.41*i-.53*j*k)
    with tempfile.TemporaryDirectory(prefix='scifort-fourier-') as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                        *map(str,SOURCES),'-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        data=(tmp/'fourier_results.bin').read_bytes()
    offset=0
    checked=0
    def check(expected,name):
        nonlocal offset,checked
        flat=np.asarray(expected,dtype=np.complex128).ravel(order='F')
        got=np.frombuffer(data,dtype=np.complex128,count=flat.size,offset=offset)
        if not np.allclose(got,flat,rtol=1e-12,atol=2e-13):
            raise AssertionError((name,np.max(np.abs(got-flat)),got[:8],flat[:8]))
        offset+=flat.nbytes
        checked+=flat.size
    for func,vec,scalar in ((ndi.fourier_gaussian,(1.2,.8),1.1),
                            (ndi.fourier_uniform,(2.1,3.3),2.4),
                            (ndi.fourier_shift,(.75,-1.35),.3)):
        check(func(a,vec),(func.__name__,'2d vector'))
        check(func(a,scalar),(func.__name__,'2d scalar'))
    for func,param in ((ndi.fourier_gaussian,(.8,1.9,.35)),
                       (ndi.fourier_uniform,(3.2,1.6,2.7)),
                       (ndi.fourier_shift,(.4,-.75,1.9))):
        check(func(c,param),(func.__name__,'3d'))
    for func,param in ((ndi.fourier_gaussian,(1.2,.7)),
                       (ndi.fourier_uniform,(2.3,3.1)),
                       (ndi.fourier_shift,(.9,-1.4))):
        check(func(r,param,n=8,axis=0),(func.__name__,'rfft axis 0'))
    for func,param in ((ndi.fourier_gaussian,(1.6,.4,2.1)),
                       (ndi.fourier_uniform,(2.1,3.1,.8)),
                       (ndi.fourier_shift,(1.25,-.4,2.2))):
        check(func(t,param,n=9,axis=-1),(func.__name__,'rfft last axis'))
    if len(data)!=offset:
        raise AssertionError((len(data),offset))
    print(f'PASS: {checked} Fourier-domain complex values agree with SciPy')


if __name__=='__main__': main()
