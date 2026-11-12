#!/usr/bin/env python3
"""Validate selected axes and monotone-deque extrema against SciPy."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]

def main():
    i,j,k=np.indices((4,5,3))
    x=np.sin(.31*((i+1)+2*(j+1)+4*(k+1)))+.3*((i+1)-(j+1))
    with TemporaryDirectory() as d:
        out=Path(d)
        files=['kinds','signal_spline','signal_spline_nd','ndimage_morphology',
               'ndimage_distance','ndimage_chamfer','ndimage_measurements',
               'ndimage_watershed','ndimage']
        cmd=['gfortran','-std=f2018','-O0','-fcheck=all','-Werror=implicit-interface',
             '-J',str(out),'-I',str(out)]+[str(ROOT/'src'/f'{f}.f90') for f in files]+[
             str(ROOT/'tools/ndimage_reference/selected_axes_driver.f90'),'-o',str(out/'driver')]
        subprocess.run(cmd,cwd=out,check=True)
        subprocess.run([str(out/'driver')],cwd=out,check=True)
        cases=[
            ndi.minimum_filter(x,size=(2,4),mode=('wrap','mirror'),axes=(2,0)),
            ndi.maximum_filter(x,size=(2,4),mode=('wrap','mirror'),axes=(2,0)),
            ndi.uniform_filter(x,size=(2,4),mode=('wrap','mirror'),axes=(2,0)),
            ndi.minimum_filter(x,size=4,mode='constant',cval=-3.5,axes=(1,)),
            ndi.maximum_filter(x,size=3,mode='nearest',axes=(0,2)),
            ndi.uniform_filter(x,size=3,mode='reflect',axes=(0,2)),
            ndi.minimum_filter1d(x,size=5000,axis=1,mode='grid-wrap',origin=1),
            ndi.maximum_filter1d(x,size=5000,axis=1,mode='constant',cval=2.5,origin=-2500),
            ndi.maximum_filter(x,size=(8,1,7),mode='nearest'),
            ndi.minimum_filter(x,size=(1,1,1),mode='constant',cval=-77.0),
        ]
        cases.extend([
            ndi.minimum_filter(x,size=(4,3),axes=(2,0),
                               mode=('constant','nearest'),cval=-3.2),
            ndi.maximum_filter(x,size=(4,3),axes=(2,0),
                               mode=('constant','mirror'),cval=3.2),
            ndi.uniform_filter(x,size=(4,3),axes=(2,0),
                               mode=('constant','nearest'),cval=-3.2),
            x,
        ])
        z=x+1j*(.3*x+.4)
        complex_cases=[
            ndi.uniform_filter(z,size=(2,3),axes=(2,0),
                               mode=('constant','wrap'),cval=1.25-.5j),
            z,
        ]
        zdata=np.fromfile(out/'selected_axes_complex_output.bin',dtype=np.complex128)
        assert zdata.size==len(complex_cases)*x.size
        for c,want in enumerate(complex_cases):
            got=zdata[c*x.size:(c+1)*x.size].reshape(x.shape,order='F')
            err=np.max(np.abs(got-want))
            if err>2e-12:raise AssertionError(('complex',c,err))
        data=np.fromfile(out/'selected_axes_output.bin',dtype=np.float64)
        assert data.size==len(cases)*x.size,(data.size,len(cases)*x.size)
        worst=0.0
        for c,want in enumerate(cases):
            got=data[c*x.size:(c+1)*x.size].reshape(x.shape,order='F')
            err=np.max(np.abs(got-want))
            if err>2e-12: raise AssertionError((c,err,got,want))
            worst=max(worst,float(err))
    print(f'PASS: {len(cases)+len(complex_cases)} selected-axis/large-window SciPy comparisons, max error {worst:.3g}')

if __name__=='__main__':main()
