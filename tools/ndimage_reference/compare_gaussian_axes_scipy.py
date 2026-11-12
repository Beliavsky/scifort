#!/usr/bin/env python3
"""Check selected Gaussian axis order, derivatives and complex outputs against SciPy."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]

def main():
    i,j,k=np.indices((4,5,3))
    x=np.sin(.31*((i+1)+2*(j+1)+4*(k+1)))+.3*((i+1)-(j+1))
    z=x+1j*(.73*x-.2)
    with TemporaryDirectory() as d:
        out=Path(d)
        files=['kinds','signal_spline','signal_spline_nd','ndimage_morphology',
               'ndimage_distance','ndimage_chamfer','ndimage_measurements',
               'ndimage_watershed','ndimage']
        cmd=['gfortran','-std=f2018','-O2','-fcheck=all','-Werror=implicit-interface',
             '-J',str(out),'-I',str(out)]+[str(ROOT/'src'/f'{f}.f90') for f in files]+[
             str(ROOT/'tools/ndimage_reference/gaussian_axes_driver.f90'),'-o',str(out/'driver')]
        subprocess.run(cmd,cwd=out,check=True)
        subprocess.run([str(out/'driver')],cwd=out,check=True)
        cases=[
            ndi.gaussian_filter(x,sigma=(.7,1.2),order=(1,0),radius=(3,2),
                                mode=('constant','mirror'),cval=-1.3,axes=(2,0)),
            ndi.gaussian_filter(x,sigma=1.1,order=2,mode='nearest',axes=(1,)),
            ndi.gaussian_filter(x,sigma=(1.0,.4),mode=('wrap','reflect'),axes=(2,0)),
            x,
            ndi.gaussian_filter(z,sigma=(.8,1.1),order=(1,2),
                                mode=('constant','nearest'),cval=.4-.8j,axes=(2,0)),
            ndi.gaussian_filter(z,sigma=1.2,mode='mirror',axes=(0,2))]
        with open(out/'gaussian_axes_output.bin','rb') as f:
            largest=0.0
            for c,want in enumerate(cases):
                dtype=np.complex128 if np.iscomplexobj(want) else np.float64
                got=np.frombuffer(f.read(x.size*np.dtype(dtype).itemsize),dtype=dtype).reshape(x.shape,order='F')
                err=np.max(np.abs(got-want))
                if err>5e-12:
                    at=np.unravel_index(np.argmax(np.abs(got-want)),x.shape)
                    raise AssertionError((c,err,at,got[at],want[at]))
                largest=max(float(err),largest)
            assert not f.read(1),'trailing data'
    print(f'PASS: {len(cases)} Gaussian selected-axis SciPy arrays, max error {largest:.3g}')

if __name__=='__main__':main()
