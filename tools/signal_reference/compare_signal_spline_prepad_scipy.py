#!/usr/bin/env python3
"""Validate padded high-order spline interpolation with singleton dimensions."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]

def main():
    x=np.empty((3,1,4),order='F')
    for k in range(1,5):
        for i in range(1,4):
            x[i-1,0,k-1]=np.sin(.3*i)+np.cos(.21*k)+.15*i*k
    z=x+1j*(.3*x+.2)
    coords=np.array([[1+.5*((k*7)%17-4),1+.4*((k*3)%9-4),
                     1+.6*((k*5)%19-6)] for k in range(1,32)],dtype=float).T
    coords[:,0]=[1,1,1]
    coords[:,1]=[3,1,4]
    coords[:,2]=[3,.5,4]
    coords[:,3]=[1,1.5,1]
    errors=[]
    with tempfile.TemporaryDirectory(prefix='scifort-spline-prepad-') as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                        *map(str,[ROOT/'src/kinds.f90',ROOT/'src/signal_spline.f90',
                                 ROOT/'src/signal_spline_nd.f90',
                                 ROOT/'tools/signal_reference/spline_prepad_reference_driver.f90']),
                        '-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        with (tmp/'spline_prepad_results.bin').open('rb') as f:
            for mode in ('nearest','grid-constant'):
                for order in range(2,6):
                    for arr,fill in ((x,-2.0),(z,-2.0+.75j)):
                        dtype=np.complex128 if np.iscomplexobj(arr) else np.float64
                        got=np.frombuffer(f.read(31*np.dtype(dtype).itemsize),dtype)
                        ref=ndi.map_coordinates(arr,coords-1,order=order,mode=mode,
                                                cval=fill,prefilter=True)
                        err=np.max(abs(got-ref))
                        errors.append(float(err))
                        if err>4e-10:
                            k=np.argmax(abs(got-ref))
                            raise AssertionError((mode,order,arr.dtype,err,k,got[k],ref[k]))
            if f.read(1): raise AssertionError('trailing data')
    print('PASS: 16 high-order singleton-axis spline comparisons,',
          'max difference',max(errors))

if __name__=='__main__':
    main()
