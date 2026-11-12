#!/usr/bin/env python3
"""Verify binary ndimage morphology/label connectivity against SciPy."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]

def canonical(labels):
    flat=labels.ravel(order='F')
    renumber={0:0}
    result=[]
    for value in flat:
        value=int(value)
        if value not in renumber:
            renumber[value]=len(renumber)
        result.append(renumber[value])
    return np.asarray(result,dtype=np.int32)


def main():
    m2=np.fromfunction(lambda i,j: (7*(i+1)+11*(j+1)+(i+1)*(j+1))%8<3,
                       (5,6),dtype=int)
    m3=np.fromfunction(lambda i,j,k:
        (3*(i+1)+7*(j+1)+5*(k+1)+(i+1)*(j+1)+(k+1)*(j+1))%9<4,
        (4,5,3),dtype=int)
    checks=0
    with tempfile.TemporaryDirectory(prefix='scifort-morphology-') as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                        *map(str,[ROOT/'src/ndimage_morphology.f90',
                                 ROOT/'tools/ndimage_reference/morphology_reference_driver.f90']),
                        '-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        with (tmp/'morphology_reference_results.bin').open('rb') as f:
            for data in (m2,m3):
                for conn in range(1,data.ndim+1):
                    st=ndi.generate_binary_structure(data.ndim,conn)
                    for iterations in (1,2):
                        for border in (0,1):
                            for name,func in (('dilation',ndi.binary_dilation),
                                              ('erosion',ndi.binary_erosion)):
                                want=func(data,structure=st,iterations=iterations,
                                          border_value=border)
                                got=np.frombuffer(f.read(4*data.size),np.int32).reshape(
                                    data.shape,order='F').astype(bool)
                                if not np.array_equal(got,want):
                                    raise AssertionError((data.shape,conn,iterations,border,
                                                          name,np.argwhere(got!=want)))
                                checks+=1
                    want,nref=ndi.label(data,structure=st)
                    n=int(np.frombuffer(f.read(4),np.int32)[0])
                    got=np.frombuffer(f.read(4*data.size),np.int32).reshape(data.shape,order='F')
                    if n!=nref or not np.array_equal(canonical(got),canonical(want)):
                        raise AssertionError(('labels',data.shape,conn,n,nref,got,want))
                    checks+=1
            if f.read(1): raise AssertionError('extra bytes')
    print(f'PASS: {checks} SciPy morphology/connected-component comparisons')

if __name__=='__main__':
    main()
