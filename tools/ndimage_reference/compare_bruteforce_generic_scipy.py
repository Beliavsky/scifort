#!/usr/bin/env python3
"""Check SciFort's brute-force distances and arbitrary derivative callbacks."""
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
    'tools/ndimage_reference/bruteforce_generic_driver.f90')]
MODES=['reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']
METRICS=['euclidean','taxicab','chessboard']


def first_derivative(x,axis,output,mode,cval):
    return ndi.correlate1d(x,[-0.5,0,0.5],axis=axis,output=output,mode=mode,cval=cval)


def second_derivative(x,axis,output,mode,cval):
    return ndi.correlate1d(x,[1,-2,1],axis=axis,output=output,mode=mode,cval=cval)


def read_values(stream,shape,dtype):
    n=np.prod(shape)
    data=np.frombuffer(stream.read(n*np.dtype(dtype).itemsize),dtype=dtype)
    if data.size!=n:
        raise AssertionError('incomplete Fortran output stream')
    return data.reshape(shape,order='F')


def check_distances(stream,mask,metric,sampling=None):
    dist=read_values(stream,mask.shape,np.float64)
    coords=read_values(stream,(mask.ndim,mask.size),np.int32)
    kwargs={'metric':metric,'return_indices':True}
    if sampling is not None: kwargs['sampling']=sampling
    expected,nearest=ndi.distance_transform_bf(mask,**kwargs)
    if mask.all():
        if not np.array_equal(dist,expected) or np.any(coords):
            raise AssertionError(('all-foreground',metric,dist,expected,coords))
        return dist.size
    err=np.max(np.abs(dist-expected))
    if err>2.e-12:
        raise AssertionError(('distances',mask.shape,metric,sampling,err,dist,expected))
    if not np.any(mask):
        one_based=(nearest.reshape(mask.ndim,-1,order='F')+1)
        if not np.array_equal(coords,one_based):
            raise AssertionError(('all-background features',coords,one_based))
    for n,position in enumerate(np.ndindex(mask.shape)):
        # Fortran flattens the same coordinates in column-major order.
        flat=np.ravel_multi_index(position,mask.shape,order='F')
        where=tuple(coords[:,flat]-1)
        if any(c<0 or c>=lim for c,lim in zip(where,mask.shape)):
            raise AssertionError(('feature bounds',metric,position,where))
        if mask[where]:
            raise AssertionError(('feature not background',metric,position,where))
        if metric=='euclidean':
            steps=np.array(position)-np.array(where)
            if sampling is not None: steps=steps*np.array(sampling)
            distance=np.linalg.norm(steps)
        elif metric=='taxicab':
            distance=sum(abs(a-b) for a,b in zip(position,where))
        else:
            distance=max(abs(a-b) for a,b in zip(position,where))
        if abs(distance-dist[position])>2.e-12:
            raise AssertionError(('feature wrong distance',metric,position,where))
    return dist.size


def main():
    mask1=np.array([(7*i+2)%9>2 for i in range(1,9)],dtype=bool)
    mask2=np.fromfunction(lambda i,j: (3*(i+1)+5*(j+1)+(i+1)*(j+1))%11>3,
                          (7,6),dtype=int)
    mask3=np.fromfunction(lambda i,j,k: (13*(i+1)+5*(j+1)+3*(k+1)+(i+1)*(k+1))%13>2,
                          (4,5,3),dtype=int)
    mask15=np.ones((2,)+(1,)*14,dtype=bool)
    mask15[(1,)+(0,)*14]=False
    a2=np.fromfunction(lambda i,j: np.sin(.27*(i+1)+.41*(j+1))+\
                       ((i+1)*(j+1)-(j+1)**2)/13.,(5,6))
    a3=np.fromfunction(lambda i,j,k: np.sin(.2*(i+1)+.5*(j+1)+.7*(k+1))+\
                       ((i+1)*(j+1)*(k+1))/11.,(4,3,2))
    total=0
    max_error=0.0
    with tempfile.TemporaryDirectory(prefix='scifort-bruteforce-generic-') as t:
        t=Path(t)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow', '-Werror=implicit-interface',
                        '-J',str(t),'-I',str(t),*(str(p) for p in SOURCES),
                        '-o',str(t/'driver')],cwd=t,check=True)
        subprocess.run([str(t/'driver')],cwd=t,check=True)
        with (t/'bruteforce_generic_results.bin').open('rb') as f:
            for mask in (mask1,mask2,mask3):
                for metric in METRICS:
                    total+=check_distances(f,mask,metric)
            total+=check_distances(f,mask2,'euclidean',[1.2,2.7])
            total+=check_distances(f,mask3,'euclidean',[.7,1.8,2.3])
            for metric in METRICS:
                total+=check_distances(f,np.ones(8,bool),metric)
            for metric in METRICS:
                total+=check_distances(f,np.zeros(8,bool),metric)
            for metric in METRICS:
                total+=check_distances(f,mask15,metric)
            for p,mode in enumerate(MODES):
                mixed=(mode,MODES[(p+1)%8],MODES[(p+4)%8])
                for x,kind,boundary,cval in (
                    (a2,'gradient',mode,-.9),(a2,'laplace',mode,-.9),
                    (a3,'gradient',mixed,.75),(a3,'laplace',mixed,.75),
                ):
                    got=read_values(f,x.shape,np.float64)
                    if kind=='gradient':
                        expected=ndi.generic_gradient_magnitude(x,first_derivative,
                                                                 mode=boundary,cval=cval)
                    else:
                        expected=ndi.generic_laplace(x,second_derivative,
                                                       mode=boundary,cval=cval)
                    err=float(np.max(np.abs(got-expected)))
                    max_error=max(max_error,err)
                    if err>2.e-12:
                        raise AssertionError((kind,x.shape,boundary,err,got,expected))
                    total+=got.size
            if f.read(1): raise AssertionError('unexpected trailing output')
    print(f'PASS: {total} brute-force distance and generic derivative values, max error={max_error:.3g}')


if __name__=='__main__':
    main()
