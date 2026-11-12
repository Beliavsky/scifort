#!/usr/bin/env python3
"""Verify mixed-mode Gaussian/uniform and complex neighborhood callbacks."""
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
    'tools/ndimage_reference/mixed_complex_driver.f90')]
MODES=['reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']


def reducer(values):
    return np.dot(values,np.arange(1,values.size+1,dtype=np.float64))


def main():
    i,j=np.mgrid[1:6,1:7]
    a=np.sin(.4*i+.7*j)+(i*i-j)/7.
    z=a+1j*np.cos(.3*i*j)
    weights=np.array([.4,-.6,1.1,.8,-.3,.2]).reshape((2,3),order='F')
    cweights=weights+1j*(.3*weights+.15)
    foot=np.array([1,0,1,1,0,1],bool).reshape((2,3),order='F')
    i,j,k=np.mgrid[1:5,1:4,1:3]
    a3=np.sin(.3*(i+2*j+3*k))
    with tempfile.TemporaryDirectory(prefix='scifort-mixed-complex-') as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                        *map(str,SOURCES),'-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        data=(tmp/'mixed_complex_results.bin').read_bytes()
    offset=0
    count=0
    def check(expected,name):
        nonlocal count,offset
        expected=np.asarray(expected)
        dtype=np.complex128 if np.iscomplexobj(expected) else np.float64
        expected=expected.astype(dtype).ravel(order='F')
        got=np.frombuffer(data,dtype=dtype,count=expected.size,offset=offset)
        if not np.allclose(got,expected,atol=2e-11,rtol=2e-11):
            error=np.max(np.abs(got-expected))
            raise AssertionError((name,error,np.argmax(np.abs(got-expected)),got[:8],expected[:8]))
        offset+=got.nbytes
        count+=expected.size
    for mode1 in MODES:
        for mode2 in MODES:
            modes=(mode1,mode2)
            cval=-.7+.6j
            check(ndi.gaussian_filter(a,(.85,1.37),order=(1,0),mode=modes,cval=-.7,truncate=3),('gauss',modes))
            check(ndi.gaussian_filter(z,(.85,1.37),order=(0,2),mode=modes,cval=cval,radius=(2,3)),('gauss complex',modes))
            check(ndi.uniform_filter(a,(2,3),mode=modes,cval=-.7),('uniform',modes))
            check(ndi.uniform_filter(z,(2,4),mode=modes,cval=cval),('uniform complex',modes))
            check(ndi.generic_filter(a,reducer,footprint=foot,origin=(-1,1),mode=mode1,cval=-.7) if mode1==mode2 else
                  generic_modes_reference(a,foot,(-1,1),modes,-.7),('generic real',modes))
            check(generic_modes_reference(z,foot,(-1,1),modes,cval),('generic complex',modes))
            check(ndi.minimum_filter(a,(2,3),mode=modes,cval=-.7),('minimum',modes))
            check(ndi.maximum_filter(a,(3,2),mode=modes,cval=-.7),('maximum',modes))
            check(full_filter_reference(a,weights,modes,-.7,False),('correlate',modes))
            check(full_filter_reference(z,cweights,modes,cval,True),('convolve',modes))
    for n in range(8):
        modes=(MODES[n],MODES[(n+1)%8],MODES[(n+4)%8])
        check(ndi.gaussian_filter(a3,(.8,.4,1.2),order=(0,1,0),mode=modes,cval=1.3),('gauss 3d',modes))
        check(ndi.uniform_filter(a3,(2,3,4),mode=modes,cval=1.3),('uniform 3d',modes))
    if offset!=len(data):
        raise AssertionError(('unconsumed data',offset,len(data)))
    print(f'PASS: {count} mixed-mode and complex-generic SciPy/reference elements')


def generic_modes_reference(arr,foot,origin,modes,cval):
    """Independent direct-window extension per axis, in C-order samples."""
    dims=arr.shape
    out=np.empty_like(arr)
    offsets=[np.array(t)-np.array(foot.shape)//2-np.array(origin)
             for t in np.argwhere(foot)]
    for loc in np.ndindex(dims):
        vals=[]
        for off in offsets:
            coord=np.array(loc)+off
            good=True
            for d,v in enumerate(coord):
                mapped=map_boundary(int(v),dims[d],modes[d])
                if mapped is None:
                    good=False
                    break
                coord[d]=mapped
            vals.append(arr[tuple(coord)] if good else cval)
        out[loc]=np.dot(np.arange(1,len(vals)+1),vals)
    return out


def full_filter_reference(arr,kernel,modes,cval,convolution):
    """Independent coordinate-space convolution/correlation under axis modes."""
    out=np.empty_like(arr)
    for loc in np.ndindex(arr.shape):
        total=0j if np.iscomplexobj(arr) else 0.
        for kk in np.ndindex(kernel.shape):
            if convolution:
                offsets=np.array(kernel.shape)//2-np.array(kk)
                coeff=kernel[kk]
            else:
                offsets=np.array(kk)-np.array(kernel.shape)//2
                coeff=np.conj(kernel[kk])
            neighbor=np.array(loc)+offsets
            good=True
            for d,v in enumerate(neighbor):
                mapped=map_boundary(int(v),arr.shape[d],modes[d])
                if mapped is None:
                    good=False
                    break
                neighbor[d]=mapped
            total+=coeff*(arr[tuple(neighbor)] if good else cval)
        out[loc]=total
    return out


def map_boundary(i,n,mode):
    if mode in ('constant','grid-constant'):
        return i if 0<=i<n else None
    if mode=='nearest': return max(0,min(n-1,i))
    if mode in ('wrap','grid-wrap'): return i%n
    if mode in ('reflect','grid-mirror'):
        p=i%(2*n)
        return 2*n-1-p if p>=n else p
    if mode=='mirror':
        if n==1: return 0
        p=i%(2*(n-1))
        return 2*(n-1)-p if p>=n else p
    raise ValueError(mode)


if __name__=='__main__':
    main()
