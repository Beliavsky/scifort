#!/usr/bin/env python3
"""Compare label extrema, editable masks and custom rank windows to SciPy."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]
SOURCES=[ROOT/z for z in (
    'src/kinds.f90', 'src/signal_spline.f90', 'src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90', 'src/ndimage_distance.f90','src/ndimage_chamfer.f90',
    'src/ndimage_measurements.f90', 'src/ndimage_watershed.f90', 'src/ndimage.f90',
    'tools/ndimage_reference/regions_masks_rank_driver.f90')]


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O0','-fcheck=all',
            '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
            *map(str,SOURCES),'-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        data=np.fromfile(tmp/'regions_masks_rank_results.bin',dtype=np.float64)
    i,j=np.mgrid[1:6,1:7]
    mask=(3*i+7*j)%5<2
    editable=(i+j)%4!=0
    image=np.sin(i+2*j)+(i*i+3*j)%17
    labels=(i+2*j)%6-2
    bfoot=np.array([1,0,1,1,0,1],dtype=bool).reshape((2,3),order='F')
    foot_a=np.array([1,0,1,1,1,0],dtype=bool).reshape((2,3),order='F')
    foot_b=np.array([1,0,1,0,1,1,1,1,0],dtype=bool).reshape((3,3),order='F')
    origins=[(-1,0),(0,-1),(0,1)]
    modes=['reflect','grid-mirror','mirror','nearest',
           'wrap','grid-wrap','constant','grid-constant']
    ids=[-2,0,1,2,3,9,1]
    cursor=0
    count=0
    def check(expected,name):
        nonlocal cursor,count
        arr=np.asarray(expected,dtype=np.float64).ravel(order='F')
        got=data[cursor:cursor+arr.size]
        assert len(got)==len(arr),(name,len(got),len(arr))
        if not np.allclose(got,arr,rtol=0,atol=1.e-12,equal_nan=True):
            bad=np.where(~np.isclose(got,arr,rtol=0,atol=1.e-12,equal_nan=True))[0]
            k=int(bad[0])
            raise AssertionError(f'{name} element {k}: got {got[k]}, expected {arr[k]}')
        cursor+=len(arr)
        count+=len(arr)
    for operation in ('binary_dilation','binary_erosion','binary_opening','binary_closing'):
        for iterations in (1,2,0,-1):
            for origin in origins:
                expected=getattr(ndi,operation)(mask,structure=bfoot,iterations=iterations,
                                                origin=origin,mask=editable)
                check(expected,(operation,iterations,origin))
    for foot in (foot_a,foot_b):
        for mode in modes:
            for origin in origins:
                args=dict(footprint=foot,origin=origin,mode=mode,cval=-3.0)
                check(ndi.median_filter(image,**args),('median',foot.shape,mode,origin))
                check(ndi.rank_filter(image,-2,**args),('rank',foot.shape,mode,origin))
                check(ndi.percentile_filter(image,63.0,**args),('percentile',foot.shape,mode,origin))
    mins,maxs,pmin,pmax=ndi.extrema(image,labels,ids)
    check(mins,'minimum labels')
    check(maxs,'maximum labels')
    # Convert zero-based C-order positions to one-based Fortran positions;
    # missing-label positions are encoded as (0,...,0).
    presence=[np.any(labels==k) for k in ids]
    check(np.array([(np.array(p)+1 if exists else np.zeros(image.ndim))
                    for p,exists in zip(pmin,presence)]).T,'min positions')
    check(np.array([(np.array(p)+1 if exists else np.zeros(image.ndim))
                    for p,exists in zip(pmax,presence)]).T,'max positions')
    bounds_lo=[]
    bounds_hi=[]
    for k in ids:
        coords=np.argwhere(labels==k)
        if len(coords):
            bounds_lo.append(coords.min(axis=0)+1)
            bounds_hi.append(coords.max(axis=0)+1)
        else:
            bounds_lo.append([0,0]);bounds_hi.append([0,0])
    check(np.array(bounds_lo).T,'box lower')
    check(np.array(bounds_hi).T,'box upper')
    slices=ndi.find_objects(labels,max_label=5)
    check(np.array([[0,0] if bbox is None else
                    [v.start+1 for v in bbox] for bbox in slices]).T,'objects lower')
    check(np.array([[0,0] if bbox is None else
                    [v.stop for v in bbox] for bbox in slices]).T,'objects upper')
    i,j,k=np.mgrid[1:5,1:4,1:3]
    mask3=(7*i+3*j+2*k)%6<3
    editable3=(2*i+j+4*k)%5!=1
    image3=np.sin(3*i-j+2*k)+(i*i+2*j+5*k)%13
    labels3=(i+2*j+3*k)%6-2
    ids3=[-2,0,1,2,5,9]
    foot3=np.array([1,0,1,0,1,1,0,1,1,0,1,1],dtype=bool).reshape((3,2,2),order='F')
    origins3=[(0,0,0),(-1,-1,0),(1,0,-1)]
    for operation in ('binary_dilation','binary_erosion','binary_opening','binary_closing'):
        for iterations in (1,2,0):
            for origin in origins3:
                expected=getattr(ndi,operation)(mask3,structure=foot3,iterations=iterations,
                                                origin=origin,mask=editable3)
                check(expected,('3D',operation,iterations,origin))
    for mode in modes:
        for origin in origins3:
            args=dict(footprint=foot3,origin=origin,mode=mode,cval=-3.0)
            check(ndi.median_filter(image3,**args),('3D median',mode,origin))
            check(ndi.rank_filter(image3,2,**args),('3D rank',mode,origin))
            check(ndi.percentile_filter(image3,75.0,**args),('3D percentile',mode,origin))
    mins3,maxs3,pmin3,pmax3=ndi.extrema(image3,labels3,ids3)
    check(mins3,'3D minima')
    check(maxs3,'3D maxima')
    present3=[np.any(labels3==label) for label in ids3]
    check(np.array([(np.array(pos)+1 if exists else np.zeros(3))
                    for pos,exists in zip(pmin3,present3)]).T,'3D min positions')
    check(np.array([(np.array(pos)+1 if exists else np.zeros(3))
                    for pos,exists in zip(pmax3,present3)]).T,'3D max positions')
    lower3=[]
    upper3=[]
    for label in ids3:
        coords=np.argwhere(labels3==label)
        lower3.append(coords.min(axis=0)+1 if len(coords) else [0,0,0])
        upper3.append(coords.max(axis=0)+1 if len(coords) else [0,0,0])
    check(np.array(lower3).T,'3D lower boxes')
    check(np.array(upper3).T,'3D upper boxes')
    assert cursor==len(data),(cursor,len(data))
    print(f'PASS: {count} labeled-extrema, editable-morphology and rank-footprint SciPy comparisons')


if __name__=='__main__':
    main()
