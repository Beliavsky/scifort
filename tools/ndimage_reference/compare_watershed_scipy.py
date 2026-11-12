#!/usr/bin/env python3
"""Compare bucket-queue watershed with SciPy, including negative marker priority."""
from pathlib import Path
import subprocess
import tempfile
from collections import deque
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]


def bounded_ift(image, markers, structure, legacy=False):
    """Slow independent priority-flood oracle enforcing every axis bound."""
    shape=image.shape
    image=np.asarray(image).ravel(order='C')
    values=np.array(markers,dtype=np.int32).ravel(order='C')
    offsets=np.array(np.nonzero(structure)).T-1
    offsets=offsets[np.any(offsets!=0,axis=1)]
    stride=np.array([np.prod(shape[k+1:],dtype=int) for k in range(len(shape))])
    cost=np.full(image.size, int(image.max())+1)
    done=np.zeros(image.size,bool)
    buckets=[deque() for _ in range(int(image.max())+1)]
    for p in np.flatnonzero(values):
        cost[p]=0
        if values[p]>0:buckets[0].appendleft(int(p))
        else:buckets[0].append(int(p))
    for bucket in buckets:
        while bucket:
            p=bucket.popleft()
            done[p]=True
            coord=np.unravel_index(p,shape)
            for delta in offsets:
                newcoord=np.array(coord)+delta
                if not legacy and (np.any(newcoord<0) or np.any(newcoord>=shape)):continue
                q=int(np.dot(stride,newcoord))
                if q<0 or q>=image.size:continue
                if done[q]:continue
                nextcost=max(int(cost[p]), abs(int(image[p])-int(image[q])))
                if nextcost>=cost[q]:continue
                if cost[q] <= int(image.max()):
                    buckets[cost[q]].remove(q)
                cost[q]=nextcost
                values[q]=values[p]
                if values[q]>0:buckets[nextcost].appendleft(q)
                else:buckets[nextcost].append(q)
    return values.reshape(shape)


def main():
    sources=[ROOT/z for z in ('src/kinds.f90','src/signal_spline.f90',
          'src/signal_spline_nd.f90','src/ndimage_morphology.f90',
          'src/ndimage_distance.f90','src/ndimage_chamfer.f90','src/ndimage_measurements.f90','src/ndimage_watershed.f90',
               'src/ndimage.f90',
          'tools/ndimage_reference/watershed_driver.f90')]
    with tempfile.TemporaryDirectory(prefix='scifort-watershed-') as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                 '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                 *map(str,sources),'-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        raw=np.fromfile(tmp/'watershed_results.bin',np.int32)
        count=0
        scipy_matched=0
        offset=0
        i,j=np.mgrid[1:7,1:8]
        for t in range(1,25):
            x=(11*i*i+17*i*j+13*j*j+9*t*i+2*t*j)%(13+t*3)
            markers=np.zeros((6,7),dtype=np.int32)
            markers[1,1]=1
            markers[4,5]=2
            if t%3==0: markers[0,6]=-3
            if t%4==0: markers[3,2]=4
            st=np.zeros((3,3),dtype=bool)
            switch=t%4
            if switch==0: st[:]=True
            elif switch==1: st[1,:]=True;st[:,1]=True
            elif switch==2: st[0,0]=st[1,1]=st[2,2]=True
            else: st[0,1]=st[1,1]=st[2,1]=True
            expected=ndi.watershed_ift(x.astype(np.uint16),markers,st)
            actual=raw[offset:offset+42].reshape((6,7),order='F')
            bounded=bounded_ift(x,markers,st)
            assert np.array_equal(actual,bounded),(t,actual,bounded)
            offset+=42
            compatible=raw[offset:offset+42].reshape((6,7),order='F')
            assert np.array_equal(compatible,bounded_ift(x,markers,st,legacy=True)),t
            scipy_matched+=int(np.array_equal(compatible,expected))
            offset+=42;count+=1
        i,j,k=np.mgrid[1:5,1:4,1:4]
        for t in range(1,13):
            x=(5*i*i+9*i*j+7*k+11*j*k+t*i)%(19+t)
            markers=np.zeros((4,3,3),dtype=np.int32)
            markers[0,0,0]=1
            markers[3,2,2]=-2
            if t%2==0:markers[1,1,1]=3
            expected=ndi.watershed_ift(x.astype(np.uint16),markers)
            actual=raw[offset:offset+36].reshape((4,3,3),order='F')
            bounded=bounded_ift(x,markers,ndi.generate_binary_structure(3,1))
            assert np.array_equal(actual,bounded),(t,actual,bounded)
            offset+=36
            compatible=raw[offset:offset+36].reshape((4,3,3),order='F')
            assert np.array_equal(compatible,bounded_ift(x,markers,ndi.generate_binary_structure(3,1),legacy=True)),t
            scipy_matched+=int(np.array_equal(compatible,expected))
            offset+=36;count+=1
        assert offset==raw.size
        print(f'PASS: watershed {count} bounded and {count} legacy reference arrays; {scipy_matched}/{count} SciPy exact')


if __name__=='__main__':main()
