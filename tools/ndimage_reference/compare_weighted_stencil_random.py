#!/usr/bin/env python3
"""Randomized weighted/directed discrete distance checks using heapq Dijkstra."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from compare_callbacks_weighted_scipy import weighted_stencil_reference

ROOT=Path(__file__).resolve().parents[2]
SOURCES=[ROOT/p for p in (
    'src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90','src/ndimage_distance.f90','src/ndimage_chamfer.f90',
    'src/ndimage_measurements.f90','src/ndimage_watershed.f90','src/ndimage.f90',
    'tools/ndimage_reference/weighted_stencil_random_driver.f90')]


def main():
    rng=np.random.default_rng(20261023)
    cases=[]
    for n in range(96):
        pair=[]
        for shape in ((4,5),(3,3,2)):
            if n==0:
                mask=np.zeros(shape,bool)
            elif n==1:
                mask=np.ones(shape,bool)
            elif n==2:
                mask=np.ones(shape,bool)
                mask.flat[0]=False
            elif n==3:
                mask=np.ones(shape,bool)
                mask.flat[-1]=False
            else:
                mask=rng.random(shape)<(.15+.17*(n%5))
            weights=rng.integers(0,7,size=(3,)*len(shape),dtype=np.int32)
            weights[(1,)*len(shape)]=0
            if n%13==0:
                weights[:]=0
            elif n%13==1:
                # Positive axes only: directional and disconnected graphs.
                weights[:]=0
                for axis in range(len(shape)):
                    edge=[1]*len(shape)
                    edge[axis]=2
                    weights[tuple(edge)]=2+axis
            pair.append((mask,weights))
        cases.append(pair)
    with tempfile.TemporaryDirectory(prefix='scifort-weighted-') as tmp:
        tmp=Path(tmp)
        with (tmp/'weighted_stencil_random_inputs.bin').open('wb') as f:
            for pair in cases:
                for mask,weights in pair:
                    f.write(np.asarray(mask,dtype='<i4').ravel(order='F').tobytes())
                    f.write(np.asarray(weights,dtype='<i4').ravel(order='F').tobytes())
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(tmp),'-I',str(tmp),
                        *map(str,SOURCES),'-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        data=np.fromfile(tmp/'weighted_stencil_random_results.bin',dtype='<i4')
    cursor=0
    count=0
    for case_no,pair in enumerate(cases):
        for mask,weights in pair:
            expected,sites=weighted_stencil_reference(mask,weights)
            n=mask.size
            nd=mask.ndim
            got=data[cursor:cursor+n].reshape(mask.shape,order='F')
            cursor+=n
            got_sites=data[cursor:cursor+n*nd].reshape((nd,n),order='F')
            cursor+=n*nd
            if not np.array_equal(got,expected) or not np.array_equal(got_sites,sites):
                raise AssertionError((case_no,mask.shape,'distance difference',
                                      got,expected,got_sites,sites))
            count+=n+nd*n
    assert cursor==data.size,(cursor,data.size)
    print(f'PASS: {count} randomized directed-weight and feature comparisons')


if __name__=='__main__':
    main()
