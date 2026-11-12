#!/usr/bin/env python3
"""SciPy comparisons for 1D kernels, generic/labeled callbacks and weighted L1."""
from pathlib import Path
import subprocess
import tempfile
import heapq
import numpy as np
from scipy import ndimage as ndi

ROOT=Path(__file__).resolve().parents[2]
SOURCES=[ROOT/p for p in (
    'src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90','src/ndimage_distance.f90','src/ndimage_chamfer.f90',
    'src/ndimage_measurements.f90','src/ndimage_watershed.f90','src/ndimage.f90',
    'tools/ndimage_reference/callbacks_weighted_driver.f90')]
MODES=['reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']


def weighted_sum(v):
    return np.dot(np.arange(1,len(v)+1,dtype=np.float64),v)


def weighted_stencil_reference(binary,weights):
    """Multi-source Dijkstra, including directed steps and Fortran-order ties."""
    dims=binary.shape
    total=binary.size
    dist=np.full(total,2**55,dtype=np.int64)
    sites=np.full(total,-1,dtype=np.int64)
    queue=[]
    for loc in np.ndindex(dims):
        pos=np.ravel_multi_index(loc,dims,order='F')
        if not binary[loc]:
            dist[pos]=0
            sites[pos]=pos
            heapq.heappush(queue,(0,pos,pos))
    offsets=[]
    for off in np.ndindex(weights.shape):
        w=weights[off]
        if w:
            offsets.append((np.array(off,dtype=int)-1,int(w)))
    while queue:
        cost,site,pos=heapq.heappop(queue)
        if cost!=dist[pos] or site!=sites[pos]:
            continue
        location=np.array(np.unravel_index(pos,dims,order='F'))
        for off,w in offsets:
            dst=location-off
            if np.any(dst<0) or np.any(dst>=dims):
                continue
            other=np.ravel_multi_index(tuple(dst),dims,order='F')
            new=cost+w
            if (new,site)<(dist[other],sites[other]):
                dist[other]=new
                sites[other]=site
                heapq.heappush(queue,(new,site,other))
    found=dist<2**55
    result=np.where(found,dist,-1).astype(np.int32).reshape(dims,order='F')
    features=np.zeros((binary.ndim,total),dtype=np.int32)
    for pos in range(total):
        if found[pos]:
            features[:,pos]=np.array(np.unravel_index(sites[pos],dims,order='F'))+1
    return result,features


def main():
    i,j=np.mgrid[1:6,1:7]
    a=np.sin(.3*i+.7*j)+i-.4*j
    z=a+1j*np.cos(.7*i-.2*j)
    labels=(i+3*j)%7-3
    mask=(i*i+2*j)%7 != 0
    i,j,k=np.mgrid[1:5,1:4,1:3]
    a3=np.cos(.6*i-j+2*k)+.7*i*j-k
    labels3=(i+2*j+3*k)%6-2
    mask3=(i+2*j+4*k)%5 != 0
    weights=np.array([1.3,-.4,.2,.7])
    cweights=np.array([.4+.7j,1.1-.3j,-.8+.9j])
    foot=np.array([1,0,1,1,0,1],bool).reshape((2,3),order='F')
    foot3=np.array([1,0,1,1,1,0,0,1,1,0,1,0],bool).reshape((3,2,2),order='F')
    ids=np.array([-2,0,3,11,-2,2])
    with tempfile.TemporaryDirectory(prefix='scifort-callbacks-') as tmp:
        tmp=Path(tmp)
        subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=all',
                        '-Werror=implicit-interface', '-J',str(tmp),'-I',str(tmp),
                        *map(str,SOURCES),'-o',str(tmp/'driver')],cwd=tmp,check=True)
        subprocess.run([str(tmp/'driver')],cwd=tmp,check=True)
        data=(tmp/'callbacks_weighted_results.bin').read_bytes()
    offset=0
    checked=0
    def check(value,name,dtype=np.float64,exact=False):
        nonlocal offset,checked
        value=np.asarray(value,dtype=dtype)
        expected=value.ravel(order='F')
        nbytes=expected.nbytes
        got=np.frombuffer(data,dtype=dtype,count=expected.size,offset=offset)
        if (exact and not np.array_equal(got,expected)) or (not exact and not np.allclose(
            got,expected,atol=5e-12,rtol=1e-12,equal_nan=True)):
            raise AssertionError((name,np.max(np.abs(got-expected)),got[:20],expected[:20]))
        offset+=nbytes
        checked+=expected.size
    for mode in MODES:
        for ax in (1,2):
            for w in (2,3,4):
                for origin in range(-(w//2),(w-1)//2+1):
                    for f in (ndi.correlate1d,ndi.convolve1d):
                        check(f(a,weights[:w],axis=ax-1,origin=origin,mode=mode,cval=-1.25),
                              (f.__name__,'real',mode,ax,w,origin))
        for ax in (1,2,3):
            for w in (2,3):
                for origin in range(-(w//2),(w-1)//2+1):
                    for f in (ndi.correlate1d,ndi.convolve1d):
                        check(f(a3,weights[:w],axis=ax-1,origin=origin,mode=mode,cval=1.75),
                              (f.__name__,'3d',mode,ax,w,origin))
        for ax in (1,2):
            for w in (2,3):
                for origin in range(-(w//2),(w-1)//2+1):
                    for f in (ndi.correlate1d,ndi.convolve1d):
                        check(f(z,cweights[:w],axis=ax-1,origin=origin,mode=mode,cval=-1.25+.3j),
                              (f.__name__,'complex',mode,ax,w,origin),dtype=np.complex128)
    for mode in MODES:
        for kw in ({'size':(2,3),'origin':(-1,1)},
                   {'footprint':foot,'origin':(0,-1)},
                   {'size':(3,3),'origin':(1,-1)}):
            check(ndi.generic_filter(a,weighted_sum,mode=mode,cval=-2.5,**kw),
                  ('generic',mode,kw['origin']))
        check(ndi.generic_filter(a3,weighted_sum,footprint=foot3,
                                 origin=(-1,0,-1),mode=mode,cval=1.125),
              ('generic3d',mode))
    for values,labeled in ((a,labels),(a3,labels3)):
        for fun in (np.max,np.var,weighted_sum):
            if fun is weighted_sum:
                # SciPy argsort() groups pixels by label but is not stable; our
                # callback preserves input C order instead of sort artifacts.
                expected=[weighted_sum(values[labeled==v]) if np.any(labeled==v) else -11.0 for v in ids]
            else:
                expected=ndi.labeled_comprehension(values,labeled,ids,fun,np.float64,-11.0)
            check(expected,('comprehension',values.shape,fun.__name__))
    for weights2,weights3 in (([1,1],[1,1,1]),([2,3],[1,3,5]),([3,2],[4,2,1]),
                               ([1,6],[6,1,2]),([4,5],[2,4,3])):
        for binary,cost in ((mask,weights2),(mask3,weights3)):
            sites=np.argwhere(~binary)
            expected=np.zeros(binary.shape,dtype=np.int32)
            for pos in np.ndindex(binary.shape):
                expected[pos]=np.min(np.sum(np.abs(sites-np.array(pos))*cost,axis=1))
            check(expected,('weighted-distance',binary.shape,cost),dtype=np.int32,exact=True)
            # Fortran indices are stored as (rank, Fortran-flattened size), one-based.
            nd=binary.ndim
            indices=np.frombuffer(data,dtype=np.int32,count=nd*binary.size,offset=offset)
            site_coords=indices.reshape((nd,binary.size),order='F')-1
            for pos in np.ndindex(binary.shape):
                flat=np.ravel_multi_index(pos,binary.shape,order='F')
                site=site_coords[:,flat]
                if np.any(site<0) or np.any(site>=binary.shape) or binary[tuple(site)]:
                    raise AssertionError(('weighted feature not background',pos,site))
                dist=np.sum(np.abs(site-np.array(pos))*cost)
                if dist!=expected[pos]:
                    raise AssertionError(('weighted feature not closest',pos,site,dist,expected[pos]))
            offset+=indices.nbytes
            checked+=indices.size
    stencils=[]
    for case in range(5):
        stencil=np.zeros((3,3),np.int32)
        if case==0:
            stencil[0,1]=stencil[2,1]=stencil[1,0]=stencil[1,2]=1
        elif case==1:
            stencil[:]=1
            stencil[1,1]=0
        elif case==2:
            stencil[:]=4
            stencil[1,1]=0
            stencil[0,1]=stencil[2,1]=stencil[1,0]=stencil[1,2]=3
        elif case==3:
            stencil[0,1],stencil[2,1],stencil[1,0],stencil[1,2]=3,2,4,1
        else:
            stencil[0,0],stencil[0,2],stencil[2,0],stencil[2,2]=1,2,3,4
        stencils.append(stencil)
    stencil3=np.fromfunction(lambda i,j,k: abs(i-1)+abs(j-1)+abs(k-1),(3,3,3),dtype=int)
    for binary,stencil in [(mask,item) for item in stencils]+[(mask3,stencil3)]:
        distances,sites=weighted_stencil_reference(binary,stencil)
        check(distances,('stencil chamfer',binary.shape),dtype=np.int32,exact=True)
        check(sites,('stencil chamfer features',binary.shape),dtype=np.int32,exact=True)
        if np.array_equal(stencil,stencils[0]):
            check_ndi=ndi.distance_transform_cdt(binary,metric='taxicab')
            assert np.array_equal(distances,check_ndi)
        elif np.array_equal(stencil,stencils[1]):
            check_ndi=ndi.distance_transform_cdt(binary,metric='chessboard')
            assert np.array_equal(distances,check_ndi)
    assert offset==len(data),(offset,len(data))
    print(f'PASS: {checked} SciPy/independent weighted-chamfer scalar comparisons')

if __name__=='__main__':
    main()
