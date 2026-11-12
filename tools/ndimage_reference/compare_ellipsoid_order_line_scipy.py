#!/usr/bin/env python3
"""Compare Fourier ellipsoid, mixed-mode rank/median and line callbacks with SciPy."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage as ndi, special
from compare_mixed_complex_scipy import map_boundary

ROOT=Path(__file__).resolve().parents[2]
SOURCES=[ROOT/p for p in (
    'src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90',
    'src/ndimage_morphology.f90','src/ndimage_distance.f90','src/ndimage_chamfer.f90',
    'src/ndimage_measurements.f90','src/ndimage_watershed.f90','src/ndimage.f90',
    'tools/ndimage_reference/ellipsoid_order_line_driver.f90')]
MODES=['reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']


def coordinate_reference(x,window,origins,modes,cval,kind,percent=50.,rank=0,footprint=None):
    """Per-axis boundary extension followed by a finite window order statistic."""
    if footprint is None: footprint=np.ones(window,dtype=bool)
    offsets=np.argwhere(footprint)-np.array(window)//2-np.array(origins)
    y=np.empty_like(x)
    for coordinate in np.ndindex(x.shape):
        values=[]
        for off in offsets:
            coords=np.array(coordinate)+off
            valid=True
            for d in range(x.ndim):
                idx=map_boundary(int(coords[d]),x.shape[d],modes[d])
                if idx is None:
                    valid=False
                    break
                coords[d]=idx
            values.append(x[tuple(coords)] if valid else cval)
        values.sort()
        if kind=='rank': chosen=rank%len(values)
        elif kind=='percentile': chosen=min(len(values)-1,int(len(values)*percent/100.))
        else: chosen=len(values)//2
        y[coordinate]=values[chosen]
    return y


def bessel_filter(a,sizes,n=-1,axis=-1):
    """Independent Fourier-ball transform for SciFort's rank > 3 extension."""
    size=np.broadcast_to(np.asarray(sizes,dtype=float),a.ndim)
    r2=np.zeros(a.shape)
    for d,(count,diameter) in enumerate(zip(a.shape,size)):
        length=n if n>0 and d==axis%a.ndim else count
        if n>0 and d==axis%a.ndim:
            frequencies=np.arange(count,dtype=float)/length
        else:
            frequencies=np.fft.fftfreq(count)
        reshape=[1]*a.ndim
        reshape[d]=count
        r2+=(diameter*frequencies.reshape(reshape))**2
    q=np.pi*np.sqrt(r2)
    m=a.ndim/2
    factor=np.ones_like(q)
    positive=q>0
    factor[positive]=special.gamma(m+1)*(2/q[positive])**m*special.jv(m,q[positive])
    return a*factor


def line_reduce(input_line,output_line):
    count=output_line.size
    width=input_line.size-count+1
    output_line[:]=input_line[:count]+2*input_line[1:count+1]-.5*input_line[width-1:width-1+count]


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-ellipsoid-order-line-') as temp:
        temp=Path(temp)
        subprocess.run(['gfortran','-std=f2018','-O0','-fcheck=all',
                        '-Werror=implicit-interface','-J',str(temp),'-I',str(temp),
                        *map(str,SOURCES),'-o',str(temp/'driver')],cwd=temp,check=True)
        subprocess.run([str(temp/'driver')],cwd=temp,check=True)
        data=(temp/'ellipsoid_order_line_results.bin').read_bytes()
    offset=0
    total=0
    def check(expected,name):
        nonlocal offset,total
        expected=np.asarray(expected)
        dtype=np.complex128 if np.iscomplexobj(expected) else np.float64
        values=expected.astype(dtype).ravel(order='F')
        got=np.frombuffer(data,dtype=dtype,count=values.size,offset=offset)
        if not np.allclose(got,values,atol=2e-12,rtol=3e-12):
            dif=np.abs(got-values)
            index=np.argmax(dif)
            raise AssertionError((name,'error',dif[index],'index',index,
                                  'actual',got[index],'expected',values[index]))
        offset+=values.nbytes
        total+=values.size
    a1_i=np.arange(1,8)
    f1=np.sin(.31*a1_i)+1j*np.cos(.47*a1_i)
    for s in (2.3,-3.1): check(ndi.fourier_ellipsoid(f1,s),('ellipsoid 1d',s))
    i,j=np.mgrid[1:6,1:7]
    f2=np.sin(.4*i+.23*j)+1j*np.cos(.27*i-.49*j)
    for s in ((1.2,3.1),2.7): check(ndi.fourier_ellipsoid(f2,s),('ellipsoid 2d',s))
    i,j,k=np.mgrid[1:5,1:4,1:3]
    f3=np.sin(.4*i+.7*j-.2*k)+1j*np.sin(.31*i*j*k)
    for s in ((1.8,.2,3.4),1.3): check(ndi.fourier_ellipsoid(f3,s),('ellipsoid 3d',s))
    i,j=np.mgrid[1:6,1:5]
    r2=np.sin(.3*i+j)+1j*np.cos(.44*i-j)
    check(ndi.fourier_ellipsoid(r2,(.9,2.4),n=8,axis=0),'ellipsoid packed first')
    i,j,k=np.mgrid[1:3,1:4,1:6]
    r3=np.sin(.5*i+.6*j+.3*k)+1j*np.cos(.41*i-.53*j*k)
    check(ndi.fourier_ellipsoid(r3,(1.1,1.6,3.8),n=9,axis=-1),'ellipsoid packed last')
    i,j,k,l=np.mgrid[1:4,1:3,1:3,1:3]
    a4=np.sin(.17*i*j+.11*k*l)+1j*np.cos(.21*i*k-.13*j*l)
    check(bessel_filter(a4,(1,2,3,4)),'ellipsoid rank 4 analytical')
    a15=np.full((2,)+(1,)*14,2-3j,dtype=np.complex128)
    check(bessel_filter(a15,1.5),'ellipsoid rank 15 analytical')
    for diameter in (5.2,12.,43.,100.):
        check(bessel_filter(a15,diameter),('ellipsoid rank15 large q',diameter))

    i,j=np.mgrid[1:6,1:7]
    a=np.sin(.4*i+.7*j)+(i*i-j)/7.
    footprint=np.array([1,0,1,1,0,1],bool).reshape((2,3),order='F')
    for mode1 in MODES:
        for mode2 in MODES:
            modes=(mode1,mode2)
            for kind,foot,rank,percent in [('median',None,0,50),('median',footprint,0,50),
                                           ('rank',footprint,-2,0),('percentile',footprint,0,35)]:
                expected=coordinate_reference(a,(2,3),(-1,1),modes,-.7,kind,percent,rank,foot)
                if mode1==mode2:
                    args={'origin':(-1,1),'mode':mode1,'cval':-.7}
                    if foot is not None: args['footprint']=foot
                    else: args['size']=(2,3)
                    if kind=='median': scipy=ndi.median_filter(a,**args)
                    elif kind=='rank': scipy=ndi.rank_filter(a,-2,**args)
                    else: scipy=ndi.percentile_filter(a,35,**args)
                    np.testing.assert_allclose(expected,scipy,atol=1e-14,
                                               err_msg=str((kind,modes)))
                check(expected,('mixed order',kind,modes))

    footprint3=np.array([1,1,0,1,1,0,0,1,1,0,1,1],bool).reshape((2,3,2),order='F')
    i,j,k=np.mgrid[1:4,1:5,1:3]
    c=np.sin(.23*i+.51*j+.37*k)+i*k*.09
    for p in range(8):
        modes=(MODES[p],MODES[(p+1)%8],MODES[(p+4)%8])
        for kind,rank,percent in [('median',0,50),('rank',2,0),('percentile',0,72)]:
            check(coordinate_reference(c,(2,3,2),(-1,0,-1),modes,1.3,
                                       kind,percent,rank,footprint3),('mixed order 3D',kind,modes))
    for mode in MODES:
        for arr,axis,width,origin,cval in ((a,0,4,1,-.7),(a,1,3,-1,-.7),(c,2,2,0,1.3)):
            expected=ndi.generic_filter1d(arr,line_reduce,width,axis=axis,mode=mode,
                                          cval=cval,origin=origin)
            check(expected,('line callback',mode,axis,width,origin))
    assert offset==len(data),(offset,len(data))
    print(f'PASS: {total} Fourier ellipsoid, mixed rank, and SciPy generic_filter1d reference values')

if __name__=='__main__': main()
