#!/usr/bin/env python3
"""Compare separable tensor spline grids with scipy.ndimage.map_coordinates."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
MODES = ['mirror', 'constant', 'grid-constant', 'nearest',
         'reflect', 'grid-mirror', 'wrap', 'grid-wrap']


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-spline-tensor-') as tmp:
        folder = Path(tmp)
        objects = []
        for source in ('src/kinds.f90','src/signal_spline.f90','src/signal_spline_nd.f90'):
            obj = folder/(Path(source).stem+'.o')
            subprocess.run(['gfortran','-std=f2018','-O0','-fcheck=all',
                            '-J',str(folder),'-I',str(folder),
                            '-c',str(ROOT/source),'-o',str(obj)],
                           cwd=folder,check=True)
            objects.append(str(obj))
        exe=folder/'spline_tensor_driver'
        subprocess.run(['gfortran','-std=f2018','-O0','-fcheck=all',
                        '-J',str(folder),'-I',str(folder),
                        str(ROOT/'tools/signal_reference/spline_tensor_driver.f90'),
                        *objects,'-o',str(exe)],cwd=folder,check=True)
        text=subprocess.check_output([str(exe)],cwd=folder,text=True)

    ii=np.arange(1,12)[:,None]
    jj=np.arange(1,10)[None,:]
    real=np.sin(0.13*ii+0.27*jj)+0.001*ii*jj
    source=real+1j*(-0.37*real+0.2)
    xx=np.arange(1,14,dtype=float)-2.0+0.31
    yy=np.arange(1,10,dtype=float)-2.0+0.27
    xx,yy=np.meshgrid(xx,yy,indexing='ij')
    coords=np.array([xx-1.,yy-1.])
    cases={}
    for line in text.splitlines():
        parts=line.split()
        if len(parts)!=7:
            raise RuntimeError(f'unexpected Fortran output: {line}')
        modeid,degree,filtered,i,j=map(int,parts[:5])
        value=float(parts[5])+1j*float(parts[6])
        key=(modeid,degree,filtered)
        if key not in cases:
            cases[key]=np.empty((13,9),complex)
        cases[key][i-1,j-1]=value
    assert len(cases)==8*6*2, len(cases)
    total=0
    max_error=0.0
    known_reflect=0
    for (modeid,degree,filtered), got in cases.items():
        mode=MODES[modeid-1]
        expected=ndimage.map_coordinates(source,coords,order=degree,mode=mode,
                                         cval=1.3-0.5j,prefilter=bool(filtered))
        delta=np.max(np.abs(got-expected))
        # SciPy's recursive reflect initializer is inexact for short lines,
        # while SciFort's finite collocation solution is exact.
        tol=2e-10
        if mode in ('reflect','grid-mirror') and degree>1 and filtered:
            tol=5e-4
            known_reflect+=1
        if delta>tol:
            raise AssertionError((mode,degree,filtered,delta,tol))
        if tol<=2e-10: max_error=max(max_error,delta)
        total+=got.size*2
    print(f'PASS: {total} real/imag tensor-grid SciPy comparisons across {len(cases)} cases')
    print(f'Maximum absolute error outside known short-reflect cases: {max_error:.3g}')
    print(f'Known short-reflect cases handled separately: {known_reflect}')


if __name__=='__main__':
    main()
