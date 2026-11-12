#!/usr/bin/env python3
"""Compare grid/affine spline gradients with independent SciPy finite differences."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
MODES = ('mirror','constant','grid-constant','nearest','reflect','grid-mirror','wrap','grid-wrap')
DEGREES = (0, 1, 3, 5)
A = np.array([[0.48, 0.19], [-0.15, 0.31]], dtype=np.float64)
SHAPE = (6,5)

def main():
    with tempfile.TemporaryDirectory(prefix='scifort-grid-affine-gradient-') as td:
        wd = Path(td)
        flags = ['gfortran','-std=f2018','-O0','-fcheck=all',
                 '-ffpe-trap=invalid,zero,overflow','-J',str(wd),'-I',str(wd)]
        objects = []
        for name in ('kinds.f90','signal_spline.f90','signal_spline_nd.f90'):
            obj=wd/(Path(name).stem+'.o')
            subprocess.run([*flags,'-c',str(ROOT/'src'/name),'-o',str(obj)],cwd=wd,check=True)
            objects.append(str(obj))
        exe=wd/'driver'
        subprocess.run([*flags,*objects,str(ROOT/'tools/signal_reference/spline_grid_affine_gradient_driver.f90'),
                        '-o',str(exe)],cwd=wd,check=True)
        out=subprocess.check_output([str(exe)],cwd=wd,text=True)

    data={}
    for line in out.splitlines():
        fields=line.split()
        if len(fields)!=20:
            raise ValueError(f'Unexpected driver record: {line!r}')
        mi,degree,kind,icase,i,j=map(int,fields[:6])
        vals=np.asarray([float(t) for t in fields[6:]]).reshape(-1,2)
        data[(mi,degree,kind,icase,i,j)] = vals[:,0] + 1j*vals[:,1]
    assert len(data)==8*len(DEGREES)*2*2*SHAPE[0]*SHAPE[1],len(data)

    ii=np.arange(1,18)[:,None]
    jj=np.arange(1,17)[None,:]
    x=np.sin(0.19*ii+0.23*jj)+0.003*ii*jj
    z=x+1j*(0.31*x-0.17)
    outcoords=np.indices(SHAPE, dtype=np.float64).reshape(2,-1)
    count=0
    maximum=0.0
    exceptional=0.0
    step=1.0e-5
    for mi in range(1,9):
        mode=MODES[mi-1]
        for degree in DEGREES:
            for kind in (0,1):
                source=x if kind==0 else z
                fill=0.43 if kind==0 else (0.43-0.29j)
                for icase in (1,2):
                    offset=np.array([1.13,1.07] if icase==1 else [-2.73,13.68])
                    def evaluate(mat,off):
                        points=off[:,None]+mat@outcoords
                        return ndimage.map_coordinates(source,points,order=degree,
                            mode=mode,cval=fill,prefilter=True)
                    expected=np.empty((7,outcoords.shape[1]),dtype=np.complex128)
                    expected[0]=evaluate(A,offset)
                    for ia in range(2):
                        for ja in range(2):
                            plus=A.copy(); minus=A.copy()
                            plus[ia,ja]+=step; minus[ia,ja]-=step
                            expected[1+ia+2*ja]=(evaluate(plus,offset)-evaluate(minus,offset))/(2*step)
                    for ia in range(2):
                        plus=offset.copy(); minus=offset.copy()
                        plus[ia]+=step; minus[ia]-=step
                        expected[5+ia]=(evaluate(A,plus)-evaluate(A,minus))/(2*step)
                    got=np.empty_like(expected)
                    for k in range(outcoords.shape[1]):
                        i,j=outcoords[:,k].astype(int)+1
                        got[:,k]=data[(mi,degree,kind,icase,i,j)]
                    error=float(np.max(np.abs(got-expected)))
                    tolerance=8e-5
                    if degree>1 and mode in ('reflect','grid-mirror'):
                        tolerance=4e-3  # documented recursive versus exact finite-line prefilter
                    if error>tolerance:
                        raise AssertionError((mode,degree,kind,icase,error,tolerance,
                                              np.unravel_index(np.argmax(abs(got-expected)),got.shape)))
                    if tolerance>1e-3:
                        exceptional=max(exceptional,error)
                    else:
                        maximum=max(maximum,error)
                    count+=got.size*(2 if kind else 1)
    print(f'PASS: {count} real/imag grid and affine-parameter SciPy comparisons')
    print(f'Max regular-case error: {maximum:.3g}; short-line reflect difference: {exceptional:.3g}')

if __name__=='__main__':
    main()
