#!/usr/bin/env python3
"""Verify analytic Fortran Hessians against independently differentiated SciPy splines."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[2]
MODES=('mirror','constant','grid-constant','nearest',
       'reflect','grid-mirror','wrap','grid-wrap')
COORDS=np.array([[2.35,7.28,14.17,-1.23,18.22,5.33,8.41,1.23],
                 [3.42,5.67,12.39,4.51,10.73,-0.74,16.32,1.47]])


def run():
    with tempfile.TemporaryDirectory(prefix='scifort-hessian-') as path:
        wd=Path(path)
        cmd=['gfortran','-std=f2018','-O0','-fcheck=all',
             '-ffpe-trap=invalid,zero,overflow','-Werror=implicit-interface',
             '-J',str(wd),'-I',str(wd)]
        objects=[]
        for name in ('kinds.f90','signal_spline.f90','signal_spline_nd.f90'):
            obj=wd/(Path(name).stem+'.o')
            subprocess.run([*cmd,'-c',str(ROOT/'src'/name),'-o',str(obj)],
                           cwd=wd,check=True)
            objects.append(str(obj))
        exe=wd/'hessian_driver'
        subprocess.run([*cmd,*objects,str(ROOT/'tools/signal_reference/spline_hessian_driver.f90'),
                        '-o',str(exe)],cwd=wd,check=True)
        output=subprocess.check_output([str(exe)],cwd=wd,text=True)
    result={}
    for line in output.splitlines():
        fields=line.split()
        if len(fields)!=16:
            raise AssertionError(f'malformed Fortran output: {line}')
        mi,degree,kind,k=map(int,fields[:4])
        nums=np.asarray([float(t) for t in fields[4:]]).reshape(6,2)
        result[(mi,degree,kind,k)]=nums[:,0]+1j*nums[:,1]
    assert len(result)==8*6*2*8
    i=np.arange(1,18)[:,None]
    j=np.arange(1,17)[None,:]
    x=np.sin(0.19*i+0.23*j)+0.003*i*j
    step=5.0e-4
    checked=0
    regular_max=0.0
    reflect_max=0.0
    for mi,mode in enumerate(MODES,1):
        for degree in range(6):
            for kind in range(2):
                source=x if kind==0 else x+1j*(0.37*x-0.12)
                fill=0.43 if kind==0 else 0.43-0.29j
                def f(delta):
                    return ndimage.map_coordinates(source,COORDS-1+delta,
                        order=degree,mode=mode,cval=fill,prefilter=True)
                delta=np.zeros_like(COORDS)
                f0=f(delta)
                plus=[]
                minus=[]
                for axis in range(2):
                    move=np.zeros_like(delta); move[axis,:]=step
                    plus.append(f(move))
                    minus.append(f(-move))
                move=np.full_like(delta,step)
                mixed=(f(move)-f(move*np.array([[1],[-1]]))-
                       f(move*np.array([[-1],[1]]))+f(-move))/(4*step**2)
                expected=np.array([f0,
                    (plus[0]-minus[0])/(2*step),
                    (plus[1]-minus[1])/(2*step),
                    (plus[0]-2*f0+minus[0])/(step**2),
                    mixed,
                    (plus[1]-2*f0+minus[1])/(step**2)])
                got=np.column_stack([result[(mi,degree,kind,k)] for k in range(1,9)])
                err=float(np.max(np.abs(got-expected)))
                # SciPy's finite line reflect prefilter and the exact
                # collocation solver need not agree on high-order short lines.
                tol=1.0e-3 if degree>1 and mode in ('reflect','grid-mirror') else 1.2e-3
                if err>tol:
                    pos=np.unravel_index(np.argmax(abs(got-expected)),got.shape)
                    raise AssertionError((mode,degree,kind,pos,err,got[pos],expected[pos]))
                if mode in ('reflect','grid-mirror'):
                    reflect_max=max(reflect_max,err)
                else:
                    regular_max=max(regular_max,err)
                checked+=got.size*(2 if kind else 1)
    print(f'PASS: {checked} real/imag SciPy value-gradient-Hessian comparisons')
    print(f'Maximum numerical differences: regular={regular_max:.3g}, reflect={reflect_max:.3g}')


if __name__=='__main__':
    run()
