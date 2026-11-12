#!/usr/bin/env python3
"""Compare streamed real affine and complex geometric plans to SciPy."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[2]
MODES=['mirror','constant','grid-constant','nearest',
       'reflect','grid-mirror','wrap','grid-wrap']
MATRIX=np.array([[1.109,-0.123],[0.187,0.943]])
OFFSET=np.array([-2.4131,1.2373])


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-plan-transform-') as name:
        tmp=Path(name)
        flags=['gfortran','-std=f2018','-O0','-fcheck=all',
               '-ffpe-trap=invalid,zero,overflow','-J',str(tmp),'-I',str(tmp)]
        objects=[]
        for source in ('kinds.f90','signal_spline.f90','signal_spline_nd.f90'):
            obj=tmp/(Path(source).stem+'.o')
            subprocess.run([*flags,'-c',str(ROOT/'src'/source),'-o',str(obj)],
                           cwd=tmp,check=True)
            objects.append(str(obj))
        exe=tmp/'spline_plan_transform_driver'
        subprocess.run([*flags,*objects,
                        str(ROOT/'tools/signal_reference/spline_plan_transform_driver.f90'),
                        '-o',str(exe)],cwd=tmp,check=True)
        output=subprocess.check_output([str(exe)],cwd=tmp,text=True)

    ii=np.arange(1,20)[:,None]
    jj=np.arange(1,18)[None,:]
    x=np.sin(0.17*ii+0.31*jj)+0.004*ii*jj
    z=x+1j*(-0.37*x+0.29)
    actual={}
    for line in output.splitlines():
        fields=line.split()
        if len(fields)!=7:
            raise ValueError('Unexpected Fortran output: '+line)
        mi,order,kind,i,j=map(int,fields[:5])
        if (mi,order,kind) not in actual:
            actual[mi,order,kind]=np.empty((12,11),dtype=np.complex128)
        actual[mi,order,kind][i-1,j-1]=float(fields[5])+1j*float(fields[6])
    assert len(actual)==8*6*2, len(actual)
    count=0
    largest=0.0
    special=[]
    for (mi,order,kind),got in actual.items():
        mode=MODES[mi-1]
        data=x if kind==0 else z
        fill=1.125 if kind==0 else 1.125-0.35j
        expected=ndimage.affine_transform(data,MATRIX,OFFSET,
                  output_shape=(12,11),order=order,mode=mode,cval=fill,prefilter=True)
        error=np.max(np.abs(got-expected))
        tolerance=3.0e-8
        if mode in ('reflect','grid-mirror') and order>1:
            # Reflect initialization in SciPy uses a shortened recursion;
            # SciFort solves the finite spline collocation problem exactly.
            tolerance=1.0e-5
            special.append((mode,order,kind,error))
        if error>tolerance:
            raise AssertionError(f'{mode} order={order} complex={kind}: {error:g} > {tolerance:g}')
        if tolerance<=3e-8:
            largest=max(largest,error)
        count+=got.size*(2 if kind else 1)
    print(f'PASS: {count} streamed spline-transform scalar SciPy comparisons')
    print(f'Maximum error outside reflect exceptions: {largest:.3g}')
    if special:
        print(f'Reflect/grid-mirror high-order cases checked: {len(special)}, '
              f'maximum error {max(row[3] for row in special):.3g}')


if __name__=='__main__':
    main()
