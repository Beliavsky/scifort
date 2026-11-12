#!/usr/bin/env python3
"""Numerically compare cached SciFort spline plans with SciPy map_coordinates.

Coordinates and output are emitted by a separate Fortran driver. SciFort uses
one-based coordinates; SciPy uses zero-based coordinates.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
CASES = [
    (0, 'mirror'), (1, 'mirror'), (3, 'mirror'), (5, 'mirror'),
    (0, 'reflect'), (1, 'reflect'),
    (3, 'nearest'), (5, 'nearest'),
    (3, 'grid-wrap'), (5, 'grid-wrap'),
    (3, 'grid-constant'), (5, 'grid-constant'),
    (1, 'constant'), (3, 'constant'), (5, 'constant'),
]
COORDS = np.empty((2,96),dtype=np.float64)
COORDS[:,:6] = np.array([[1.0,3.25,6.75,8.8,-0.8,3.9],
                          [1.0,2.5,5.2,7.3,5.0,-1.2]], dtype=np.float64)
for point in range(7,97):
    COORDS[0,point-1]=-2.0+(point*17 % 117)/10.0
    COORDS[1,point-1]=-1.75+(point*31 % 108)/10.0

X = np.array([[np.sin(0.23*i+0.57*j)+0.1*i*j for j in range(1,8)]
              for i in range(1,9)],dtype=np.float64)


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-plan-') as temp:
        tmp=Path(temp)
        command=['gfortran','-std=f2018','-O0','-fcheck=all','-J',str(tmp),'-I',str(tmp)]
        objects=[]
        for src in ('kinds.f90','signal_spline.f90','signal_spline_nd.f90'):
            obj=tmp/(src.removesuffix('.f90')+'.o')
            subprocess.run([*command,'-c',str(ROOT/'src'/src),'-o',str(obj)],check=True,cwd=tmp)
            objects.append(str(obj))
        exe=tmp/'spline_plan_reference_driver'
        subprocess.run([*command,*objects,str(ROOT/'tools/signal_reference/spline_plan_reference_driver.f90'),
                        '-o',str(exe)],check=True,cwd=tmp)
        output=subprocess.check_output([str(exe)],text=True,cwd=tmp)
    lines=output.strip().splitlines()
    assert len(lines)==len(CASES)*len(COORDS.T), (len(lines),len(CASES)*len(COORDS.T))
    biggest=0.0
    for case_idx,(order,mode) in enumerate(CASES):
        actual=[]
        for point_idx in range(COORDS.shape[1]):
            s=lines[case_idx*COORDS.shape[1]+point_idx].split()
            assert int(s[0])==case_idx+1 and int(s[1])==point_idx+1
            actual.append(float(s[2]))
        expected=ndimage.map_coordinates(X, COORDS-1.0, order=order,mode=mode,
                                          cval=-2.5,prefilter=True)
        error=np.max(np.abs(expected-np.array(actual)))
        if error>3e-8:
            raise AssertionError(f'mode={mode} order={order} error={error}: actual={actual} expected={expected}')
        biggest=max(biggest,error)
    print(f'PASS: {len(lines)} spline plan vs SciPy values; maximum error {biggest:.3g}')


if __name__=='__main__':
    main()
