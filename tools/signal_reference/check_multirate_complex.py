#!/usr/bin/env python3
"""Cross-check complex FIR resampling and padding modes against SciPy 1.17.

Usage: python check_multirate_complex.py path/to/multirate_complex_driver
"""
import re
import subprocess
import sys
import numpy as np
from scipy import signal


def arrays():
    i1 = np.arange(1, 8)
    r1 = np.sin(.31*i1) + .07*i1
    c1 = r1 + 1j*np.cos(.41*i1)
    i2,j2 = np.indices((7,4))+1
    r2 = r1[:,None] + .09*j2 + .03*i2*j2
    c2 = r2 + 1j*np.cos(.19*(i2+3*j2))
    i3,j3,k3 = np.indices((2,7,3))+1
    r3 = r1[j3-1]+.11*k3+.04*i3
    c3 = r3+1j*np.cos(.23*(i3+j3+2*k3))
    rh=np.array([.08,.17,.51,.17,.08])
    h=rh+1j*np.array([.04,-.08,.13,-.03,.06])
    return {(typ,rank): val for typ,rank,val in (
        ('r',1,r1),('c',1,c1),('r',2,r2),('c',2,c2),
        ('r',3,r3),('c',3,c3))},rh,h


def reference(name, values, rh, h):
    if name.startswith('ctaps_') or name.startswith('realtaps_'):
        match=re.fullmatch(r'(ctaps|realtaps)_(r|c)([123])_a([123])_([a-z]+)',name)
        if match is None: raise ValueError(name)
        kind,xtype,rank,axis,mode=match.groups()
        x=values[(xtype,int(rank))]
        return signal.resample_poly(x,3,2,window=h if kind=='ctaps' else rh,
                                    axis=int(axis)-1,padtype=mode)
    if name.startswith('outside_'):
        match=re.fullmatch(r'outside_(r|c)2_a([12])_constant',name)
        if match is None: raise ValueError(name)
        typ,axis=match.groups()
        return signal.resample_poly(values[(typ,2)],3,2,window=h,
                                    axis=int(axis)-1,padtype='constant',cval=.2-.12j)
    if name.startswith('up_'):
        match=re.fullmatch(r'up_(rr|rc|cr|cc)1_a1_([a-z]+)',name)
        if match is None: raise ValueError(name)
        xtype,mode=match.groups()
        return signal.upfirdn(rh if xtype[0]=='r' else h, values[(xtype[1],1)],
                              up=2,down=3,mode=mode)
    raise ValueError(name)


def main(executable):
    lines=subprocess.run([executable],check=True,text=True,capture_output=True).stdout.splitlines()
    values,rh,h=arrays()
    count=0
    value_count=0
    worst=0.0
    worst_case=''
    ptr=0
    while ptr<len(lines):
        fields=lines[ptr].split()
        if len(fields)<4 or fields[0]!='CASE':
            raise AssertionError(f'unexpected header line {ptr+1}: {lines[ptr]}')
        name=fields[1]
        rank=int(fields[2])
        shape=tuple(int(t) for t in fields[3:])
        if len(shape)!=rank: raise AssertionError((name,shape))
        length=int(np.prod(shape))
        block=lines[ptr+1:ptr+1+length]
        if len(block)!=length: raise AssertionError('short data '+name)
        actual=np.array([complex(*map(float,item.split())) for item in block]).reshape(shape,order='F')
        expected=reference(name,values,rh,h)
        if actual.shape != expected.shape: raise AssertionError((name,actual.shape,expected.shape))
        err=float(np.max(np.abs(actual-expected)))
        if err>worst: worst=err;worst_case=name
        if err>5e-11*max(1,float(np.max(np.abs(expected)))):
            raise AssertionError(f'{name} difference={err:.3e}')
        count+=1
        value_count+=length
        ptr+=length+1
    print(f'PASS: {count} SciPy cases, {value_count} numerical values, max error={worst:.3e} ({worst_case})')


if __name__=='__main__':
    main(sys.argv[1])
