#!/usr/bin/env python3
"""Compare rank-three lfilter along all axes with SciPy."""
from __future__ import annotations
import pathlib
import subprocess
import tempfile
import numpy as np
from scipy import signal

ROOT=pathlib.Path(__file__).resolve().parents[2]
RNG=np.random.default_rng(20261018)


def encode(x):
    return [f'{z.real:.17g} {z.imag:.17g}' for z in np.ravel(x,order='F')]


def main():
    lines,cases=[],[]
    for kind in (0,1):
        for case in range(150):
            n=RNG.integers(1,9,size=3).tolist()
            nb=int(RNG.integers(1,5))
            na=int(RNG.integers(1,5))
            axis=case%3
            haszi=case%2
            b=0.2*RNG.normal(size=nb)
            a=0.1*RNG.normal(size=na)
            a[0]=1.0+RNG.random()
            x=RNG.normal(size=n)
            if kind:
                b=b+0.1j*RNG.normal(size=nb)
                a=a+0.1j*RNG.normal(size=na)
                x=x+0.5j*RNG.normal(size=n)
            sh=n.copy()
            sh[axis]=max(nb,na)-1
            zi=np.zeros(sh,dtype=x.dtype) if not haszi else RNG.normal(size=sh)*0.1
            if kind and haszi:
                zi=zi+0.1j*RNG.normal(size=sh)
            lines.append(' '.join(map(str,[kind,*n,nb,na,axis+1,haszi])))
            lines.extend(encode(b))
            lines.extend(encode(a))
            lines.extend(encode(x))
            if haszi:
                lines.extend(encode(zi))
            y,zf=signal.lfilter(b,a,x,axis=axis,zi=zi)
            cases.append((y,zf))
    with tempfile.TemporaryDirectory() as tmp:
        w=pathlib.Path(tmp)
        cmd=['gfortran','-std=f2018','-O2','-Werror=implicit-interface',f'-I{w}',f'-J{w}']
        cmd += [str(ROOT/'src'/f'{src}.f90') for src in ('kinds','signal_lfilter')]
        cmd += [str(ROOT/'tools/signal_reference/lfilter_rank3_driver.f90'),'-o',str(w/'probe')]
        subprocess.run(cmd,cwd=w,check=True)
        proc=subprocess.run([str(w/'probe')],input='\n'.join(lines)+'\n',text=True,
                            capture_output=True,check=True,timeout=120)
        output=iter(proc.stdout.splitlines())
        worst=0.0
        for i,(y,zf) in enumerate(cases):
            status,*sizes=map(int,next(output).split())
            if status!=0 or list(y.shape)+list(zf.shape)!=sizes:
                raise AssertionError(f'case {i}: bad status or shapes {status}, {sizes}')
            def read_array(shape):
                return np.array([complex(*map(float,next(output).split()))
                    for _ in range(np.prod(shape))]).reshape(shape,order='F')
            fy,fz=read_array(y.shape),read_array(zf.shape)
            diff=max(np.max(np.abs(fy-y)),np.max(np.abs(fz-zf)) if zf.size else 0)
            worst=max(worst,diff)
            if not np.allclose(fy,y,rtol=3e-12,atol=3e-12) or not np.allclose(fz,zf,rtol=3e-12,atol=3e-12):
                raise AssertionError(f'case {i}: max difference {diff}')
        print(f'PASS: {len(cases)} rank-three lfilter cases; max abs difference {worst:.3g}')

if __name__=='__main__':
    main()
