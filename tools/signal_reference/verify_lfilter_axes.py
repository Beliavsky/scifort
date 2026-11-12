#!/usr/bin/env python3
"""Compare rank-two real/complex lfilter and its final states against SciPy."""
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
        for case in range(200):
            nx,ny=RNG.integers(1,15,size=2)
            na=int(RNG.integers(1,6))
            nb=int(RNG.integers(1,6))
            axis=case%2
            haszi=(case//2)%2
            b=RNG.standard_normal(nb)*0.24
            a=RNG.standard_normal(na)*0.15
            a[0]=1.0+RNG.random()
            x=RNG.standard_normal((nx,ny))
            if kind:
                b=b+0.1j*RNG.standard_normal(nb)
                a=a+0.1j*RNG.standard_normal(na)
                x=x+1j*RNG.standard_normal((nx,ny))
            zi=None
            if haszi:
                shape=[nx,ny]
                shape[axis]=max(na,nb)-1
                zi=0.15*RNG.standard_normal(shape)
                if kind:
                    zi=zi+0.1j*RNG.standard_normal(shape)
            lines.append(f'{kind} {nx} {ny} {nb} {na} {axis+1} {haszi}')
            lines.extend(encode(b))
            lines.extend(encode(a))
            lines.extend(encode(x))
            if zi is not None:
                lines.extend(encode(zi))
            if zi is None:
                y,zf=signal.lfilter(b,a,x,axis=axis,zi=np.zeros(
                    [max(na,nb)-1 if i==axis else n for i,n in enumerate((nx,ny))],dtype=x.dtype))
            else:
                y,zf=signal.lfilter(b,a,x,axis=axis,zi=zi)
            cases.append((y,zf,axis))
    with tempfile.TemporaryDirectory() as tmp:
        work=pathlib.Path(tmp)
        cmd=['gfortran','-std=f2018','-O2','-Werror=implicit-interface',f'-I{work}',f'-J{work}']
        cmd += [str(ROOT/'src'/f'{src}.f90') for src in ('kinds','signal_lfilter')]
        cmd += [str(ROOT/'tools/signal_reference/lfilter_axes_driver.f90'),'-o',str(work/'probe')]
        subprocess.run(cmd,cwd=work,check=True)
        run=subprocess.run([str(work/'probe')],input='\n'.join(lines)+'\n',
                           text=True,capture_output=True,check=True,timeout=60)
        output=iter(run.stdout.splitlines())
        worst=0.0
        for case,(y,zf,axis) in enumerate(cases):
            status,*shapes=map(int,next(output).split())
            if status!=0 or list(y.shape)+list(zf.shape)!=shapes:
                raise AssertionError(f'case {case}: status={status} shapes={shapes}, '
                                     f'expected {list(y.shape)+list(zf.shape)}')
            def read_array(shape):
                return np.array([complex(*map(float,next(output).split()))
                    for _ in range(np.prod(shape))]).reshape(shape,order='F')
            got=read_array(y.shape)
            got_zf=read_array(zf.shape)
            diff=max(np.max(np.abs(got-y)),np.max(np.abs(got_zf-zf)) if zf.size else 0)
            worst=max(worst,diff)
            if not np.allclose(got,y,atol=3e-12,rtol=3e-12) or not np.allclose(got_zf,zf,atol=3e-12,rtol=3e-12):
                raise AssertionError(f'case {case}, axis={axis}, error={diff}')
        print(f'PASS: {len(cases)} rank-2 lfilter cases; max abs error {worst:.3g}')

if __name__=='__main__':
    main()
