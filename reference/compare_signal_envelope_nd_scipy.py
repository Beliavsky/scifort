#!/usr/bin/env python3
"""Deterministic rank-four envelope_into comparisons against SciPy 1.17.

Usage: python reference/compare_signal_envelope_nd_scipy.py path/to/driver
"""
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.signal import envelope


def main() -> None:
    executable = Path(sys.argv[1]) if len(sys.argv)>1 else Path('build/signal_reference/ref_envelope_nd')
    rng = np.random.default_rng(773229)
    cases = []
    for t in range(360):
        shp = tuple(int(v) for v in rng.integers(1, 5, size=4))
        ax = int(rng.integers(0,4))
        n=shp[ax]
        nout=int(rng.choice([1, 2, n, n+1, 2*n+1, 13]))
        lo=int(rng.integers(-(n//2), (n+1)//2))
        hi=int(rng.integers(lo+1,(n+1)//2+1))
        cx=int(rng.integers(0,2))
        squared=int(rng.integers(0,2))
        res=int(rng.integers(0,3))
        z=rng.normal(size=shp)+1j*rng.normal(size=shp)
        if cx == 0:
            z = z.real
        cases.append((shp,ax,nout,lo,hi,cx,squared,res,z))
    lines = [str(len(cases))]
    for shp,ax,nout,lo,hi,cx,squared,res,z in cases:
        lines.append(' '.join(map(str,(*shp,ax+1,nout,lo,hi,cx,squared,res))))
        lines.extend(f'{v.real:.17g} {v.imag:.17g}' for v in np.ravel(z,order='F'))
    proc = subprocess.run([str(executable.resolve())], input='\n'.join(lines)+'\n',
                          text=True,capture_output=True,check=True)
    output=iter(proc.stdout.splitlines())
    max_err=0.0
    total=0
    for index,(shp,ax,nout,lo,hi,cx,squared,res,z) in enumerate(cases):
        status,*dims = map(int,next(output).split())
        if status != 0:
            raise AssertionError((index,'Fortran status',status))
        out_shape=tuple(dims)
        if out_shape[ax]!=nout:
            raise AssertionError((index,'Fortran shape',out_shape))
        got=np.array([list(map(float,next(output).split())) for _ in range(np.prod(out_shape))])
        got=np.stack([got[:,j].reshape(out_shape,order='F') for j in range(3)],axis=-1)
        reference=envelope(z,bp_in=(lo,hi),n_out=nout,squared=bool(squared),
                           residual=[None,'lowpass','all'][res],axis=ax)
        ref_env = reference if res == 0 else reference[0]
        ref_res = np.zeros_like(ref_env) if res == 0 else reference[1]
        ref_env=np.asarray(ref_env)
        ref_res=np.asarray(ref_res)
        err=max(float(np.max(np.abs(got[...,0]-ref_env))),
                float(np.max(np.abs(got[...,1]+1j*got[...,2]-ref_res))))
        max_err=max(max_err,err)
        total+=np.prod(out_shape)
        if err>3e-10:
            raise AssertionError((index,shp,ax,nout,lo,hi,cx,squared,res,err))
    print(f'PASS: {len(cases)} SciPy envelope cases, {total} values, max_abs={max_err:.3e}')


if __name__=='__main__':
    main()
