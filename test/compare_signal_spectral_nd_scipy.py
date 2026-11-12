#!/usr/bin/env python3
"""SciPy 1.17 independent references for rank-independent spectral estimators.

Usage: python test/compare_signal_spectral_nd_scipy.py /path/to/fortran/exe
Generates deterministic cases, runs the Fortran driver, and compares spectra.
"""
import itertools
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np
from scipy import signal

DIMS = {2:(12,8),3:(12,8,10),6:(1,12,1,8,1,10),15:(1,)*12+(12,8,10)}
AXES = {2:(1,2),3:(1,2,3),6:(2,4,6),15:(13,14,15)}
WINDOWS = ('boxcar','hann','hamming')
DETREND = ('constant','linear','none')
SCALING = ('density','spectrum')
AVERAGE = ('mean','median')


def input_arrays(cid, dims, xt, yt):
    ix=np.arange(1,np.prod(dims)+1,dtype=np.float64)
    a=np.sin(0.071*ix+0.013*cid)+0.2*np.cos(0.047*ix-0.011*cid)
    b=np.cos(0.043*ix-0.021*cid)+0.3*np.sin(0.059*ix+0.009*cid)
    if xt:
        a=a+1j*0.22*np.cos(0.081*ix+0.017*cid)
    if yt:
        b=b+1j*0.17*np.sin(0.065*ix-0.012*cid)
    return a.reshape(dims,order='F'),b.reshape(dims,order='F')


def scalar_spectra(case):
    cid, rank, meth, xt, yt, ax, wi, di, si, ai=case
    dims=DIMS[rank]
    x,y=input_arrays(cid,dims,xt,yt)
    axis=ax-1
    # Flatten to independent channels and use SciPy's public API on the last
    # axis, avoiding SciPy 1.17 high-rank padding/shape limitations.
    xx=np.moveaxis(x,axis,-1)
    yy=np.moveaxis(y,axis,-1)
    lead=xx.shape[:-1]
    vx=xx.reshape(-1,xx.shape[-1])
    vy=yy.reshape(-1,yy.shape[-1])
    kwargs=dict(fs=7.0,window=WINDOWS[wi],nfft=16,
                detrend=(False if di == 2 else DETREND[di]),axis=-1)
    if meth >= 3 and (xt or yt):
        # SciPy 1.17 decides sidedness from the first array and may reject
        # mixed dtypes. Explicit promotion requests the correct full spectrum.
        vx=np.asarray(vx,dtype=np.complex128)
        vy=np.asarray(vy,dtype=np.complex128)
    if meth==1:
        f,p=signal.periodogram(vx,scaling=SCALING[si],**kwargs)
    elif meth==2:
        f,p=signal.welch(vx,nperseg=6,noverlap=3,average=AVERAGE[ai],
                         scaling=SCALING[si],**kwargs)
    elif meth==3:
        f,p=signal.csd(vx,vy,nperseg=6,noverlap=3,average=AVERAGE[ai],
                       scaling=SCALING[si],**kwargs)
    elif meth==4:
        f,p=signal.coherence(vx,vy,nperseg=6,noverlap=3,**kwargs)
    else:
        raise ValueError(meth)
    p=np.moveaxis(p.reshape(lead+(len(f),)),-1,axis)
    mask=np.ones_like(p,dtype=bool)
    if meth==4 and di!=2:
        # Detrending makes DC power roundoff-sized. Coherence there is an
        # ill-conditioned 0/0 ratio and differs across FFT implementations.
        selector=[slice(None)]*rank
        selector[axis]=0
        mask[tuple(selector)]=False
    return np.asarray(f),np.asarray(p).flatten(order='F'),mask.flatten(order='F')


def main(executable):
    cases=[]
    cid=0
    for rank, axes in AXES.items():
        for ax,meth in itertools.product(axes,range(1,5)):
            pairs=((0,0),(1,0)) if meth<3 else ((0,0),(1,0),(0,1),(1,1))
            for xt,yt in pairs:
                for wi in range(3):
                    # Change detrending, scaling and averaging independently.
                    cid+=1
                    cases.append((cid,rank,meth,xt,yt,ax,wi,(cid//2)%3,cid%2,(cid//3)%2))
    with tempfile.TemporaryDirectory(prefix='scifort_spectral_nd_') as tmp:
        casefile=Path(tmp)/'cases.txt'
        outfile=Path(tmp)/'out.bin'
        casefile.write_text(''.join(' '.join(map(str,c))+'\n' for c in cases))
        subprocess.run([executable,str(casefile),str(outfile)],check=True)
        maxerror=0.0
        values=0
        with outfile.open('rb') as handle:
            for case in cases:
                freq, expected, mask=scalar_spectra(case)
                status=np.fromfile(handle,dtype=np.int32,count=1)
                actualf=np.fromfile(handle,dtype=np.float64,count=len(freq))
                actual=np.fromfile(handle,dtype=np.complex128,count=len(expected))
                if len(status)!=1 or status[0]!=0 or len(actual)!=len(expected):
                    raise AssertionError(f'case {case}: missing or failed data')
                maxerror=max(maxerror,float(np.max(np.abs(freq-actualf))))
                diff=np.abs(actual-expected)
                ok=np.isclose(actual,expected,atol=4e-10,rtol=4e-10,equal_nan=True)
                ok[~mask]=True
                if not np.all(ok):
                    k=np.flatnonzero(~ok)[0]
                    raise AssertionError(f'case {case}, bin {k}: SciPy {expected[k]} vs Fortran {actual[k]}, error={diff[k]}')
                err=np.nanmax(diff[mask]) if np.isfinite(diff[mask]).any() else 0.0
                maxerror=max(maxerror,float(err))
                values+=len(expected)
            if handle.read(1):
                raise AssertionError('unexpected extra Fortran output')
    print(f'PASS: {len(cases)} SciPy cases, {values} spectral values, '
          f'max abs difference {maxerror:.3g}')

if __name__=='__main__':
    if len(sys.argv)!=2:
        raise SystemExit('usage: compare_signal_spectral_nd_scipy.py FORTRAN_EXE')
    main(sys.argv[1])
