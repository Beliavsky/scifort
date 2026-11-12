#!/usr/bin/env python3
"""Check rank-independent STFT and spectrogram against SciPy 1.17.

Usage: python test/compare_signal_timefreq_nd_scipy.py PATH_TO_FORTRAN_REFERENCE
"""
import itertools
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy import signal

DIMS = {
    2: (12,8),
    3: (12,8,10),
    6: (1,12,1,8,1,10),
    14: (1,)*11+(12,8,10),
}
AXES = {2:(1,2),3:(1,2,3),6:(2,4,6),14:(12,13,14)}
WINDOWS = ('boxcar','hann','tukey')
DETREND = ('constant','linear',False)
SCALING = ('density','spectrum')
MODES = {2:'psd',3:'complex',4:'magnitude',5:'angle',6:'phase'}


def input_data(cid,dims,xtype):
    ix=np.arange(1,np.prod(dims)+1,dtype=np.float64)
    x=np.sin(0.071*ix+0.013*cid)+0.2*np.cos(0.047*ix-0.011*cid)
    if xtype:
        x=x+1j*0.22*np.cos(0.081*ix+0.017*cid)
    return x.reshape(dims,order='F')


def scipy_reference(case, mode_override=None):
    cid,rk,meth,xtype,ax,wi,di,si,bi,pa=case
    d=DIMS[rk]
    x=input_data(cid,d,xtype)
    # Perform independent one-dimensional SciPy computations, avoiding
    # SciPy 1.17 higher-rank padding bugs in the legacy STFT interface.
    moved=np.moveaxis(x,ax-1,-1)
    lead=moved.shape[:-1]
    flat=moved.reshape(-1,moved.shape[-1])
    vals=[]
    kwargs=dict(fs=7.0,window=WINDOWS[wi],nperseg=6,noverlap=2,
                nfft=16,detrend=DETREND[di])
    if meth==1:
        kwargs.update(boundary=('zeros' if bi else None),padded=bool(pa),
                      scaling=('psd' if si==0 else 'spectrum'))
    else:
        kwargs.update(scaling=SCALING[si],mode=(mode_override or MODES[meth]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for lane in flat:
            if meth==1:
                f,t,part=signal.stft(lane,**kwargs)
            else:
                f,t,part=signal.spectrogram(lane,**kwargs)
            vals.append(part)
    values=np.asarray(vals).reshape(lead+(len(f),len(t)))
    values=np.moveaxis(values,-2,ax-1)
    return f,t,values.flatten(order='F')


def main(executable):
    cases=[]
    for rk,axes in AXES.items():
        for ax,meth,xt,wi in itertools.product(axes,range(1,7),(0,1),range(3)):
            cid=len(cases)+1
            cases.append((cid,rk,meth,xt,ax,wi,(cid//2)%3,cid%2,(cid//4)%2,(cid//3)%2))
    with tempfile.TemporaryDirectory(prefix='scifort_timefreq_nd_') as tmp:
        casefile=Path(tmp)/'cases.txt'
        outfile=Path(tmp)/'values.bin'
        casefile.write_text(''.join(' '.join(map(str,c))+'\n' for c in cases))
        subprocess.run([executable,str(casefile),str(outfile)],check=True)
        values=0
        maxdiff=0.0
        with outfile.open('rb') as fd:
            for case in cases:
                freqs,times,expected=scipy_reference(case)
                status=np.fromfile(fd,dtype=np.int32,count=1)
                actualf=np.fromfile(fd,dtype=np.float64,count=len(freqs))
                actualt=np.fromfile(fd,dtype=np.float64,count=len(times))
                actual=np.fromfile(fd,dtype=np.complex128,count=len(expected))
                if len(status)!=1 or status[0]!=0 or len(actual)!=len(expected):
                    raise AssertionError(f'case {case}: missing output or invalid status')
                ferr=np.max(np.abs(freqs-actualf))
                terr=np.max(np.abs(times-actualt))
                diff=np.abs(actual-expected)
                if case[2]==6:
                    # Unwrapped phase can differ by integral 2*pi turns
                    # after a numerically undefined zero-power frequency bin.
                    diff=np.abs(np.remainder((actual.real-expected.real)+np.pi,
                                             2*np.pi)-np.pi)
                mask=np.ones_like(expected,dtype=bool)
                if case[2] in (5,6):
                    # Phase is not numerically defined where the STFT value
                    # vanishes (e.g. detrended DC for a boxcar window).
                    # Different FFT sum order can change 0/pi at such bins.
                    _,_,complex_values=scipy_reference(case,mode_override='complex')
                    scale=np.max(np.abs(complex_values))
                    mask=np.abs(complex_values)>1e-12*scale
                maxdiff=max(maxdiff,float(ferr),float(terr),
                            float(np.max(diff[mask])) if mask.any() else 0.0)
                if case[2]==6:
                    ok=np.asarray(diff<3e-10)
                else:
                    ok=np.isclose(actual,expected,atol=3e-10,rtol=2e-10,equal_nan=True)
                ok[~mask]=True
                if not np.all(ok):
                    bad=np.flatnonzero(~ok)[0]
                    raise AssertionError(f'case {case}, coefficient {bad}: '
                                         f'SciPy={expected[bad]}, Fortran={actual[bad]}, '
                                         f'difference={diff[bad]}')
                if not np.allclose(actualf,freqs,rtol=1e-12,atol=1e-12):
                    raise AssertionError(f'case {case}: bad frequencies')
                if not np.allclose(actualt,times,rtol=1e-12,atol=1e-12):
                    raise AssertionError(f'case {case}: bad time coordinates')
                values+=len(expected)
            if fd.read(1):
                raise AssertionError('unexpected trailing Fortran output')
    print(f'PASS: {len(cases)} SciPy STFT/spectrogram cases, {values} coefficients, '
          f'max absolute error {maxdiff:.3g}')


if __name__=='__main__':
    if len(sys.argv)!=2:
        raise SystemExit('usage: compare_signal_timefreq_nd_scipy.py FORTRAN_REFERENCE_EXE')
    main(sys.argv[1])
