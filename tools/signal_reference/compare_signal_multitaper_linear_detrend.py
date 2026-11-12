#!/usr/bin/env python3
"""Cross-check linear-detrended multitaper spectral matrices against SciPy.

SciPy provides the DPSS tapers and independent segment-linear detrending.
The spectral-matrix estimator and equal weighting are assembled explicitly
because SciPy does not have a corresponding multitaper spectral-matrix API.
Usage: python compare_signal_multitaper_linear_detrend.py DRIVER_EXE
"""
import subprocess
import sys
import numpy as np
from scipy import signal


def reference(data, nw, seg, ov, nk, nfft, rate, is_complex):
    n, nc = data.shape
    tapers = signal.windows.dpss(seg,nw,Kmax=nk,sym=False,norm=2)
    starts = list(range(0,n-seg+1,seg-ov))
    nfreq = nfft if is_complex else nfft//2+1
    result = np.zeros((nfreq,len(starts),nc,nc),dtype=np.complex128)
    for frame, start in enumerate(starts):
        segment = signal.detrend(data[start:start+seg],axis=0,type='linear')
        spectrum = np.empty((nk,nfreq,nc),dtype=complex)
        for k in range(nk):
            taper = tapers[k,:,None]*segment
            transform = (np.fft.fft(taper,n=nfft,axis=0) if is_complex else
                         np.fft.rfft(taper,n=nfft,axis=0))
            spectrum[k] = transform/np.sqrt(rate*np.dot(tapers[k],tapers[k]))
        for c in range(nc):
            for d in range(nc):
                result[:,frame,c,d] = np.mean(np.conj(spectrum[:,:,c])*spectrum[:,:,d],axis=0)
    if not is_complex:
        result[1: -1 if nfft%2==0 else None] *= 2
    return result


def main():
    if len(sys.argv) != 2:
        raise SystemExit('usage: compare_signal_multitaper_linear_detrend.py DRIVER_EXE')
    rng=np.random.default_rng(20261026)
    cases=[]
    for k in range(64):
        n=int(rng.integers(20,73))
        seg=int(rng.integers(8,min(n,36)+1))
        ov=int(rng.integers(0,seg))
        nc=int(rng.integers(1,4))
        nk=int(rng.integers(1,min(6,seg)+1))
        nw=float(rng.uniform(1.2,min(3.5,seg/2-0.2)))
        nfft=seg+int(rng.integers(0,15))
        rate=float(rng.uniform(10,200))
        is_complex=k%2==1
        data=rng.normal(size=(n,nc))
        if is_complex:
            data=data+1j*rng.normal(size=(n,nc))
        t=np.arange(n)[:,None]
        data=data+(rng.normal(size=(1,nc))*20+rng.normal(size=(1,nc))*t)
        if is_complex:
            data=data+1j*(rng.normal(size=(1,nc))*10+rng.normal(size=(1,nc))*t)
        if k%11==0:
            data=data*0+2.+.025*t+(1j if is_complex else 0)*t
        cases.append((data,nw,seg,ov,nk,nfft,rate,is_complex))
    def pair(z):
        return f'{np.real(z):.17e} {np.imag(z):.17e}\n'
    lines=[str(len(cases))+'\n']
    for data,nw,seg,ov,nk,nfft,rate,is_complex in cases:
        n,nc=data.shape
        lines.append(f'{n} {nc} {seg} {ov} {nk} {nfft} {int(is_complex)} {nw:.17e} {rate:.17e}\n')
        lines.extend(pair(v) for ch in range(nc) for v in data[:,ch])
    process=subprocess.run([sys.argv[1]],input=''.join(lines),text=True,
                           capture_output=True,check=True)
    output=iter(process.stdout.splitlines())
    max_diff=0.
    total=0
    for k,(data,nw,seg,ov,nk,nfft,rate,is_complex) in enumerate(cases):
        expected=reference(data,nw,seg,ov,nk,nfft,rate,is_complex)
        code,nf,nframe,nc=map(int,next(output).split())
        assert code==0,(k,code)
        assert (nf,nframe,nc)==expected.shape[:3],(k,nf,nframe,nc,expected.shape)
        actual=np.zeros_like(expected)
        for frame in range(nframe):
            for j in range(nc):
                for i in range(nc):
                    for f in range(nf):
                        actual[f,frame,i,j]=complex(*map(float,next(output).split()))
        error=float(np.max(abs(actual-expected)))
        max_diff=max(max_diff,error)
        np.testing.assert_allclose(actual,expected,rtol=1e-8,atol=3e-9,
                                   err_msg=f'case {k}')
        total+=actual.size
    print(f'PASS {len(cases)} SciPy multitaper linear-detrend cases, '
          f'{total} spectral matrix coefficients, max abs={max_diff:.3e}')


if __name__=='__main__':
    main()
