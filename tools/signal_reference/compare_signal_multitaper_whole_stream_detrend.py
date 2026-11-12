#!/usr/bin/env python3
"""Independently verify SciFort whole-record and streaming multitaper linear detrending.

SciPy 1.17 supplies DPSS tapers and least-squares linear detrending; NumPy
performs the FFT and weighted cross-spectral matrix calculation. SciPy has no
matching multitaper matrix estimator, so the reference combines its primitives.
Usage: python compare_signal_multitaper_whole_stream_detrend.py DRIVER_EXE
"""
import subprocess
import sys
import numpy as np
from scipy import signal


def spectral_matrix(x, nw, k, nfft, rate, complex_input, mode):
    n,nc=x.shape
    tapers,ratios=signal.windows.dpss(n,nw,Kmax=k,sym=False,norm=2,return_ratios=True)
    clean=signal.detrend(x,type='linear',axis=0)
    if complex_input:
        ft=np.fft.fft(clean[None,:,:]*tapers[:,:,None],n=nfft,axis=1)
    else:
        ft=np.fft.rfft(clean[None,:,:]*tapers[:,:,None],n=nfft,axis=1)
    q=np.ones(k) if mode=='equal' else ratios
    weights=q/q.sum()
    ft=ft/np.sqrt(rate*np.sum(tapers**2,axis=1))[:,None,None]
    mat=np.einsum('k,kfc,kfd->fcd',weights,ft.conj(),ft,optimize=True)
    if not complex_input:
        mat[1:-1 if nfft%2==0 else None]*=2.
    return mat


def main():
    if len(sys.argv)!=2:
        raise SystemExit('usage: compare_signal_multitaper_whole_stream_detrend.py DRIVER_EXE')
    rng=np.random.default_rng(20261026)
    cases=[]
    lines=[]
    for i in range(96):
        n=int(rng.integers(15,83))
        nc=int(rng.integers(1,4))
        seg=int(rng.integers(8,min(36,n)+1))
        ov=int(rng.integers(0,seg))
        k=int(rng.integers(2,min(6,seg)+1))
        nw=float(rng.uniform(1.2,min(3.5,seg/2-0.1)))
        nfft=n+int(rng.integers(0,17))
        complex_input=bool(i%2)
        mode=('equal','eigen')[i%2 if i%7 else (i//7)%2]
        rate=float(rng.uniform(1.,160.))
        x=rng.standard_normal((n,nc))
        if complex_input:
            x=x+1j*rng.standard_normal((n,nc))
        t=np.arange(n)[:,None]
        x+=rng.uniform(-12.,12.,(1,nc))+rng.uniform(-2.,2.,(1,nc))*t
        if complex_input:
            x+=1j*(rng.uniform(-12.,12.,(1,nc))+rng.uniform(-2.,2.,(1,nc))*t)
        if i%13==0:
            x=np.ones((n,nc))*8.3 + t*0.125
            if complex_input: x=x+1j*(-3.+t*.57)
        cases.append((x,nw,seg,ov,k,nfft,complex_input,mode,rate))
        lines.append(f'{n} {nc} {seg} {ov} {k} {nfft} {int(complex_input)} '
                     f'{1 if mode=="equal" else 2} {nw:.17g} {rate:.17g}')
        lines.extend(f'{v.real:.17g} {v.imag:.17g}' for c in range(nc) for v in x[:,c])
    out=subprocess.run([sys.argv[1]],input=f'{len(cases)}\n'+'\n'.join(lines)+'\n',
                       capture_output=True,text=True,check=True)
    rows=iter(out.stdout.splitlines())
    maximum=0.
    total=0
    for idx,(x,nw,seg,ov,k,nfft,iscomplex,mode,rate) in enumerate(cases):
        sfull,sinit,stf,spsd,sci,nf,nframe=map(int,next(rows).split())
        assert (sfull,sinit,stf,spsd,sci)==(0,0,0,0,0), (idx,sfull,sinit,stf,spsd,sci)
        expected=spectral_matrix(x,nw,k,nfft,rate,iscomplex,mode)
        assert nf==expected.shape[0], (idx,nf,expected.shape)
        coordinates=np.array([list(map(float,next(rows).split())) for _ in range(nf)])
        freq=np.fft.fftfreq(nfft,1./rate) if iscomplex else np.fft.rfftfreq(nfft,1./rate)
        np.testing.assert_allclose(coordinates[:,0],freq,rtol=1e-13,atol=1e-12)
        full=np.array([complex(*map(float,next(rows).split()))
                       for _ in range(nf*x.shape[1]*x.shape[1])])
        full=full.reshape(x.shape[1],x.shape[1],nf).transpose(2,1,0)
        np.testing.assert_allclose(full,expected,rtol=3e-8,atol=1e-9,err_msg=f'whole matrix {idx}')
        np.testing.assert_allclose(coordinates[:,1],expected[:,0,0].real,rtol=3e-8,atol=1e-9)
        np.testing.assert_allclose(coordinates[:,2],expected[:,0,0].real,rtol=3e-8,atol=1e-9)
        step=seg-ov
        positions=list(range(0,len(x)-seg+1,step))
        assert nframe==len(positions),(idx,nframe,len(positions))
        stream=np.array([complex(*map(float,next(rows).split()))
                         for _ in range(nframe*(seg+3 if iscomplex else (seg+3)//2+1)*x.shape[1]**2)])
        nfreq=(seg+3 if iscomplex else (seg+3)//2+1)
        stream=stream.reshape(nframe,x.shape[1],x.shape[1],nfreq).transpose(3,0,2,1)
        expected_frames=np.stack([spectral_matrix(x[start:start+seg],nw,k,seg+3,rate,iscomplex,mode)
                                  for start in positions],axis=1)
        np.testing.assert_allclose(stream,expected_frames,rtol=3e-8,atol=1e-9,err_msg=f'stream {idx}')
        maximum=max(maximum,float(np.max(abs(full-expected))),float(np.max(abs(stream-expected_frames))))
        total+=full.size+stream.size
    assert next(rows,None) is None, 'extra driver output'
    print(f'PASS: {len(cases)} independent SciPy comparisons, {total} cross-spectral coefficients')
    print(f'Maximum absolute difference {maximum:.3e}')

if __name__=='__main__':
    main()
