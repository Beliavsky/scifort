#!/usr/bin/env python3
"""Compile spectral driver and compare deterministic randomized cases with SciPy."""
from __future__ import annotations

import pathlib
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = pathlib.Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261015)
PROGRAMS = ("kinds.f90", "signal_dpss.f90", "signal_windows_extra.f90", "signal_spectral.f90", "signal_fir_advanced.f90", "signal_spectral_complex.f90", "signal_convolve.f90", "signal_convolve2d.f90", "signal_fftconvolve.f90", "signal_detrend.f90", "signal_peaks.f90", "signal_savgol.f90", "signal_lfilter.f90", "signal_sos.f90", "signal_filtfilt.f90", "signal_filtfilt_nd.f90", "signal_analog_design.f90", "signal_design.f90", "signal_elliptic.f90", "signal_analog_sos.f90", "signal_iir_extra.f90", "signal_fourier_resample.f90", "signal_resample.f90", "signal_analog.f90", "signal_zpk_conversion.f90", "signal_state_space.f90", "signal_cont2discrete.f90", "signal_discrete_simulation.f90", "signal_continuous_simulation.f90", "signal_lti_frequency.f90", "signal_partial_fraction.f90", "signal_sequences.f90", "signal_hilbert.f90", "signal_waveforms.f90", "signal_periodic.f90", "signal_irregular.f90", "signal_remez.f90", "signal_czt.f90", "signal.f90")


def run(exe, method, x, y=None, *, window='hann', nw=16, no=8, nf=16, fs=13.,
        trend='constant', scaling='density', average='mean', boundary='none', padded=False, flag=0):
    if y is None:
        y = np.array([], dtype=float)
    data = '\n'.join([method, f'{len(x)} {len(y)} {nw} {no} {nf} {fs:.17g} {flag or int(padded)}',
                      window, trend, scaling, average, boundary,
                      *(f'{z:.17g}' for z in x), *(f'{z:.17g}' for z in y)]) + '\n'
    p = subprocess.run([str(exe)], input=data, capture_output=True, text=True, timeout=10)
    if p.returncode:
        raise RuntimeError(f'{method}: {p.stderr}\n{data}')
    lines = p.stdout.strip().splitlines()
    status, m, t = map(int, lines[0].split())
    if status:
        raise ValueError(f'{method} status={status}: {data}')
    arr = np.array([[float(q.replace('D','E')) for q in line.split()] for line in lines[1:]])
    if arr.shape[0] != m*t:
        raise AssertionError((method, arr.shape, (m,t), lines[:10]))
    if method in ('stft', 'spectrogram'):
        f = arr[:m,0]
        times = arr[::m,1]
        if method == 'stft':
            out = (arr[:,2] + 1j*arr[:,3]).reshape((t,m)).T
        else:
            out = arr[:,2].reshape((t,m)).T
        return f, times, out
    if method == 'csd':
        return arr[:,0], arr[:,1]+1j*arr[:,2]
    return arr[:,0], arr[:,1]


def eq(a, b, label):
    if a.shape != b.shape or not np.allclose(a,b,rtol=1e-9,atol=5e-10,equal_nan=True):
        diff = np.max(np.abs(a-b)) if a.shape == b.shape else None
        raise AssertionError(f'{label}: {a.shape} vs {b.shape}, err={diff}\ngot={a!r}\nref={b!r}')


def verify(exe):
    count = 0
    for kind in ('get_window', 'periodogram', 'welch', 'csd', 'coherence', 'stft', 'spectrogram'):
        for case in range(200):
            n = int(RNG.integers(10, 110))
            x = RNG.standard_normal(n) + 0.1*np.arange(n)
            y = RNG.standard_normal(n) - 0.05*np.arange(n)
            window = ('hann', 'boxcar', 'hamming', 'blackman', 'bartlett', 'tukey', 'flattop')[case%7]
            scipywindow = ('tukey',0.25) if window == 'tukey' else window
            nw = int(RNG.integers(1,min(n,30)+1))
            no = int(RNG.integers(0,nw))
            nf = nw + int(RNG.integers(0,12))
            fs = float(RNG.uniform(0.5,200))
            trend = ('constant', 'linear', 'none')[case%3]
            scaling = ('density', 'spectrum')[case%2]
            average = ('mean', 'median')[case%2]
            try:
                if kind == 'get_window':
                    mode = (case%2==0)
                    f,p = run(exe,kind,x,window=window,nw=nw,fs=fs,flag=int(mode))
                    ref = signal.get_window(scipywindow,nw,fftbins=mode)
                    eq(p,ref,f'{kind} {case}')
                elif kind == 'periodogram':
                    f,p = run(exe,kind,x,window=window,nw=nw,nf=nf,trend=trend,scaling=scaling,fs=fs)
                    rf,rp = signal.periodogram(x,fs=fs,window=scipywindow,nfft=nf,detrend=trend if trend!='none' else False,scaling=scaling)
                    eq(f,rf,f'{kind} f {case}');eq(p,rp,f'{kind} p {case}')
                elif kind == 'welch':
                    f,p = run(exe,kind,x,window=window,nw=nw,no=no,nf=nf,trend=trend,scaling=scaling,fs=fs,average=average)
                    rf,rp = signal.welch(x,fs=fs,window=scipywindow,nperseg=nw,noverlap=no,nfft=nf,detrend=trend if trend!='none' else False,scaling=scaling,average=average)
                    eq(f,rf,f'{kind} f {case}');eq(p,rp,f'{kind} p {case}')
                elif kind == 'csd':
                    f,p = run(exe,kind,x,y,window=window,nw=nw,no=no,nf=nf,trend=trend,scaling=scaling,fs=fs,average=average)
                    rf,rp = signal.csd(x,y,fs=fs,window=scipywindow,nperseg=nw,noverlap=no,nfft=nf,detrend=trend if trend!='none' else False,scaling=scaling,average=average)
                    eq(f,rf,f'{kind} f {case}');eq(p,rp,f'{kind} p {case}')
                elif kind == 'coherence':
                    f,p = run(exe,kind,x,y,window=window,nw=nw,no=no,nf=nf,trend=trend,fs=fs)
                    rf,rp = signal.coherence(x,y,fs=fs,window=scipywindow,nperseg=nw,noverlap=no,nfft=nf,detrend=trend if trend!='none' else False)
                    ap = signal.welch(x,fs=fs,window=scipywindow,nperseg=nw,noverlap=no,nfft=nf,detrend=trend if trend!='none' else False)[1]
                    bp = signal.welch(y,fs=fs,window=scipywindow,nperseg=nw,noverlap=no,nfft=nf,detrend=trend if trend!='none' else False)[1]
                    valid = (ap > np.max(ap)*1e-18) & (bp > np.max(bp)*1e-18)
                    if nw <= 2 and trend == "linear":
                        valid[:] = False
                    eq(f,rf,f'{kind} f {case}');eq(p[valid],rp[valid],f'{kind} p {case}')
                elif kind == 'stft':
                    boundary = ('zeros','none')[case%2]
                    pad = bool(case%3 == 0)
                    scale = ('spectrum','psd')[case%2]
                    f,t,z = run(exe,kind,x,window=window,nw=nw,no=no,nf=nf,trend=trend,scaling=scale,fs=fs,boundary=boundary,padded=pad)
                    rf,rt,rz = signal.stft(x,fs=fs,window=scipywindow,nperseg=nw,noverlap=no,nfft=nf,detrend=trend if trend!='none' else False,boundary=None if boundary=='none' else 'zeros',padded=pad,scaling=scale)
                    eq(f,rf,f'{kind} f {case}');eq(t,rt,f'{kind} t {case}');eq(z,rz,f'{kind} z {case}')
                elif kind == 'spectrogram':
                    f,t,p = run(exe,kind,x,window=window,nw=nw,no=no,nf=nf,trend=trend,scaling=scaling,fs=fs)
                    rf,rt,rp = signal.spectrogram(x,fs=fs,window=scipywindow,nperseg=nw,noverlap=no,nfft=nf,detrend=trend if trend!='none' else False,scaling=scaling)
                    eq(f,rf,f'{kind} f {case}');eq(t,rt,f'{kind} t {case}');eq(p,rp,f'{kind} p {case}')
                count+=1
            except Exception:
                print('Failed',kind,'case',case,'n',n,'nw',nw,'no',no,'nf',nf,'window',window,'trend',trend,'scaling',scaling,'average',average)
                raise
        print(f'{kind}: 200 matched')
    for method in ('csd', 'coherence'):
        for case in range(80):
            n = int(RNG.integers(16,128))
            m = int(RNG.integers(16,128))
            x = RNG.normal(size=n)
            y = RNG.normal(size=m)
            nw = int(RNG.integers(4,min(n,m,35)))
            no = int(RNG.integers(0,nw))
            nf = nw + int(RNG.integers(0,15))
            f, v = run(exe,method,x,y,nw=nw,no=no,nf=nf,fs=13.)
            if method == 'csd':
                rf, rv = signal.csd(x,y,nperseg=nw,noverlap=no,nfft=nf,fs=13.)
            else:
                rf, rv = signal.coherence(x,y,nperseg=nw,noverlap=no,nfft=nf,fs=13.)
            eq(f,rf,f'{method} unequal frequencies')
            eq(v,rv,f'{method} unequal spectra')
            count += 1
        print(f'{method}: 80 unequal-length cases matched')
    for method in ('welch', 'stft'):
        for case in range(40):
            n = int(RNG.integers(280,600))
            x = RNG.normal(size=n)
            nw = [64,128,256,512][case%4]
            if nw > n:
                nw = 256
            no = nw//2
            nf = nw
            if method == 'welch':
                f, p = run(exe,method,x,nw=nw,no=no,nf=nf,fs=21.)
                rf,rp = signal.welch(x,nperseg=nw,noverlap=no,nfft=nf,fs=21.)
                eq(f,rf,'radix2 Welch f');eq(p,rp,'radix2 Welch PSD')
            else:
                f,t,p = run(exe,method,x,nw=nw,no=no,nf=nf,fs=21.,trend='none',scaling='spectrum',boundary='zeros',padded=True)
                rf,rt,rp = signal.stft(x,nperseg=nw,noverlap=no,nfft=nf,fs=21.)
                eq(f,rf,'radix2 STFT f');eq(t,rt,'radix2 STFT t');eq(p,rp,'radix2 STFT Z')
            count += 1
        print(f'{method}: 40 large power-of-two FFT cases matched')
    return count


def main():
    with tempfile.TemporaryDirectory() as tmp:
        temp = pathlib.Path(tmp)
        args=['gfortran','-std=f2018','-O2','-J'+str(temp),'-I'+str(temp)]
        sources=[str(ROOT/'src'/p) for p in PROGRAMS]
        sources.append(str(ROOT/'tools/signal_reference/spectral_driver.f90'))
        subprocess.run(args+sources+['-o',str(temp/'spectral_driver')],check=True,cwd=temp)
        n=verify(temp/'spectral_driver')
        print(f'{n} spectral comparisons passed against SciPy {signal.__version__ if hasattr(signal,"__version__") else "installed"}')


if __name__=='__main__':
    main()
