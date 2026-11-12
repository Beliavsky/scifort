#!/usr/bin/env python3
"""Compile SciFort complex spectral driver and compare to SciPy reference values."""
from __future__ import annotations

import pathlib
import subprocess
import tempfile
import warnings

import numpy as np
from scipy import signal

ROOT = pathlib.Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261015)
SOURCES = (
    'kinds.f90', 'signal_dpss.f90', 'signal_windows_extra.f90', 'signal_spectral.f90', 'signal_fir_advanced.f90', 'signal_spectral_complex.f90',
    'signal_convolve.f90', 'signal_convolve2d.f90', 'signal_fftconvolve.f90', 'signal_detrend.f90', 'signal_peaks.f90',
    'signal_savgol.f90', 'signal_lfilter.f90', 'signal_sos.f90', 'signal_filtfilt.f90', 'signal_filtfilt_nd.f90', 'signal_analog_design.f90', 'signal_design.f90', 'signal_elliptic.f90', 'signal_analog_sos.f90', 'signal_iir_extra.f90', 'signal_fourier_resample.f90', 'signal_resample.f90', 'signal_analog.f90', 'signal_zpk_conversion.f90', 'signal_state_space.f90', 'signal_cont2discrete.f90', 'signal_discrete_simulation.f90', 'signal_continuous_simulation.f90', 'signal_lti_frequency.f90', 'signal_partial_fraction.f90', 'signal_sequences.f90', 'signal_hilbert.f90', 'signal_waveforms.f90', 'signal_periodic.f90', 'signal_irregular.f90', 'signal_remez.f90', 'signal_czt.f90', 'signal.f90',
)


def driver(exe, method, x, y, *, fs, window, nw, no, nf, trend, scaling,
           average, boundary, padded):
    lines = [method, f'{len(x)} {len(y)} {nw} {no} {nf} {fs:.17g} {int(padded)}',
             window, trend, scaling, average, boundary]
    lines.extend(f'{z.real:.17g} {z.imag:.17g}' for z in x)
    lines.extend(f'{z.real:.17g} {z.imag:.17g}' for z in y)
    proc = subprocess.run([str(exe)], input='\n'.join(lines)+'\n', text=True,
                          capture_output=True, timeout=20)
    if proc.returncode:
        raise RuntimeError(f'{method}: {proc.stderr}')
    out = proc.stdout.strip().splitlines()
    status, m, t = map(int, out[0].split())
    if status:
        raise AssertionError(f'{method}: unexpected status {status}')
    values = np.array([[float(y) for y in s.split()] for s in out[1:]])
    assert values.shape[0] == m * t, (method, values.shape, m, t)
    if method in ('stft','spectrogram'):
        freqs = values[:m,0]
        times = values[::m,1]
        if method == 'stft':
            z = (values[:,2]+1j*values[:,3]).reshape(t,m).T
        else:
            z = values[:,2].reshape(t,m).T
        return freqs, times, z
    if method == 'csd':
        return values[:,0], values[:,1]+1j*values[:,2]
    return values[:,0], values[:,1]


def same(actual, ref, label):
    if actual.shape != ref.shape or not np.allclose(actual,ref,rtol=4e-9,atol=1e-9,equal_nan=True):
        err = np.max(np.abs(actual-ref)) if actual.shape == ref.shape else None
        raise AssertionError(f'{label}: shapes {actual.shape}/{ref.shape}; error {err}\n{actual}\n{ref}')


def verify(exe):
    count = 0
    for method in ('periodogram','welch','csd','coherence','stft','spectrogram'):
        for case in range(160):
            n = int(RNG.integers(15,112))
            x = (RNG.standard_normal(n)+1j*RNG.standard_normal(n)) * (0.5+case%5)
            ny = int(RNG.integers(15,112)) if method in ('csd','coherence') and case%4==0 else n
            y = RNG.standard_normal(ny)+1j*RNG.standard_normal(ny)
            nw = int(RNG.integers(3,min(n,ny,34)+1))
            no = int(RNG.integers(0,nw))
            nf = nw+int(RNG.integers(0,20))
            if method == 'periodogram':
                nf = n+int(RNG.integers(-4,16))
                nf = max(1,nf)
            window = ('boxcar','hann','hamming','blackman','bartlett','flattop','tukey')[case%7]
            swindow = ('tukey',0.25) if window == 'tukey' else window
            trend = ('constant','linear','none')[case%3]
            scaling = ('density','spectrum')[case%2]
            average = ('mean','median')[case%2]
            boundary = ('zeros','none')[case%2]
            padded = (case%3 == 0)
            if method == 'stft':
                scaling = ('psd','spectrum')[case%2]
            fs = float(RNG.uniform(0.8,53))
            reftrend = False if trend=='none' else trend
            try:
                result = driver(exe,method,x,y,fs=fs,window=window,nw=nw,no=no,nf=nf,
                                trend=trend,scaling=scaling,average=average,
                                boundary=boundary,padded=padded)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',UserWarning)
                    if method == 'periodogram':
                        ref = signal.periodogram(x,fs=fs,window=swindow,nfft=nf,
                                                 detrend=reftrend,scaling=scaling)
                    elif method == 'welch':
                        ref = signal.welch(x,fs=fs,window=swindow,nperseg=nw,
                                           noverlap=no,nfft=nf,detrend=reftrend,
                                           scaling=scaling,average=average)
                    elif method == 'csd':
                        ref = signal.csd(x,y,fs=fs,window=swindow,nperseg=nw,
                                         noverlap=no,nfft=nf,detrend=reftrend,
                                         scaling=scaling,average=average)
                    elif method == 'coherence':
                        ref = signal.coherence(x,y,fs=fs,window=swindow,nperseg=nw,
                                               noverlap=no,nfft=nf,detrend=reftrend)
                    elif method == 'stft':
                        ref = signal.stft(x,fs=fs,window=swindow,nperseg=nw,
                                          noverlap=no,nfft=nf,detrend=reftrend,
                                          boundary=None if boundary=='none' else 'zeros',
                                          padded=padded,scaling=scaling)
                    else:
                        ref = signal.spectrogram(x,fs=fs,window=swindow,nperseg=nw,
                                                 noverlap=no,nfft=nf,detrend=reftrend,
                                                 scaling=scaling)
                for item in range(len(result)):
                    actual, target = result[item], ref[item]
                    if method == 'coherence' and item==1:
                        # The ratio is inherently undefined where both spectral estimates
                        # are near machine zero after linear detrending.
                        if nw <= 2 and trend=='linear':
                            continue
                        # Tiny spectral bins after detrending are dominated by
                        # floating-point cancellation: coherence is undefined.
                        pxx = signal.welch(x, fs=fs, window=swindow, nperseg=nw,
                                           noverlap=no, nfft=nf, detrend=reftrend)[1]
                        pyy = signal.welch(y, fs=fs, window=swindow, nperseg=nw,
                                           noverlap=no, nfft=nf, detrend=reftrend)[1]
                        mask = np.isfinite(target) & (pxx > np.max(pxx)*1e-18) & \
                               (pyy > np.max(pyy)*1e-18)
                        actual, target = actual[mask], target[mask]
                    same(actual,target,f'{method} {case} field {item}')
            except Exception:
                print('Failed',method,'case',case,'n',n,'ny',ny,'nw',nw,
                      'no',no,'nf',nf,'window',window,'trend',trend)
                raise
            count += 1
        print(f'{method}: 160 matched')
    return count


def main():
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        cmd = ['gfortran','-std=f2018','-O2',f'-J{work}',f'-I{work}']
        cmd += [str(ROOT/'src'/source) for source in SOURCES]
        cmd += [str(ROOT/'tools/signal_reference/spectral_complex_driver.f90'),'-o',str(work/'driver')]
        subprocess.run(cmd,cwd=work,check=True)
        count = verify(work/'driver')
        print(f'{count} complex/two-sided signal cases matched SciPy')

if __name__ == '__main__':
    main()
