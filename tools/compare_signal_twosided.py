#!/usr/bin/env python3
"""Independently validate two-sided legacy SciFort spectral estimators against SciPy.

Run using: python tools/compare_signal_twosided.py /path/to/reference_signal_twosided
The Fortran fixture is deterministic; no SciPy output is used in its generation.
"""
from __future__ import annotations
import argparse
import subprocess
import warnings
import numpy as np
from scipy import signal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    args = parser.parse_args()
    output = subprocess.check_output([args.driver], text=True)
    lines = iter(output.splitlines())
    n = 37
    k = np.arange(1, n+1, dtype=np.float64)
    x = .18*k+.7*np.sin(.48*k)+.36*np.cos(.17*k*k)
    y = -.06*k+.8*np.cos(.26*k)+.25*np.sin(.37*k)
    z = x+1j*y
    worst = 0.0
    nvalues = 0
    ncases = 0
    for row in lines:
        name, *arguments = row.split()
        assert name == 'CASE', row
        family,nfft,wi,si,oi,nfreq,ntime = map(int, arguments)
        onesided = oi == 1
        window = 'hann' if wi == 1 else ('gaussian', 2.25)
        scale = 'density' if si == 1 else 'spectrum'
        trend = 'constant' if oi == 1 else 'linear'
        average = 'mean' if wi == 1 else 'median'
        f_obs = np.array([float(next(lines)) for _ in range(nfreq)])
        common = dict(fs=25.,window=window,nfft=nfft,detrend=trend,
                      return_onesided=onesided)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if family == 1:
                f, result = signal.periodogram(x,scaling=scale,**common)
            elif family == 2:
                f, result = signal.welch(x,nperseg=8,noverlap=3,scaling=scale,
                                         average=average,**common)
            elif family == 3:
                f, result = signal.csd(x,y,nperseg=8,noverlap=3,scaling=scale,
                                       average=average,**common)
            elif family == 4:
                f, t, result = signal.stft(x,nperseg=8,noverlap=3,scaling=('psd' if si==1 else 'spectrum'),
                                           boundary=('even' if wi==1 else 'odd'),padded=(wi==2),**common)
            elif family in (5,6):
                f, t, result = signal.spectrogram(x if family==5 else z,
                    nperseg=8,noverlap=3,scaling=scale,mode='psd',**common)
            else:
                raise AssertionError(family)
        if family in (1,2,3):
            if family == 3:
                values = np.array([complex(*map(float,next(lines).split())) for _ in range(nfreq)])
            else:
                values = np.array([float(next(lines)) for _ in range(nfreq)])
        else:
            t_obs = np.array([float(next(lines)) for _ in range(ntime)])
            if family == 4:
                values = np.array([complex(*map(float,next(lines).split()))
                                   for _ in range(nfreq*ntime)]).reshape(ntime,nfreq).T
            else:
                values = np.array([float(next(lines)) for _ in range(nfreq*ntime)]).reshape(ntime,nfreq).T
            assert np.array_equal(t_obs.shape,t.shape),(family,nfft,wi,si,oi)
            assert np.allclose(t_obs,t,rtol=1.e-11,atol=1.e-11), (family,nfft,wi,si,oi,'times',np.max(abs(t_obs-t)))
        identifier=(family,nfft,wi,si,oi)
        assert f_obs.shape==f.shape, (identifier,'frequency shape',f_obs.shape,f.shape)
        assert np.allclose(f_obs,f,rtol=1.e-11,atol=1.e-11),(identifier,'frequencies',f_obs,f)
        assert values.shape==result.shape, (identifier,'value shape',values.shape,result.shape)
        deviation=np.max(np.abs(values-result)) if values.size else 0.0
        assert np.allclose(values,result,rtol=3.e-11,atol=5.e-12),(identifier,'value mismatch',deviation,
            np.unravel_index(np.argmax(abs(values-result)),values.shape))
        worst=max(worst,deviation)
        nvalues += values.size
        ncases += 1
    assert ncases == 192,ncases
    print(f'PASS: {ncases} SciPy 1.17 two-sided spectral comparisons, '
          f'{nvalues} coefficients, max absolute difference {worst:.3e}')

if __name__=='__main__':
    main()
