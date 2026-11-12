#!/usr/bin/env python3
"""Compare mixed-dtype SciFort cross-spectrum outputs to SciPy 1.17.

For legacy CSD/coherence, both SciPy operands are explicitly promoted to
complex128. SciPy 1.17's direct mixed-dtype dispatch is order-dependent.
For ShortTimeFFT.spectrogram, SciPy is called with the original mixed types.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy import signal


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    args = parser.parse_args()
    data = iter(args.results.read_text().splitlines())

    index = np.arange(1.0, 80.0)
    x = 0.7*np.cos(0.19*index) + np.sin(0.43*index) + 0.013*index
    z = (0.4*np.sin(0.23*index) + x + 1j*(np.cos(0.37*index) - 0.15*np.sin(0.09*index)))
    ncase = 0
    ncoeff = 0
    largest = 0.0
    mode_counts: dict[str, int] = {'C': 0, 'H': 0, 'S': 0}

    for line in data:
        fields = line.split()
        tag = fields[0]
        if tag == 'C':
            _, pair, wid, nw, no, nf, avg_id, trend_id, scale_id, nout, status = fields
            pair, wid, nw, no, nf, avg_id, trend_id, scale_id, nout, status = map(
                int, (pair,wid,nw,no,nf,avg_id,trend_id,scale_id,nout,status)
            )
            values = np.array([list(map(float, next(data).split())) for _ in range(nout)])
            assert status == 0 and values.shape == (nf, 3)
            w = 'hann' if wid == 1 else ('gaussian', 3.5)
            a, b = ((x.astype(complex), z) if pair == 1 else (z, x.astype(complex)))
            frequencies, reference = signal.csd(
                a, b, fs=20.0, window=w, nperseg=nw, noverlap=no, nfft=nf,
                detrend='linear' if trend_id == 2 else False,
                scaling='spectrum' if scale_id == 2 else 'density',
                average='median' if avg_id == 2 else 'mean', return_onesided=False,
            )
            actual = values[:,1] + 1j*values[:,2]
            err = np.max(np.abs(reference - actual))
            grid = np.max(np.abs(frequencies - values[:,0]))
        elif tag == 'H':
            _, pair, wid, nw, no, nf, trend_id, nout, status = fields
            pair, wid, nw, no, nf, trend_id, nout, status = map(
                int, (pair,wid,nw,no,nf,trend_id,nout,status)
            )
            values = np.array([list(map(float, next(data).split())) for _ in range(nout)])
            assert status == 0 and values.shape == (nf, 2)
            w = 'hann' if wid == 1 else ('gaussian', 3.5)
            a,b = ((x.astype(complex),z) if pair == 1 else (z,x.astype(complex)))
            frequencies, reference = signal.coherence(
                a,b,fs=20.0,window=w,nperseg=nw,noverlap=no,nfft=nf,
                detrend='linear' if trend_id == 2 else False,
            )
            actual = values[:,1]
            err = np.max(np.abs(reference - actual))
            grid = np.max(np.abs(frequencies-values[:,0]))
        elif tag == 'S':
            _, pair, mode_id, pad_id, phase, scale_id, nw, hop, nf, n, nfreq, ntime = fields
            pair, mode_id, pad_id, phase, scale_id, nw, hop, nf, n, nfreq, ntime = map(
                int, (pair,mode_id,pad_id,phase,scale_id,nw,hop,nf,n,nfreq,ntime)
            )
            values = np.array([list(map(float, next(data).split())) for _ in range(nfreq*ntime)])
            assert values.shape == (nfreq*ntime, 2)
            actual = (values[:,0]+1j*values[:,1]).reshape((ntime,nfreq)).T
            mode = 'twosided' if mode_id == 1 else 'centered'
            pad = {1:'zeros',2:'edge',3:'even',4:'odd'}[pad_id]
            s = signal.ShortTimeFFT(
                signal.get_window('hann',nw),hop,20.0,mfft=nf,fft_mode=mode,
                scale_to='psd' if scale_id == 1 else 'magnitude',phase_shift=phase,
            )
            a,b = ((x,z) if pair == 1 else (z,x))
            reference = s.spectrogram(a[:n],b[:n],padding=pad)
            if reference.shape != actual.shape:
                raise AssertionError(f'ShortTimeFFT shape mismatch {reference.shape} vs {actual.shape}')
            err = np.max(np.abs(reference-actual))
            grid = 0.0
        else:
            raise AssertionError(f'Unexpected record: {line[:100]!r}')
        if not np.isfinite(err) or err > 2e-10*max(1.0,float(np.max(np.abs(reference)))):
            raise AssertionError(f'Case {ncase+1} ({tag}) error {err}; {fields}')
        if grid > 2e-12:
            raise AssertionError(f'Case {ncase+1} ({tag}) frequency error {grid}; {fields}')
        largest = max(largest,float(err))
        ncoeff += int(actual.size)
        ncase += 1
        mode_counts[tag] += 1
    if not (ncase >= 250 and all(v > 0 for v in mode_counts.values())):
        raise AssertionError(f'Insufficient comparison cases: {mode_counts}')
    print(f'PASS: {ncase} mixed spectral cases ({mode_counts}), '
          f'{ncoeff} coefficients, max error {largest:.3e}')


if __name__ == '__main__':
    main()
