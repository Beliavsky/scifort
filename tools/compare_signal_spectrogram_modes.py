#!/usr/bin/env python3
"""Compare legacy spectrogram PSD/complex/magnitude/angle/phase with SciPy 1.17.

Invoked with the path to the compiled independent Fortran driver. The reference
is computed from SciPy on seeded inputs, not from Fortran-generated expectations.
"""
from __future__ import annotations

import sys
import subprocess
import warnings
import numpy as np
from scipy import signal
from scipy.signal import windows


def make_cases():
    rng = np.random.default_rng(20261028)
    cases = []
    for kind in (0, 1):
        for mode in range(1, 6):
            for win in range(4):
                for choice in range(8):
                    nw = (8, 9, 16, 17, 23, 24, 30, 31)[choice]
                    n = nw * (2 + choice % 3) + choice + 9
                    no = [0, 2, nw // 4, nw-2][(choice+mode) % 4]
                    nf = nw + [0, 1, 3][(choice+mode) % 3]
                    onesided = (choice + mode) % 2
                    trend = 1 + (choice+mode) % 3
                    scale = 1 + (choice//2) % 2
                    fs = [2.0, 125.0, 360.0][(choice + win + kind) % 3]
                    x = rng.normal(size=n) + 0.2*np.cos(np.arange(n)*0.19)
                    if kind:
                        x = x + 1j*(rng.normal(size=n) + 0.1*np.sin(np.arange(n)*0.49))
                    cases.append((kind,mode,win,n,nw,no,nf,onesided,trend,scale,fs,x))
    # Exercise omitted nperseg, noverlap and nfft defaults independently.
    for kind in (0, 1):
        for mode in range(1, 6):
            for win in range(3):
                for choice in range(4):
                    n = [22, 33, 57, 71][choice]
                    nw = [0, 15, 0, 13][choice]
                    no = [-1, -1, 1, -1][choice]
                    nf = [0, 0, 0, 17][choice]
                    onesided = (choice+mode+win) % 2
                    trend = 1 + (choice+mode) % 3
                    scale = 1 + choice % 2
                    fs = [2.0, 100.0][(choice+win) % 2]
                    x = rng.normal(size=n) + 0.4*np.cos(np.arange(n)*0.41)
                    if kind:
                        x = x + 1j*(rng.normal(size=n)+0.3*np.sin(np.arange(n)*0.2))
                    cases.append((kind,mode,win,n,nw,no,nf,onesided,trend,scale,fs,x))
    return cases


def compare(executable: str):
    cases = make_cases()
    payload = [str(len(cases))]
    for kind,mode,win,n,nw,no,nf,onesided,trend,scale,fs,x in cases:
        payload.append(f'{kind} {mode} {win} {n} {nw} {no} {nf} {onesided} {trend} {scale} {fs:.17e}')
        payload.append(' '.join(f'{v:.17e}' for v in x.real))
        payload.append(' '.join(f'{v:.17e}' for v in x.imag))
        if win == 3:
            payload.append(' '.join(f'{v:.17e}' for v in windows.blackmanharris(nw, sym=False)))
    result = subprocess.run([executable], input='\n'.join(payload)+'\n', text=True,
                            capture_output=True, check=True)
    lines = iter(result.stdout.splitlines())
    coefficient_count = 0
    maximum = 0.0
    for idx, (kind,mode,win,n,nw,no,nf,onesided,trend,scale,fs,x) in enumerate(cases):
        mname = ('psd', 'complex', 'magnitude', 'angle', 'phase')[mode-1]
        wt = (('tukey', .25), 'hann', ('gaussian', 4.3), windows.blackmanharris(nw,sym=False))[win]
        kwargs = dict(window=wt, fs=fs,
                      scaling=('density','spectrum')[scale-1],
                      detrend=('constant','linear',False)[trend-1], mode=mname)
        if nw > 0:
            kwargs['nperseg'] = nw
        if no >= 0:
            kwargs['noverlap'] = no
        if nf > 0:
            kwargs['nfft'] = nf
        kwargs['return_onesided'] = bool(onesided) if kind == 0 else False
        with warnings.catch_warnings():
            # Deliberately omit nperseg in some cases to verify SciPy's
            # documented shortening of the 256-sample default window.
            warnings.filterwarnings('ignore', message=r'nperseg = 256 is greater.*')
            freq, time, truth = signal.spectrogram(x,**kwargs)
        s, nfreq, ntime = map(int,next(lines).split())
        if s != 0 or nfreq != len(freq) or ntime != len(time):
            raise AssertionError(f'case {idx} status/shape {(s,nfreq,ntime)} != {(len(freq),len(time))}')
        ff = np.array([float(next(lines)) for _ in range(nfreq)])
        tt = np.array([float(next(lines)) for _ in range(ntime)])
        out = np.array([complex(*map(float,next(lines).split())) for _ in range(nfreq*ntime)])
        out = out.reshape(ntime,nfreq).T
        if not np.allclose(ff,freq,rtol=1e-12,atol=1e-11):
            raise AssertionError(f'case {idx} frequency error {max(abs(ff-freq))}')
        if not np.allclose(tt,time,rtol=1e-12,atol=1e-11):
            raise AssertionError(f'case {idx} time error {max(abs(tt-time))}')
        differences = abs(out-truth)
        delta = float(np.max(differences))
        maximum = max(delta,maximum)
        if delta > 5e-10:
            mxidx = np.unravel_index(differences.argmax(),truth.shape)
            raise AssertionError(f'case {idx}: {kind=} {mname=} {win=} {scale=} {nf=} '
                                 f'delta={delta:.4g} at {mxidx}, expected={truth[mxidx]}, '
                                 f'observed={out[mxidx]}')
        if mname == 'phase' and np.any(abs(np.diff(out.real,axis=0)) > np.pi+1e-9):
            raise AssertionError(f'case {idx}: phase was not unwrapped')
        coefficient_count += truth.size
    print(f'PASS: {len(cases)} SciPy spectrogram cases; {coefficient_count} coefficients; '
          f'maximum absolute error {maximum:.3e}')


if __name__ == '__main__':
    compare(sys.argv[1])
