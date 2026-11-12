#!/usr/bin/env python3
"""Compare SciFort inverse STFT for tuple/window-vector families to SciPy 1.17.

Called by tools/run_signal_istft_parameter_windows.py with one argument: the
path to the previously compiled Fortran reference program.  All STFT input
coefficients originate in SciPy; the Fortran program never generates them.
"""
from __future__ import annotations

import io
from pathlib import Path
import subprocess
import sys
import warnings

import numpy as np
from scipy import signal


def test_cases():
    rng = np.random.default_rng(884011)
    windows = [
        ('hann', 0, 0, 0, None),
        ('gaussian', 1, 2.4, 0, None),
        ('kaiser', 1, 4.5, 0, None),
        ('tukey', 1, 0.35, 0, None),
        ('exponential', 1, 2.1, 0, None),
        ('chebwin', 1, 65.0, 0, None),
        ('dpss', 1, 2.0, 0, None),
        ('general_hamming', 1, 0.63, 0, None),
        ('general_gaussian', 2, 1.6, 2.7, None),
        ('taylor', 2, 4.0, 37.0, None),
        ('general_cosine', 3, 0, 0, [0.48, 0.36, 0.16]),
        ('explicit', 4, 0, 0, None),
    ]
    cases=[]
    for j in range(192):
        name, mode, a, b, cs = windows[j % len(windows)]
        nw = [8, 9, 10, 11, 12, 13][(j // len(windows)) % 6]
        nfft = nw + [0, 1, 5, 6][(j // 3) % 4]
        no = nw - [2, 3, 4][(j // 2) % 3]
        fs = [1.0, 10.0, 48.0][(j // 9) % 3]
        real = j % 2 == 0
        scaling = ['spectrum', 'psd'][(j // 5) % 2]
        size = nw+13+(j % 11)
        x = rng.normal(size=size)
        if not real:
            x=x+1j*0.7*rng.normal(size=size)
        if mode==0:
            w=name
        elif mode==1:
            if name=='exponential':
                w=(name,None,a)
            else:
                w=(name,a)
        elif mode==2:
            w=(name,int(a) if name=="taylor" else a,b)
        elif mode==3:
            w=(name,np.asarray(cs))
        else:
            w=np.exp(-0.4*(np.arange(nw)-(nw-1)/2)**2)+0.08
            cs=w
        # The legacy API pads trailing data, and ISTFT retains that padding.
        with warnings.catch_warnings():
            warnings.simplefilter('error', category=RuntimeWarning)
            f,t,z=signal.stft(x,fs=fs,window=w,nperseg=nw,noverlap=no,nfft=nfft,
                              detrend=False,boundary=['zeros','even','odd','constant'][j%4],
                              padded=True,scaling=scaling,return_onesided=real)
            tt,yy=signal.istft(z,fs=fs,window=w,nperseg=nw,noverlap=no,nfft=nfft,
                                input_onesided=real,boundary=True,scaling=scaling)
        cases.append((j+1,name,mode,a,b,np.asarray(cs if cs is not None else [],dtype=float),
                      real,nw,no,nfft,fs,scaling,z,tt,yy))
    return cases


def compare(executable:Path):
    cases=test_cases()
    stdin=io.StringIO()
    stdin.write(f'{len(cases)}\n')
    for _,name,mode,a,b,co,real,nw,no,nfft,fs,scaling,z,tt,yy in cases:
        nf,nt=z.shape
        stdin.write(f'{nf} {nt} {nw} {no} {nfft} {1 if real else 2} {mode} {fs:.17e} {a:.17e} {b:.17e} {len(co)}\n')
        stdin.write(f'{name if name != "explicit" else "boxcar"} {scaling}\n')
        for val in co:
            stdin.write(f'{val:.17e}\n')
        for val in z.ravel(order='F'):
            stdin.write(f'{val.real:.17e} {val.imag:.17e}\n')
    result=subprocess.run([str(executable)],input=stdin.getvalue(),capture_output=True,text=True)
    if result.returncode:
        raise RuntimeError(f'Fortran failed: {result.stderr}\nstdout: {result.stdout[-1000:]}')
    lines=iter(result.stdout.splitlines())
    worst_error=0.0
    nvalues=0
    for case in cases:
        ix,name,mode,a,b,co,real,nw,no,nfft,fs,scaling,z,tt,yy=case
        got_case,code,n=map(int,next(lines).split())
        if (got_case,code,n)!=(ix,0,len(yy)):
            raise AssertionError((ix,name,'status/size',got_case,code,n,len(yy)))
        data=np.asarray([[float(v) for v in next(lines).split()] for k in range(n)])
        got_t=data[:,0]
        got_y=data[:,1]+1j*data[:,2]
        delta=np.max(np.abs(got_y-yy)) if n else 0.
        time_error=np.max(np.abs(got_t-tt)) if n else 0.
        worst_error=max(delta,worst_error)
        if not np.allclose(got_y,yy,rtol=3e-10,atol=2e-10):
            raise AssertionError((ix,name,mode,real,scaling,nw,nfft,'max error',delta))
        if not np.allclose(got_t,tt,rtol=4e-13,atol=4e-13):
            raise AssertionError((ix,name,mode,'time grid',time_error))
        nvalues+=n
    print(f'PASS: {len(cases)} independent SciPy 1.17 inverse STFT cases; '
          f'{nvalues:,} reconstructed values; max absolute error {worst_error:.3g}')


if __name__=='__main__':
    compare(Path(sys.argv[1]).absolute())
