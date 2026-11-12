#!/usr/bin/env python3
"""Numerical cross-check of parameterized Fourier resample against SciPy 1.17."""
from __future__ import annotations
import subprocess
import sys
import numpy as np
from scipy import signal

WINDOWS = {
    1: ('kaiser', 4.2, 0.0, ('kaiser',4.2)),
    2: ('gaussian', 3.1, 0.0, ('gaussian',3.1)),
    3: ('tukey', .37, 0.0, ('tukey',.37)),
    4: ('general_gaussian', 1.9, 3.7, ('general_gaussian',1.9,3.7)),
    5: ('taylor', 3., 34., ('taylor',3,34)),
    6: ('general_cosine', 0., 0., ('general_cosine',[.6,.3,.1])),
    7: ('dpss', 2.2, 0., ('dpss',2.2)),
}


def main(executable: str) -> None:
    rng = np.random.default_rng(20261028)
    cases = []
    for family, (name, p, p2, window) in WINDOWS.items():
        for mode in (0,1,2):
            for nx in (7, 12, 19, 28):
                for ny in (5, 14, 25):
                    x = rng.normal(size=nx) + .16*np.sin(np.arange(nx)*.31)
                    if mode != 0:
                        x = x + 1j*rng.normal(size=nx)
                    if mode == 2:
                        x = np.fft.fft(x)
                    cases.append((mode,family,nx,ny,p,p2,x,window))
    lines = [str(len(cases))]
    for mode,family,nx,ny,p,p2,x,window in cases:
        lines.append(f'{mode} {family} {nx} {ny} {p:.17e} {p2:.17e}')
        lines.append(' '.join(f'{v:.17e}' for v in np.real(x)))
        lines.append(' '.join(f'{v:.17e}' for v in np.imag(x)))
    proc = subprocess.run([executable], input='\n'.join(lines)+'\n', text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          timeout=180,check=True)
    it = iter(proc.stdout.splitlines())
    max_error, count = 0., 0
    for index, (mode,family,nx,ny,p,p2,x,window) in enumerate(cases):
        reference_window = window
        if family == 6:
            # SciPy 1.17 refuses general_cosine tuples in resample's array-API
            # path. Independently make its periodic window and shift to FFT order.
            reference_window = np.fft.fftshift(
                signal.windows.general_cosine(nx, [.6,.3,.1], sym=False))
        ref=signal.resample(x, ny, window=reference_window,
                            domain='freq' if mode==2 else 'time')
        st,num=map(int,next(it).split())
        if st!=0 or num!=ny:
            raise AssertionError(f'case {index} status={st}, n={num}, want={ny}')
        got=np.array([complex(*map(float,next(it).split())) for _ in range(num)])
        err=float(np.max(np.abs(got-ref)))
        max_error=max(max_error,err)
        count += num
        if err > 3.e-10:
            k=int(np.argmax(np.abs(got-ref)))
            raise AssertionError(f'case {index}, window={window} mode={mode}, n={nx}->{ny}: '
                                 f'bad sample {k}: {got[k]} vs {ref[k]}, error {err}')
    print(f'PASS: {len(cases)} SciPy Fourier resampling comparisons / {count} output samples; '
          f'maximum absolute error {max_error:.3e}',flush=True)

if __name__=='__main__':
    main(sys.argv[1])
