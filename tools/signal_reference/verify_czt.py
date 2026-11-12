#!/usr/bin/env python3
"""Cross-check Fortran CZT/zoom FFT against scipy.signal (seeded cases)."""
from __future__ import annotations
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from scipy import signal
import scipy

ROOT = Path(__file__).resolve().parents[2]

def main():
    rng = np.random.default_rng(20261017)
    cases = []
    for i in range(480):
        mode = i % 4 + 1
        n = int(rng.integers(1, 128 if mode != 2 else 65))
        m = int(rng.integers(1, 161 if mode != 2 else 91))
        x = rng.standard_normal(n) + 1j*rng.standard_normal(n)
        fs = float(rng.choice([1.,2.,10.,100.]))
        fn = np.sort(rng.uniform(-.4*fs, .75*fs,2))
        endpoint = bool(rng.integers(2))
        a = 1.+0.j
        w = 1.+0.j
        if mode == 1:
            expected = signal.czt(x,m=m)
        elif mode == 2:
            theta = rng.uniform(-np.pi,np.pi)
            phi = rng.uniform(-np.pi,np.pi)
            if i%8 == 1:
                # Direct Horner for a spiral just off the unit circle.
                n = min(n,20)
                m = min(m,25)
                x = x[:n]
                w = rng.uniform(.985,1.015)*np.exp(1j*theta)
                a = rng.uniform(.985,1.015)*np.exp(1j*phi)
            else:
                w = np.exp(1j*theta)
                a = np.exp(1j*phi)
            expected = signal.czt(x,m=m,w=w,a=a)
        elif mode == 3:
            expected = signal.zoom_fft(x,fn,m=m,fs=fs,endpoint=endpoint)
        else:
            expected = signal.zoom_fft(x,fn[1],m=m,fs=fs,endpoint=endpoint)
        cases.append((mode,x,m,fn,fs,int(endpoint),a,w,expected))

    source = [str(len(cases))]
    for mode,x,m,fn,fs,endpt,a,w,_ in cases:
        source.append(f'{len(x)} {m} {mode}')
        source.append(f'{fn[0]:.17g} {fn[1]:.17g} {fs:.17g} {endpt}')
        source.append(f'{a.real:.17g} {a.imag:.17g} {w.real:.17g} {w.imag:.17g}')
        source.append(' '.join(f'{z:.17g}' for z in x.real))
        source.append(' '.join(f'{z:.17g}' for z in x.imag))

    with tempfile.TemporaryDirectory(prefix='scifort-czt-') as td:
        root=Path(td)
        cmd=['gfortran','-std=f2018','-O2','-Werror=implicit-interface','-fcheck=all',
             '-I',str(root),'-J',str(root)]
        cmd += [str(ROOT/'src'/f'{s}.f90') for s in ('kinds','signal_dpss', 'signal_windows_extra', 'signal_spectral','signal_czt')]
        cmd += [str(ROOT/'tools/signal_reference/czt_driver.f90'),'-o',str(root/'driver')]
        subprocess.run(cmd,check=True,cwd=root)
        process = subprocess.run([str(root/'driver')],input='\n'.join(source)+'\n',
                                 text=True,capture_output=True,cwd=root,check=True,timeout=180)
    data = iter(process.stdout.splitlines())
    worst = 0.
    for i,(mode,x,m,fn,fs,endpt,a,w,expected) in enumerate(cases):
        st,n = map(int,next(data).split())
        if st!=0 or n!=m:
            raise AssertionError(f'case {i} mode {mode}: status {st}, length {n}')
        y = np.fromstring(next(data),sep=' ') + 1j*np.fromstring(next(data),sep=' ')
        norm = max(1.,float(np.max(np.abs(expected))))
        scaled = float(np.max(np.abs(y-expected)))/norm
        worst = max(worst,scaled)
        if scaled > 2.e-9:
            raise AssertionError(f'case {i} mode={mode}, n={len(x)}, m={m}, error={scaled:.3e}')
    if next(data,None) is not None:
        raise AssertionError('Unexpected trailing output')
    print(f'PASS: {len(cases)}/{len(cases)} czt/zoom_fft comparisons against SciPy {scipy.__version__}; '
          f'maximum scaled error {worst:.3e}')

if __name__ == '__main__':
    main()
