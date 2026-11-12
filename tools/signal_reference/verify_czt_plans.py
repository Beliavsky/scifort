#!/usr/bin/env python3
"""Compare reusable Fortran CZT/ZoomFFT plans and batched axes with SciPy."""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import numpy as np
import scipy
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]


def main():
    rng = np.random.default_rng(20261024)
    cases = []
    for i in range(520):
        mode = i % 5 + 1
        n = int(rng.integers(1, 115 if mode != 3 else 25))
        m = int(rng.integers(1, 151 if mode != 3 else 35))
        channels = int(rng.integers(1, 6))
        axis = i % 2 + 1  # Fortran one-based dimension.
        real_input = bool((i // 5) % 2)
        endpoint = bool((i // 13) % 2)
        if mode >= 4 and m == 1:  # SciPy 1.17 ZoomFFT divides by zero for endpoint=True, m=1.
            endpoint = False
        fs = float(rng.choice([1., 2., 4., 10., 48.]))
        fn = np.sort(rng.uniform(-.4*fs, .8*fs, size=2))
        a = complex(1.)
        w = complex(1.)
        if mode == 2:
            a = np.exp(1j*rng.uniform(-2, 2))
            w = np.exp(1j*rng.uniform(-2, 2))
        if mode == 3:
            a = rng.uniform(.99, 1.01)*np.exp(1j*rng.uniform(-2, 2))
            w = rng.uniform(.99, 1.01)*np.exp(1j*rng.uniform(-2, 2))
        shape = (n,channels) if axis == 1 else (channels,n)
        x = rng.standard_normal(shape)
        if not real_input:
            x = x + 1j*rng.standard_normal(shape)
        if mode == 1:
            plan = signal.CZT(n,m=m)
        elif mode in (2,3):
            plan = signal.CZT(n,m=m,w=w,a=a)
        elif mode == 4:
            plan = signal.ZoomFFT(n,fn,m=m,fs=fs,endpoint=endpoint)
        else:
            plan = signal.ZoomFFT(n,float(fn[1]),m=m,fs=fs,endpoint=endpoint)
        expected = plan(x,axis=axis-1)
        if not np.all(np.isfinite(expected)):
            raise AssertionError(f'SciPy produced nonfinite reference values in case {i}')
        cases.append((mode,n,m,channels,axis,real_input,endpoint,fn,fs,w,a,x,expected))

    lines = [str(len(cases))]
    for mode,n,m,channels,axis,real_input,endpoint,fn,fs,w,a,x,_ in cases:
        lines.append(f'{mode} {n} {m} {channels} {axis} {int(real_input)} {int(endpoint)}')
        lines.append(' '.join(f'{r:.17g}' for r in (*fn,fs,w.real,w.imag,a.real,a.imag)))
        flat = x.ravel(order='F')
        lines.append(' '.join(f'{z:.17g}' for z in flat.real))
        lines.append(' '.join(f'{z:.17g}' for z in flat.imag))
    with tempfile.TemporaryDirectory(prefix='scifort-czt-plan-') as td:
        build = Path(td)
        flags = ['gfortran','-std=f2018','-O2','-fcheck=all','-Werror=implicit-interface',
                 '-J',str(build),'-I',str(build)]
        sources = ('kinds','signal_dpss','signal_windows_extra','signal_spectral','signal_czt')
        cmd = flags+[str(ROOT/'src'/f'{s}.f90') for s in sources]
        cmd += [str(ROOT/'tools/signal_reference/czt_plan_driver.f90'),'-o',str(build/'driver')]
        subprocess.run(cmd,cwd=build,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        process = subprocess.run([str(build/'driver')],input='\n'.join(lines)+'\n',
                                 cwd=build,text=True,capture_output=True,check=True,timeout=180)
    data = iter(process.stdout.splitlines())
    worst = 0.
    for i,case in enumerate(cases):
        mode,n,m,channels,axis,real_input,endpoint,fn,fs,w,a,x,expected = case
        st,nrow,ncol = map(int,next(data).split())
        if st != 0 or (nrow,ncol)!=expected.shape:
            raise AssertionError(f'case {i}, mode {mode}: status={st}, shape={(nrow,ncol)}, wanted={expected.shape}')
        real = np.fromstring(next(data),sep=' ')
        imag = np.fromstring(next(data),sep=' ')
        got = (real+1j*imag).reshape((nrow,ncol),order='F')
        if not np.all(np.isfinite(got)):
            raise AssertionError(f'Fortran produced nonfinite values in case {i}')
        scaled = float(np.max(np.abs(got-expected)))/max(1.,float(np.max(np.abs(expected))))
        worst = max(worst,scaled)
        if scaled > 2e-9:
            raise AssertionError(f'case {i}, mode={mode}, shape={x.shape}, axis={axis}: {scaled:.3e}')
    if next(data,None) is not None:
        raise AssertionError('Unexpected output from Fortran')
    print(f'PASS: {len(cases)}/{len(cases)} reusable CZT/ZoomFFT cases, '
          f'{sum(c[-1].size for c in cases):,} complex values, '
          f'max scaled difference {worst:.3e} vs SciPy {scipy.__version__}')


if __name__ == '__main__':
    main()
