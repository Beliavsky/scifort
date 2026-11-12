#!/usr/bin/env python3
"""Validate SciFort FFT and overlap-add convolution against SciPy."""
from __future__ import annotations

import pathlib
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = pathlib.Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261015)
SOURCES = (
    'kinds.f90', 'signal_dpss.f90', 'signal_windows_extra.f90', 'signal_spectral.f90', 'signal_fir_advanced.f90', 'signal_spectral_complex.f90',
    'signal_convolve.f90', 'signal_convolve2d.f90', 'signal_fftconvolve.f90', 'signal_detrend.f90',
    'signal_peaks.f90', 'signal_savgol.f90', 'signal_lfilter.f90', 'signal_sos.f90', 'signal_filtfilt.f90', 'signal_filtfilt_nd.f90', 'signal_analog_design.f90', 'signal_design.f90', 'signal_elliptic.f90', 'signal_analog_sos.f90', 'signal_iir_extra.f90', 'signal_fourier_resample.f90', 'signal_resample.f90', 'signal_analog.f90', 'signal_zpk_conversion.f90', 'signal_state_space.f90', 'signal_cont2discrete.f90', 'signal_discrete_simulation.f90', 'signal_continuous_simulation.f90', 'signal_lti_frequency.f90', 'signal_partial_fraction.f90', 'signal_sequences.f90', 'signal_hilbert.f90', 'signal_waveforms.f90', 'signal_periodic.f90', 'signal_irregular.f90', 'signal_remez.f90', 'signal_czt.f90', 'signal.f90',
)


def run(exe, method, mode, x, y):
    complex_mode = np.iscomplexobj(x) or np.iscomplexobj(y)
    lines = [method,mode,f'{len(x)} {len(y)} {int(complex_mode)}']
    lines += [f'{z.real:.17g} {z.imag:.17g}' for z in x]
    lines += [f'{z.real:.17g} {z.imag:.17g}' for z in y]
    result = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',
                            capture_output=True,text=True,timeout=20)
    if result.returncode:
        raise RuntimeError(result.stderr)
    output = result.stdout.strip().splitlines()
    status,n = map(int,output[0].split())
    if status != 0:
        raise AssertionError(f'{method} returned status {status}')
    arr = np.array([[float(q) for q in row.split()] for row in output[1:]])
    if len(arr) != n:
        raise AssertionError(f'{method}: count {len(arr)} vs {n}')
    return arr[:,0]+1j*arr[:,1]


def verify(exe):
    count = 0
    for method in ('fftconvolve','oaconvolve'):
        for complex_mode in (False, True):
            for case in range(250):
                mode = ('full','same','valid')[case%3]
                n = int(RNG.integers(1,270))
                m = int(RNG.integers(1,34))
                if case%5==0:
                    n,m = m,n
                if case%7==0:
                    n = int(RNG.integers(400,1050))
                    m = int(RNG.integers(1,20))
                x = RNG.standard_normal(n)
                y = RNG.standard_normal(m)
                if complex_mode:
                    x = x+1j*RNG.standard_normal(n)
                    y = y+1j*RNG.standard_normal(m)
                got = run(exe,method,mode,x,y)
                reference = getattr(signal,method)(x,y,mode=mode)
                if got.shape != reference.shape or not np.allclose(got,reference,atol=2e-9,rtol=4e-9):
                    err = np.max(np.abs(got-reference)) if got.shape==reference.shape else None
                    raise AssertionError(f'{method} {complex_mode} {case}, n={n} m={m} '
                                         f'mode={mode} err={err}')
                count += 1
            print(f'{method}, complex={complex_mode}: 250 cases matched')
    return count


def main():
    with tempfile.TemporaryDirectory() as tmp:
        work=pathlib.Path(tmp)
        cmd=['gfortran','-std=f2018','-O2',f'-J{work}',f'-I{work}']
        cmd += [str(ROOT/'src'/src) for src in SOURCES]
        cmd += [str(ROOT/'tools/signal_reference/fftconvolve_driver.f90'),
                '-o',str(work/'driver')]
        subprocess.run(cmd,cwd=work,check=True)
        count=verify(work/'driver')
        print(f'{count} FFT convolution cases matched SciPy')

if __name__=='__main__':
    main()
