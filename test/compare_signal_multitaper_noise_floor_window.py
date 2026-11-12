"""Independent SciPy comparisons for windowed multitaper noise-floor estimation.

The SciFort API is an extension rather than a SciPy API. The reference uses
SciPy DPSS windows, signal.detrend, ndimage.median_filter, and NumPy FFT.

Usage: python test/compare_signal_multitaper_noise_floor_window.py [gfortran] [O0|O3]
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import detrend as scipy_detrend
from scipy.signal.windows import dpss

from compare_signal_multitaper_time_frequency import multitaper_weights

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / (name + '.f90') for name in
           ('kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral', 'signal_multitaper')]
DRIVER = ROOT / 'reference' / 'signal_multitaper_noise_floor_window_reference_driver.f90'


def reference(data, nw, seg, overlap, k, nf, fs, trend, weighting, width):
    real_input = not np.iscomplexobj(data)
    n, channels = data.shape
    tapers, eig = dpss(seg, nw, Kmax=k, sym=False, return_ratios=True)
    freq = np.fft.rfftfreq(nf, 1/fs) if real_input else np.fft.fftfreq(nf, 1/fs)
    starts = range(0, n-seg+1, seg-overlap)
    times = (np.asarray(list(starts)) + seg//2) / fs
    power = np.zeros((len(freq), len(times), channels))
    for frame, start in enumerate(starts):
        for ch in range(channels):
            x = data[start:start+seg,ch].copy()
            if trend == 1:
                x -= x.mean()
            elif trend == 2:
                x = scipy_detrend(x, type='linear')
            y = x[None,:]*tapers
            ft = np.fft.rfft(y, n=nf, axis=1) if real_input else np.fft.fft(y, n=nf, axis=1)
            per = abs(ft)**2/(fs*np.sum(tapers*tapers,axis=1)[:,None])
            q = multitaper_weights(per, eig, np.mean(abs(x)**2)/fs, weighting)
            power[:,frame,ch] = np.sum(q*per, axis=0)
    if real_input:
        if nf % 2 == 0:
            power[1:-1] *= 2
        else:
            power[1:] *= 2
    bg = median_filter(power, size=(width,1,1), mode=('reflect' if real_input else 'wrap'))
    normalized = np.divide(power,bg,out=np.zeros_like(power),where=bg>0)
    normalized[(power>0)&(bg==0)] = np.finfo(float).max
    return freq, times, power, bg, normalized


def main():
    fc = sys.argv[1] if len(sys.argv)>1 else os.environ.get('FC', 'gfortran')
    mode = sys.argv[2] if len(sys.argv)>2 else 'O0'
    if mode not in ('O0','O3'):
        raise ValueError('mode must be O0 or O3')
    rng = np.random.default_rng(20261026)
    cases, lines = [], []
    for index in range(120):
        n = int(rng.integers(30,94))
        seg = int(rng.integers(12,min(41,n)+1))
        overlap = int(rng.integers(0,seg))
        k = int(rng.integers(2,min(seg,7)+1))
        nw = float(rng.uniform(1.2,min(seg/2-0.1,4.5)))
        nf = seg+int(rng.integers(0,22))
        fs = float(rng.uniform(2,90))
        nc = int(rng.integers(1,4))
        real_input = index%2 == 0
        trend = index%3
        weighting = ('equal','eigen','adaptive')[(index//3)%3]
        wf = ('equal','eigen','adaptive').index(weighting)
        nfout = nf//2+1 if real_input else nf
        width = int(rng.choice(np.arange(3,nfout+1,2)))
        axis = 1+index%2
        data = rng.normal(size=(n,nc))
        if not real_input:
            data = data+1j*rng.normal(size=(n,nc))
        if index%21 == 0:
            data *= 0
        elif index%19 == 0:
            data = data+2.0*np.sin(np.arange(n)[:,None]*0.37)
        cases.append((data,nw,seg,overlap,k,nf,fs,trend,weighting,width))
        lines.append(f'{n} {nc} {seg} {overlap} {k} {nf} {int(not real_input)} '
                     f'{trend} {wf} {width} {axis} {nw:.17g} {fs:.17g}')
        lines.extend(f'{v.real:.17g} {v.imag:.17g}' for ch in range(nc) for v in data[:,ch])
    input_text = str(len(cases))+'\n'+'\n'.join(lines)+'\n'
    with tempfile.TemporaryDirectory(prefix='scifort-noise-window-') as directory:
        build = Path(directory)
        objs = []
        flags = ['-O3'] if mode == 'O3' else ['-O0', '-g', '-fcheck=all']
        for src in SOURCES:
            obj = build/(src.stem+'.o')
            subprocess.run([fc,'-std=f2018','-Werror=implicit-interface',*flags,
                            '-I',str(build),'-J',str(build),'-c',str(src),'-o',str(obj)],
                           cwd=build,check=True)
            objs.append(str(obj))
        exe=build/'noise_window_reference'
        subprocess.run([fc,'-std=f2018',*flags,'-I',str(build),str(DRIVER),
                        *objs,'-o',str(exe)],cwd=build,check=True)
        output=subprocess.run([str(exe)],input=input_text,text=True,capture_output=True,
                              cwd=build,check=True).stdout.splitlines()
    lines_iter=iter(output)
    max_power=max_background=max_ratio=0.0
    values=0
    for idx,case in enumerate(cases):
        data,nw,seg,overlap,k,nf,fs,trend,weighting,width=case
        status,bins,frames,channels=map(int,next(lines_iter).split())
        assert status == 0,(idx,status)
        actual_freq=np.array([float(next(lines_iter)) for _ in range(bins)])
        actual_times=np.array([float(next(lines_iter)) for _ in range(frames)])
        raw=np.array([[float(v) for v in next(lines_iter).split()]
                      for _ in range(bins*frames*channels)])
        actual=raw.reshape((channels,frames,bins,3)).transpose((2,1,0,3))
        fr,t,p,b,r=reference(*case)
        assert actual.shape == (*p.shape,3),(idx,actual.shape,p.shape)
        np.testing.assert_allclose(actual_freq,fr,rtol=3e-13,atol=3e-13,err_msg=f'freq {idx}')
        np.testing.assert_allclose(actual_times,t,rtol=3e-13,atol=3e-13,err_msg=f'time {idx}')
        np.testing.assert_allclose(actual[:,:,:,0],p,rtol=4e-8,atol=5e-10,
                                   err_msg=f'power {idx}')
        np.testing.assert_allclose(actual[:,:,:,1],b,rtol=4e-8,atol=5e-10,
                                   err_msg=f'background {idx}')
        np.testing.assert_allclose(actual[:,:,:,2],r,rtol=4e-8,atol=2e-9,
                                   err_msg=f'whitened {idx}')
        max_power=max(max_power,float(np.max(abs(actual[:,:,:,0]-p))))
        max_background=max(max_background,float(np.max(abs(actual[:,:,:,1]-b))))
        max_ratio=max(max_ratio,float(np.max(abs(actual[:,:,:,2]-r))))
        values+=p.size
    assert next(lines_iter,None) is None
    print(f'PASS: {len(cases)} independent windowed SciPy comparisons '
          f'({values} PSD bins); max abs power {max_power:.3g}, '
          f'background {max_background:.3g}, ratio {max_ratio:.3g}; {mode}')

if __name__ == '__main__':
    main()
