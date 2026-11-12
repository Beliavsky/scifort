"""Independent SciPy/NumPy checks for colored-noise multitaper background.

The median-normalized PSD is a SciFort extension, not a SciPy API.
The reference uses scipy.signal.windows.dpss, NumPy FFT, scipy.signal.detrend,
and scipy.ndimage.median_filter.

Usage: python test/compare_signal_multitaper_noise_floor.py [gfortran] [O0|O3]
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

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / (n + '.f90') for n in
           ('kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral', 'signal_multitaper')]
DRIVER = ROOT / 'reference' / 'signal_multitaper_noise_floor_reference_driver.f90'


def reference(x: np.ndarray, nw: float, k: int, nf: int, fs: float,
              trend: int, weight: int, width: int):
    if trend == 1:
        x = x - x.mean()
    elif trend == 2:
        x = scipy_detrend(x, type='linear')
    tapers, eig = dpss(len(x), nw, Kmax=k, sym=False, return_ratios=True)
    ffts = (np.fft.fft(tapers*x[None, :], n=nf, axis=1)
            if np.iscomplexobj(x) else np.fft.rfft(tapers*x[None, :], n=nf, axis=1))
    raw = np.abs(ffts)**2 / (fs*np.sum(tapers*tapers, axis=1)[:, None])
    power = np.average(raw, axis=0, weights=eig if weight else None)
    if np.isrealobj(x):
        if nf % 2 == 0:
            power[1:-1] *= 2
        else:
            power[1:] *= 2
        freq = np.fft.rfftfreq(nf, 1/fs)
    else:
        freq = np.fft.fftfreq(nf, 1/fs)
    background = median_filter(power, size=width, mode='reflect' if np.isrealobj(x) else 'wrap')
    ratio = np.divide(power, background, out=np.zeros_like(power), where=background>0)
    ratio[(power>0) & (background==0)] = np.finfo(float).max
    return freq, power, background, ratio


def main():
    fc = sys.argv[1] if len(sys.argv) > 1 else os.getenv('FC', 'gfortran')
    mode = sys.argv[2] if len(sys.argv) > 2 else 'O0'
    assert mode in ('O0', 'O3')
    rng = np.random.default_rng(20261026)
    cases = []
    input_lines = []
    for case_id in range(144):
        n = int(rng.integers(25, 112))
        k = int(rng.integers(2, 7))
        nw = float(rng.uniform(1.1, min(4.8, n/2-0.5)))
        nf = n + int(rng.integers(0, 70))
        fs = float(rng.uniform(1.0, 100.0))
        complex_flag = case_id % 2
        trend = case_id % 3
        weight = (case_id//3) % 2
        nbins = nf if complex_flag else nf//2+1
        width = int(rng.integers(1, min(12, nbins//2) + 1))*2+1
        if width > nbins:
            width = nbins if nbins % 2 else nbins-1
        # Stable seeded colored AR(1) background with several occasional lines.
        excitation = rng.normal(size=n) + (1j*rng.normal(size=n) if complex_flag else 0)
        x = np.empty(n, dtype=np.complex128 if complex_flag else np.float64)
        x[0] = excitation[0]
        for i in range(1,n):
            x[i] = 0.82*x[i-1] + excitation[i]
        if case_id % 3 == 0:
            line = 0.8*np.exp(2j*np.pi*(case_id%19+3)*np.arange(n)/nf)
            x = x+line if complex_flag else x+line.real
        if case_id % 23 == 0:
            x[:] = 0
        cases.append((x,nw,k,nf,fs,trend,weight,width))
        input_lines.append(f'{n} {k} {nf} {complex_flag} {trend} {weight} {width} '
                           f'{nw:.17g} {fs:.17g}')
        input_lines.extend(f'{v.real:.17g} {v.imag:.17g}' for v in x)
    stdin = str(len(cases))+'\n'+'\n'.join(input_lines)+'\n'
    with tempfile.TemporaryDirectory(prefix='scifort-noise-floor-') as tmpdir:
        tmp = Path(tmpdir)
        objects = []
        flags = ['-O3'] if mode == 'O3' else ['-O0', '-g', '-fcheck=all']
        for src in SOURCES:
            obj = tmp / (src.stem+'.o')
            subprocess.run([fc, '-std=f2018', '-Werror=implicit-interface', *flags,
                            '-I', str(tmp), '-J', str(tmp), '-c', str(src), '-o', str(obj)],
                           check=True, cwd=tmp)
            objects.append(str(obj))
        exe = tmp / 'noise_floor_reference'
        subprocess.run([fc, '-std=f2018', *flags, '-I', str(tmp), str(DRIVER),
                        *objects, '-o', str(exe)], check=True, cwd=tmp)
        result = subprocess.run([str(exe)], input=stdin, text=True, capture_output=True,
                                check=True, cwd=tmp)
    lines = iter(result.stdout.splitlines())
    max_psd_error = max_bg_error = max_ratio_error = 0.0
    values = 0
    for index,(x,nw,k,nf,fs,trend,weight,width) in enumerate(cases):
        status, size = (int(v) for v in next(lines).split())
        assert status == 0, (index, status)
        actual = np.array([[float(s) for s in next(lines).split()] for _ in range(size)])
        expected = np.column_stack(reference(x,nw,k,nf,fs,trend,weight,width))
        assert actual.shape == expected.shape, (index,actual.shape,expected.shape)
        np.testing.assert_allclose(actual[:,0],expected[:,0],rtol=2e-13,atol=2e-13)
        np.testing.assert_allclose(actual[:,1:],expected[:,1:],rtol=3e-8,atol=2e-9,
                                   err_msg=f'colored spectrum {index}')
        max_psd_error = max(max_psd_error,float(np.max(np.abs(actual[:,1]-expected[:,1]))))
        max_bg_error = max(max_bg_error,float(np.max(np.abs(actual[:,2]-expected[:,2]))))
        max_ratio_error = max(max_ratio_error,float(np.max(np.abs(actual[:,3]-expected[:,3]))))
        values += size
    assert next(lines,None) is None
    print(f'PASS: {len(cases)} SciPy colored-noise comparisons ({values} bins); '
          f'max abs PSD error {max_psd_error:.3g}, background error {max_bg_error:.3g}, '
          f'whitened error {max_ratio_error:.3g}; {mode}')


if __name__ == '__main__':
    main()
