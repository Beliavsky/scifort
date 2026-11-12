"""Compare multitaper spectrogram/F-test against independent SciPy/NumPy.

    python test/compare_signal_multitaper_time_frequency.py [gfortran]

The SciFort calls are made via a compiled Fortran stream driver; reference
calculations use scipy.signal.windows.dpss, NumPy FFT and scipy.stats.f.sf.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal.windows import dpss
from scipy.stats import f as fisher

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / f'{part}.f90' for part in (
    'kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral', 'signal_multitaper')]
DRIVER = ROOT / 'reference' / 'signal_multitaper_time_frequency_reference_driver.f90'


def multitaper_weights(raw, eig, noise, method):
    """Vectorized frequency-dependent Thomson weights from eigenspectra."""
    k, bins = raw.shape
    if method == 'equal':
        return np.ones((k, bins)) / k
    base = eig / eig.sum()
    if method == 'eigen' or k < 3:
        return np.broadcast_to(base[:, None], raw.shape).copy()
    spectral = (eig[:2] @ raw[:2]) / eig[:2].sum()
    small = np.maximum(np.finfo(float).tiny * 1000,
                       np.finfo(float).eps * np.maximum(spectral, noise))
    for _ in range(200):
        den = eig[:, None] * spectral + (1 - eig[:, None]) * noise
        d2 = np.where(den > small, eig[:, None] * (spectral[None, :] /
                                                 np.maximum(den, small)) ** 2, 0)
        total = d2.sum(axis=0)
        q = np.where(total > 0, d2 / np.maximum(total, 1e-300), base[:, None])
        updated = np.sum(q * raw, axis=0)
        if np.all(abs(updated - spectral) / np.maximum(np.maximum(updated, spectral), small) < 1e-10):
            return q
        spectral = updated
    return q


def spectrogram_reference(x, nw, seg, overlap, k, nf, rate, center, mode):
    n, nc = x.shape
    tapers, concentrations = dpss(seg, nw, Kmax=k, sym=False, return_ratios=True)
    real_input = not np.iscomplexobj(x)
    freq = (np.fft.rfftfreq(nf, 1 / rate) if real_input else
            np.fft.fftfreq(nf, 1 / rate))
    starts = list(range(0, n-seg+1, seg-overlap))
    times = (np.asarray(starts) + seg//2) / rate
    result = np.empty((len(freq), len(starts), nc))
    for iframe, start in enumerate(starts):
        for ch in range(nc):
            y = x[start:start+seg, ch].copy()
            if center:
                y -= np.mean(y)
            taper_data = tapers * y[None, :]
            spec = (np.fft.rfft(taper_data, n=nf, axis=1) if real_input else
                    np.fft.fft(taper_data, n=nf, axis=1))
            per = abs(spec)**2 / (rate * (tapers**2).sum(axis=1))[:, None]
            q = multitaper_weights(per, concentrations, np.mean(abs(y)**2) / rate, mode)
            result[:, iframe, ch] = np.sum(per*q, axis=0)
    if real_input:
        if nf % 2 == 0:
            result[1:-1, :, :] *= 2
        else:
            result[1:, :, :] *= 2
    return freq, times, result


def harmonic_reference(x, nw, k, nf, rate, center):
    n = len(x)
    real_input = not np.iscomplexobj(x)
    tapers = dpss(n, nw, Kmax=k, sym=False)
    y = x - np.mean(x) if center else x.copy()
    spectrum = (np.fft.rfft(tapers*y, n=nf, axis=1) if real_input else
                np.fft.fft(tapers*y, n=nf, axis=1))
    response = tapers.sum(axis=1)
    norm = np.sum(response**2)
    amplitude = np.sum(response[:, None]*spectrum, axis=0)/norm
    residual = np.sum(abs(spectrum-response[:, None]*amplitude[None, :])**2, axis=0)
    explained = abs(amplitude)**2 * norm
    # In exact, noiseless line cases the sample residual may underflow to zero.
    fstat = np.divide((k-1)*explained, residual, out=np.zeros_like(explained), where=residual > 0)
    fstat = np.where((residual == 0) & (explained > 0), np.finfo(float).max, fstat)
    pvalue = fisher.sf(fstat,2,2*k-2)
    freq = np.fft.rfftfreq(nf, 1 / rate) if real_input else np.fft.fftfreq(nf, 1 / rate)
    return freq, fstat, pvalue, amplitude


def main():
    fc = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('FC', 'gfortran')
    rng = np.random.default_rng(20261020)
    spec_cases, harmonic_cases = [], []
    spec_lines, harmonic_lines = [], []
    for idx in range(90):
        n = int(rng.integers(24, 93))
        seg = int(rng.integers(8, min(n, 49)+1))
        overlap = int(rng.integers(0, seg))
        nc = int(rng.integers(1, 4))
        k = int(rng.integers(2, min(seg,7)+1))
        nw = float(rng.uniform(1.2, min(3.9, seg/2-0.1)))
        nf = seg+int(rng.integers(0, 17))
        rate = float(rng.uniform(5, 120))
        center = idx % 2 == 0
        mode = ('equal', 'eigen', 'adaptive')[idx % 3]
        ax = 1 if idx % 4 < 2 else 2
        real_input = idx % 2 == 0
        data = rng.normal(size=(n,nc))
        if not real_input:
            data = data + 1j*rng.normal(size=(n,nc))
        if idx % 23 == 0:
            data *= 0
        spec_cases.append((data,nw,seg,overlap,k,nf,rate,center,mode))
        spec_lines.append((n,nc,seg,overlap,k,nf,int(not real_input),int(center),ax,nw,rate,mode))
        for c in range(nc):
            for v in data[:,c]:
                spec_lines.append((v.real, v.imag))
    for idx in range(100):
        n = int(rng.integers(15, 91))
        k = int(rng.integers(2,min(n,8)+1))
        nw = float(rng.uniform(1.2, min(4.5,n/2-0.1)))
        nf = n+int(rng.integers(0, 15))
        rate = float(rng.uniform(1, 80))
        real_input = idx % 2 == 0
        center = idx % 3 == 0
        data = rng.normal(size=n)
        if not real_input:
            data = data+1j*rng.normal(size=n)
        if idx % 10 == 0:
            freq_bin = min(n//3,8)
            data = data*0.01 + 2*np.exp(2j*np.pi*freq_bin*np.arange(n)/n)
            if real_input:
                data = data.real
            nf = n
        if idx % 17 == 0:
            data *= 0
        harmonic_cases.append((data,nw,k,nf,rate,center))
        harmonic_lines.append((n,k,nf,int(not real_input),int(center),nw,rate))
        for v in data:
            harmonic_lines.append((v.real,v.imag))
    # Each input case is followed by its samples, in Fortran column order.
    input_text = f'{len(spec_cases)}\n'+''.join(' '.join(map(str, row))+'\n' for row in spec_lines)
    input_text += f'{len(harmonic_cases)}\n'+''.join(' '.join(map(str, row))+'\n' for row in harmonic_lines)
    with tempfile.TemporaryDirectory(prefix='scifort-tf-') as tmp:
        tmp_path = Path(tmp)
        objects = []
        for src in SOURCES:
            obj = tmp_path / (src.stem + '.o')
            subprocess.run([fc, '-std=f2018', '-O2', '-I', str(tmp_path), '-J', str(tmp_path),
                            '-c',str(src),'-o',str(obj)], cwd=tmp,check=True)
            objects.append(obj)
        exe = tmp_path/'tf_reference'
        subprocess.run([fc, '-std=f2018', '-O2', '-I', str(tmp_path), str(DRIVER),
                        *map(str, objects), '-o',str(exe)], cwd=tmp, check=True)
        process = subprocess.run([str(exe)], input=input_text, text=True,
                                 capture_output=True, check=True, cwd=tmp)
    rows = process.stdout.splitlines()
    index = 0
    errors = []
    for case, (data,nw,seg,overlap,k,nf,rate,center,mode) in enumerate(spec_cases):
        status, bins, frames, nc = map(int,rows[index].split()); index += 1
        assert status == 0, ('spectrogram status', case, status)
        f = np.array([float(v) for v in rows[index:index+bins]]); index += bins
        t = np.array([float(v) for v in rows[index:index+frames]]); index += frames
        power = np.array([float(v) for v in rows[index:index+bins*frames*nc]])
        index += bins*frames*nc
        power = power.reshape((nc,frames,bins)).transpose(2,1,0)
        f0,t0,p0 = spectrogram_reference(data,nw,seg,overlap,k,nf,rate,center,mode)
        np.testing.assert_allclose(f,f0,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(t,t0,rtol=1e-12,atol=1e-12)
        diff = np.max(abs(power-p0))
        errors.append(diff)
        np.testing.assert_allclose(power,p0,rtol=5e-9,atol=3e-10,
                                   err_msg=f'spectrogram case {case}')
    p_errors = []
    f_errors = []
    for case,(data,nw,k,nf,rate,center) in enumerate(harmonic_cases):
        status,bins = map(int,rows[index].split()); index += 1
        assert status == 0, ('F status',case,status)
        values = np.array([[float(v) for v in row.split()] for row in rows[index:index+bins]])
        index += bins
        freq, f0,p0,a0 = harmonic_reference(data,nw,k,nf,rate,center)
        np.testing.assert_allclose(values[:,0],freq,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(values[:,3]+1j*values[:,4],a0,rtol=2e-8,atol=3e-9,
                                   err_msg=f'amplitude case {case}')
        mask = (f0 < 1e7) & (values[:,1] < 1e7)
        if mask.any():
            np.testing.assert_allclose(values[mask,1],f0[mask],rtol=4e-6,atol=2e-8,
                                       err_msg=f'F statistic case {case}')
            f_errors.append(max(abs(values[mask,1]-f0[mask])))
        np.testing.assert_allclose(values[:,2],p0,rtol=2e-6,atol=4e-9,
                                   err_msg=f'F pvalue case {case}')
        p_errors.append(max(abs(values[:,2]-p0)))
    assert index == len(rows), ('unparsed driver output',index,len(rows))
    print(f'PASS: {len(spec_cases)} multitaper spectrogram + {len(harmonic_cases)} harmonic F-test cases')
    print(f'Maximum absolute spectrogram PSD error: {max(errors):.3g}')
    print(f'Maximum absolute harmonic F p-value error: {max(p_errors):.3g}')
    print(f'Maximum absolute finite harmonic F error (F < 1e7): {max(f_errors):.3g}')


if __name__ == '__main__':
    main()
