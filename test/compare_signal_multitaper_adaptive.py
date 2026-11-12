"""Independent SciPy DPSS, NumPy FFT and Student-t oracle for adaptive spectra.

    python test/compare_signal_multitaper_adaptive.py [gfortran]

Uses fixed seeded signals and its own frequency-dependent Thomson iterative
solver and delete-one log-jackknife.  No external Fortran dependencies.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.signal.windows import dpss
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / s for s in ('kinds.f90', 'signal_dpss.f90', 'signal_windows_extra.f90', 'signal_spectral.f90',
          'signal_multitaper.f90')]
DRIVER = ROOT / 'reference' / 'signal_multitaper_adaptive_reference_driver.f90'


def weights(p, lam, noise):
    """Thomson's adaptive per-frequency weights, independently vectorized."""
    nt, nf = p.shape
    if nt < 3:
        return np.repeat((lam / lam.sum())[:, None], nf, axis=1)
    s = (lam[0] * p[0] + lam[1] * p[1]) / (lam[0] + lam[1])
    cutoff = np.maximum(np.finfo(float).tiny * 1000,
                        np.finfo(float).eps * np.maximum(s, noise))
    for iteration in range(200):
        den = lam[:, None] * s + (1-lam[:, None]) * noise
        d2 = np.where(den > cutoff,
                      lam[:, None] * (s / np.maximum(den, cutoff))**2, 0.0)
        total = d2.sum(axis=0)
        q = np.where(total > 0, d2 / np.maximum(total, 1e-300),
                     lam[:, None] / lam.sum())
        snext = (p*q).sum(axis=0)
        if np.all(np.abs(snext-s)/np.maximum(np.maximum(snext,s),cutoff) < 1e-10):
            return q
        s = snext
    return q


def oracle(data, nw, k, fs, nf, center, conf):
    n, nc = data.shape
    real_input = not np.iscomplexobj(data)
    taper, lam = dpss(n, nw, Kmax=k, sym=False, return_ratios=True)
    if center:
        data = data - data.mean(axis=0, keepdims=True)
    windows = data[None,:,:] * taper[:,:,None]
    X = (np.fft.rfft(windows, n=nf, axis=1) if real_input else
         np.fft.fft(windows, n=nf, axis=1))
    # X has (taper, frequency, channel) axes, normalized per taper.
    X /= np.sqrt(fs * (taper**2).sum(axis=1))[:,None,None]
    nk = X.shape[1]
    q = np.empty((k,nk,nc))
    for ch in range(nc):
        p = abs(X[:,:,ch])**2
        var = np.mean(abs(data[:,ch])**2) / fs
        q[:,:,ch] = weights(p, lam, var)
    W = X * np.sqrt(q)
    matrix = np.einsum('kfc,kfd->fcd', np.conj(W), W)
    if real_input:
        if nf % 2 == 0:
            matrix[1:-1] *= 2
        else:
            matrix[1:] *= 2
    freq = (np.fft.rfftfreq(nf, d=1/fs) if real_input else
            np.fft.fftfreq(nf, d=1/fs))
    psd = matrix[:,0,0].real
    raw = abs(X[:,:,0])**2
    noise = np.mean(abs(data[:,0])**2)/fs
    jk = []
    for drop in range(k):
        keep = np.arange(k) != drop
        qleft = weights(raw[keep], lam[keep], noise)
        jk.append(np.log(np.maximum((qleft * raw[keep]).sum(axis=0),
                                    np.finfo(float).tiny)))
    jk = np.array(jk)
    jkvar = (k-1) / k * ((jk-jk.mean(axis=0))**2).sum(axis=0)
    stderr = np.sqrt(jkvar)
    critical = student_t.ppf((1+conf)/2, k-1)
    lower = psd * np.exp(-critical * stderr)
    upper = psd * np.exp(critical * stderr)
    dof = 2 / (q[:,:,0]**2).sum(axis=0)
    return freq, matrix, (psd, lower, upper, stderr, dof)


def main():
    fc = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('FC', 'gfortran')
    rng = np.random.default_rng(20261020)
    cases, lines = [], []
    for j in range(140):
        n = int(rng.integers(13,69))
        nc = int(rng.integers(1,4))
        nf = n + int(rng.integers(0,24))
        nw = float(rng.uniform(1.2,min(4.5,n/2-0.1)))
        k = int(rng.integers(3,min(n,7)+1))
        fs = float(rng.uniform(1.,80.))
        real_input = j % 2 == 0
        center = j % 3 != 0
        ax = 1 if j % 4 < 2 else 2
        conf = [0.8,0.9,0.95,0.99][j%4]
        data = rng.normal(size=(n,nc))
        if not real_input:
            data = data + 1j*rng.normal(size=(n,nc))
        if j%11==0:
            data[:,nc-1] = data[:,0] * (1.5 if real_input else (0.5-1.7j))
        if j%27==0:
            data = data * 0
        cases.append((data,nw,k,fs,nf,center,conf))
        lines.append(f'{n} {nc} {nf} {k} {nw:.17g} {fs:.17g} '
                     f'{int(not real_input)} {ax} {int(center)} {conf}\n')
        for c in range(nc):
            for v in data[:,c]:
                lines.append(f'{v.real:.17g} {v.imag:.17g}\n')
    from contextlib import nullcontext
    external = os.getenv('SCIFORT_MT_REF_EXE')
    with (nullcontext() if external else tempfile.TemporaryDirectory(prefix='scifort-adaptive-')) as temp:
        tmp = Path(temp) if temp else None
        exe = Path(external) if external else tmp / ('reference.exe' if os.name=='nt' else 'reference')
        if not external:
            subprocess.run([fc,'-std=f2018','-O2','-fcheck=all',
                            '-Werror=implicit-interface', *map(str,SOURCES),
                            str(DRIVER),'-o',str(exe)], cwd=tmp,check=True)
        p = subprocess.run([str(exe)],input=''.join(lines),text=True,
                           capture_output=True,check=True)
    got = iter(p.stdout.splitlines())
    max_matrix = max_interval = max_dof = 0.
    for i,(data,nw,k,fs,nf,center,conf) in enumerate(cases):
        freq, matrix, interval = oracle(data,nw,k,fs,nf,center,conf)
        header = next(got).split()
        assert header == ['CASE',str(len(freq)),str(data.shape[1]),'0','0'],(i,header)
        f = np.fromstring(next(got),sep=' ')
        assert np.allclose(f,freq,atol=1e-12)
        result = np.zeros_like(matrix)
        for d in range(data.shape[1]):
            for c in range(data.shape[1]):
                raw = np.fromstring(next(got),sep=' ')
                result[:,c,d] = raw[::2]+1j*raw[1::2]
        vals = [np.fromstring(next(got),sep=' ') for _ in range(5)]
        err = np.max(abs(result-matrix))
        max_matrix = max(max_matrix,err)
        if not np.allclose(result,matrix,atol=3e-10,rtol=4e-8):
            raise AssertionError(f'case {i} matrix error {err}')
        for index,(observed,expected) in enumerate(zip(vals,interval)):
            diff = np.max(abs(observed-expected))
            if index==4: max_dof=max(max_dof,diff)
            else: max_interval=max(max_interval,diff)
            if not np.allclose(observed,expected,atol=1e-8,rtol=5e-7):
                raise AssertionError(f'case {i} interval[{index}] error {diff}, maxval {np.max(abs(expected))}')
    print(f'Adaptive multitaper SciPy-DPSS/NumPy FFT/Student-t comparisons: {len(cases)} passed; '
          f'max matrix abs error {max_matrix:.3g}, CI/SE max error {max_interval:.3g}, '
          f'equivalent-DOF max error {max_dof:.3g}')


if __name__ == '__main__':
    main()
