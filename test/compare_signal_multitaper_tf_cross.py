"""Independent SciPy/NumPy verification of time-frequency multitaper CSD/coherence.

    python test/compare_signal_multitaper_tf_cross.py [gfortran]

SciFort is compiled as a stream driver. The reference uses scipy.signal.windows.dpss
and numpy FFT and independently calculates per-channel Thomson weights.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal.windows import dpss

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / 'src' / f'{name}.f90' for name in (
    'kinds','signal_dpss', 'signal_windows_extra', 'signal_spectral','signal_multitaper')]
DRIVER = ROOT / 'reference' / 'signal_multitaper_tf_cross_reference_driver.f90'


def weights(raw: np.ndarray, eig: np.ndarray, noise: float, mode: str) -> np.ndarray:
    """Calculate taper weights from eigenspectra (not from the Fortran code)."""
    k, bins = raw.shape
    if mode == 'equal':
        return np.ones_like(raw) / k
    base = eig/eig.sum()
    if mode == 'eigen' or k < 3:
        return np.broadcast_to(base[:, None], raw.shape).copy()
    spectrum = (eig[:2] @ raw[:2]) / eig[:2].sum()
    cutoff = np.maximum(np.finfo(float).tiny*1000,
                        np.finfo(float).eps*np.maximum(spectrum, noise))
    q = np.empty_like(raw)
    for _ in range(200):
        den = eig[:, None]*spectrum + (1-eig[:, None])*noise
        w = np.where(den > cutoff, eig[:, None]*(spectrum[None, :]/np.maximum(den,cutoff))**2, 0)
        total = w.sum(axis=0)
        q = np.where(total > 0, w / np.maximum(total,1e-300), base[:,None])
        updated = (q*raw).sum(axis=0)
        if np.all(abs(updated-spectrum)/np.maximum(np.maximum(updated,spectrum),cutoff) < 1e-10):
            return q
        spectrum = updated
    return q


def reference(x,y,nw,seg,ov,k,nf,rate,center,mode):
    real_pair = np.isrealobj(x) and np.isrealobj(y)
    taper,eig = dpss(seg,nw,Kmax=k,sym=False,return_ratios=True)
    freqs = (np.fft.rfftfreq(nf,1/rate) if real_pair else np.fft.fftfreq(nf,1/rate))
    starts = list(range(0,len(x)-seg+1,seg-ov))
    times = (np.array(starts)+seg//2)/rate
    out = np.zeros((len(freqs),len(starts),5))
    for frame,start in enumerate(starts):
        xx = x[start:start+seg].copy()
        yy = y[start:start+seg].copy()
        if center:
            xx -= xx.mean()
            yy -= yy.mean()
        if real_pair:
            fx = np.fft.rfft(taper*xx[None,:],n=nf,axis=1)
            fy = np.fft.rfft(taper*yy[None,:],n=nf,axis=1)
        else:
            fx = np.fft.fft(taper*xx[None,:],n=nf,axis=1)
            fy = np.fft.fft(taper*yy[None,:],n=nf,axis=1)
        norm = np.sqrt(rate*np.sum(taper**2,axis=1))[:,None]
        fx = fx/norm
        fy = fy/norm
        qx = weights(abs(fx)**2,eig,np.mean(abs(xx)**2)/rate,mode)
        qy = weights(abs(fy)**2,eig,np.mean(abs(yy)**2)/rate,mode)
        csd = (np.sqrt(qx*qy)*np.conj(fx)*fy).sum(axis=0)
        pxx = (qx*abs(fx)**2).sum(axis=0)
        pyy = (qy*abs(fy)**2).sum(axis=0)
        if real_pair:
            # Double positive non-DC frequencies, but never the Nyquist bin.
            stop = -1 if nf % 2 == 0 else None
            csd[1:stop] *= 2
            pxx[1:stop] *= 2
            pyy[1:stop] *= 2
        den = pxx*pyy
        coh = np.divide(abs(csd)**2,den,out=np.zeros_like(den),where=den>np.finfo(float).tiny)
        out[:,frame,0] = csd.real
        out[:,frame,1] = csd.imag
        out[:,frame,2] = pxx
        out[:,frame,3] = pyy
        out[:,frame,4] = np.minimum(1,coh)
    return freqs,times,out


def main():
    rng = np.random.default_rng(20261020)
    cases = []
    lines = ['120']
    for i in range(120):
        n = int(rng.integers(25,81))
        seg = int(rng.integers(9,min(n,32)+1))
        ov = int(rng.integers(0,max(1,3*seg//4)))
        k = int(rng.integers(2,min(seg,7)+1))
        nf = seg+int(rng.integers(0,17))
        nw = float(rng.uniform(1.2,min(4.3,seg/2-0.15)))
        rate = float(rng.uniform(1,200))
        kind = i%4
        center = i%2 == 0
        mode = ('equal','eigen','adaptive')[i%3]
        x = rng.normal(size=n)
        y = 0.5*x+rng.normal(size=n)
        if kind in (1,3):
            x = x + 1j*rng.normal(size=n)
        if kind in (1,2):
            y = y + 1j*rng.normal(size=n)
        if i%19 == 0:
            x = x*0
        if i%23 == 0:
            y = y*0
        cases.append((x,y,nw,seg,ov,k,nf,rate,center,mode))
        lines.append(f'{n} {seg} {ov} {k} {nf} {kind} {int(center)} {nw:.17g} {rate:.17g} {mode}')
        for a,b in zip(x,y):
            lines.append(f'{a.real:.17g} {a.imag:.17g} {b.real:.17g} {b.imag:.17g}')
    fc = sys.argv[1] if len(sys.argv)>1 else os.environ.get('FC','gfortran')
    with tempfile.TemporaryDirectory(prefix='scifort-cross-tf-') as tmp:
        work = Path(tmp)
        objects = []
        for src in SOURCES:
            obj = work / (src.stem+'.o')
            subprocess.run([fc,'-std=f2018','-O2','-I',tmp,'-J',tmp,'-c',str(src),'-o',str(obj)],
                           cwd=tmp,check=True)
            objects.append(obj)
        exe = work/'cross_reference'
        subprocess.run([fc,'-std=f2018','-O2','-I',tmp,str(DRIVER),*map(str,objects),'-o',str(exe)],
                       cwd=tmp,check=True)
        process = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,
                                 capture_output=True,cwd=tmp,check=True)
    rows = process.stdout.splitlines()
    cursor = 0
    worst_power = 0
    worst_coherence = 0
    for i,(x,y,nw,seg,ov,k,nf,rate,center,mode) in enumerate(cases):
        status,nfreq,nframe = map(int,rows[cursor].split()); cursor+=1
        assert status == 0, ('invalid status',i,status)
        f = np.array([float(v) for v in rows[cursor:cursor+nfreq]]); cursor+=nfreq
        t = np.array([float(v) for v in rows[cursor:cursor+nframe]]); cursor+=nframe
        vals = np.array([[float(v) for v in line.split()] for line in rows[cursor:cursor+nfreq*nframe]])
        cursor += nfreq*nframe
        vals = vals.reshape((nframe,nfreq,5)).transpose(1,0,2)
        ref_freq,ref_time,ref_values = reference(x,y,nw,seg,ov,k,nf,rate,center,mode)
        np.testing.assert_allclose(f,ref_freq,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(t,ref_time,rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(vals[...,:4],ref_values[...,:4],rtol=5e-8,atol=1e-9,
                                   err_msg=f'cross spectrum {i}')
        np.testing.assert_allclose(vals[...,4],ref_values[...,4],rtol=5e-8,atol=1e-9,
                                   err_msg=f'time-frequency coherence {i}')
        worst_power = max(worst_power,float(np.max(abs(vals[...,:4]-ref_values[...,:4]))))
        worst_coherence = max(worst_coherence,float(np.max(abs(vals[...,4]-ref_values[...,4]))))
    assert cursor == len(rows), ('remaining output',len(rows)-cursor)
    print(f'PASS: {len(cases)} cross-spectrogram/coherence cases against independent SciPy/NumPy reference')
    print(f'Maximum absolute spectral discrepancy: {worst_power:.3g}')
    print(f'Maximum absolute coherence discrepancy: {worst_coherence:.3g}')


if __name__ == '__main__':
    main()
