"""Independent SciPy checks for multichannel weighted/missing Lomb-Scargle.

Run: python test/compare_signal_irregular_nd_weights_scipy.py checked|optimized
The Fortran driver is compiled afresh to avoid using stale module files.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import lombscargle

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / 'reference' / 'signal_irregular_nd_weights_reference_driver.f90'


def cases():
    return [(n, nf, nc, v) for n in (8, 17, 51, 89)
            for nf in (1, 5, 21) for nc in (1, 3) for v in range(14)]


def build_data(n, nf, nc, variant):
    i = np.arange(1,n+1, dtype=float)
    x = .31*(i-1)+.006*i*i+.17*np.sin(.23*i)
    w = .20+.11*((i.astype(int)*7)%11)
    freqs = .27+.31*np.arange(1,nf+1)+.0007*np.arange(1,nf+1)**2
    y = np.empty((2,n,nc))
    wn = np.empty_like(y)
    for channel in range(1,nc+1):
        for k in (1,2):
            y[k-1,:,channel-1] = (np.cos(.85*x+.23*k*channel)
                +.37*np.sin(1.74*x-.13*channel)+.05*np.cos(.19*i*k)
                +.4*channel+.13*k)
            wn[k-1,:,channel-1] = .18+.10*((i.astype(int)*5+k*3+channel*7)%13)
            if variant in (6,7,8,9):
                missing = (i.astype(int)+3*k+channel)%7==0
                y[k-1,missing,channel-1] = np.nan
                wn[k-1,(i.astype(int)+2*k+channel)%6==0,channel-1]=0
    if variant in (8,9):
        w[i.astype(int)%5==0]=0
    return x, y, freqs, w, wn


def expected(n,nf,nc,variant):
    x,y,freq,w,wn=build_data(n,nf,nc,variant)
    mode = 'amplitude' if variant in (2,5,9,11) else ('normalize' if variant in (1,4,7,12) else 'power')
    out = np.empty((2,nf,nc),dtype=complex if mode=='amplitude' else float)
    for ch in range(nc):
        for k in range(2):
            wchan = wn[k,:,ch] if variant in (3,4,5,6,7,10,12) else (None if variant==13 else w)
            xx=x
            yy=y[k,:,ch]
            omit = variant in (6,7,8,9)
            if omit:
                keep = np.isfinite(yy)
                if wchan is not None:
                    keep &= wchan>0
                xx=xx[keep]
                yy=yy[keep]
                if wchan is not None:
                    wchan=wchan[keep]
            out[k,:,ch]=lombscargle(xx, yy, freq, weights=wchan,
                precenter=(variant in (10,11,12)),
                floating_mean=(variant in (1,2,4,5,7,8,9,12)),normalize=mode)
    return out


def main():
    mode=sys.argv[1] if len(sys.argv)>1 else 'checked'
    if mode not in ('checked','optimized'):
        raise SystemExit('checked or optimized required')
    flags=['-std=f2018','-Werror=implicit-interface']
    flags+=(['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3'])
    workcases=cases()
    with tempfile.TemporaryDirectory(prefix='scifort-lomb-batch-') as tmp:
        wd=Path(tmp)
        command=[os.environ.get('FC','gfortran'),*flags,
                 str(ROOT/'src/kinds.f90'),str(ROOT/'src/signal_irregular.f90'),
                 str(ROOT/'src/signal_irregular_nd.f90'),str(DRIVER),'-o',str(wd/'ref')]
        subprocess.run(command,cwd=wd,capture_output=True,check=True)
        proc=subprocess.run([str(wd/'ref')],cwd=wd,
                input=''.join(' '.join(map(str,c))+'\n' for c in workcases),
                text=True,capture_output=True,check=True)
    lines=iter(proc.stdout.splitlines())
    max_abs=max_scaled=0.0
    count=0
    for case in workcases:
        n,nf,nc,variant=case
        header=next(lines).split()
        assert header==['CASE',str(n),str(nf),str(nc),str(variant),'0'],(case,header)
        is_complex=variant in (2,5,9,11)
        vals=[]
        for _ in range(2*nf*nc):
            line=next(lines)
            vals.append(complex(*map(float,line.split())) if is_complex else float(line))
        actual=np.array(vals).reshape((2,nf,nc),order='F')
        target=expected(*case)
        error=np.max(np.abs(actual-target))
        scaled=error/max(1.,np.max(np.abs(target)))
        if not np.isfinite(scaled) or scaled>1e-10:
            raise AssertionError((case,error,scaled,actual.ravel()[:4],target.ravel()[:4]))
        max_abs=max(max_abs,error)
        max_scaled=max(max_scaled,scaled)
        count+=actual.size
    assert next(lines,None) is None
    print(f'PASS {mode}: {len(workcases)} independent SciPy comparisons, '
          f'{count} spectral values; max abs={max_abs:.3g}, max scaled={max_scaled:.3g}')


if __name__=='__main__':
    main()
