"""Compare SciFort harmonic-line multiplicity corrections against SciPy.

    python test/compare_signal_multitaper_harmonic_scan.py [gfortran]

The scan is a SciFort extension. The underlying F-statistic is independently
checked against SciPy DPSS + NumPy FFT; corrections use scipy.stats and
statsmodels.stats.multitest.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal.windows import dpss
from scipy.stats import f as fisher, false_discovery_control
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT/'src'/f'{s}.f90' for s in (
    'kinds','signal_dpss', 'signal_windows_extra', 'signal_spectral','signal_multitaper')]
DRIVER = ROOT/'reference'/'signal_multitaper_harmonic_scan_reference_driver.f90'
METHODS = ('none','bonferroni','bh','by','holm')


def harmonic_reference(x,nw,k,nf,rate,center):
    n = len(x)
    tapers = dpss(n,nw,Kmax=k,sym=False)
    y = x-np.mean(x) if center else x.copy()
    spectra = np.fft.rfft(tapers*y,n=nf,axis=1) if np.isrealobj(y) else np.fft.fft(tapers*y,n=nf,axis=1)
    response = tapers.sum(axis=1)
    norm = sum(response**2)
    amplitude = np.sum(response[:,None]*spectra,axis=0)/norm
    explained = abs(amplitude)**2*norm
    residual = np.sum(abs(spectra-response[:,None]*amplitude[None,:])**2,axis=0)
    fstat = np.divide((k-1)*explained,residual,out=np.zeros_like(explained),where=residual>0)
    fstat[(residual==0)&(explained>0)] = np.finfo(float).max
    pvalue = fisher.sf(fstat,2,2*k-2)
    freq = (np.fft.rfftfreq(nf,1/rate) if np.isrealobj(y) else np.fft.fftfreq(nf,1/rate))
    return freq,pvalue,fstat,amplitude


def adjust(p,mask,method):
    adjusted = np.ones(len(p))
    m = np.sum(mask)
    if m == 0:
        return adjusted
    if method == 'none':
        adjusted[mask] = p[mask]
    elif method == 'bonferroni':
        adjusted[mask] = np.minimum(1,m*p[mask])
    elif method == 'holm':
        adjusted[mask] = multipletests(p[mask],method='holm')[1]
    else:
        adjusted[mask] = false_discovery_control(p[mask],method=method)
    return adjusted


def main():
    fc = sys.argv[1] if len(sys.argv)>1 else os.getenv('FC','gfortran')
    rng = np.random.default_rng(20261020)
    cases = []
    lines = []
    for case_id in range(225):
        n = int(rng.integers(16,88))
        k = int(rng.integers(2,min(7,n)+1))
        nw = float(rng.uniform(1.1,min(4.0,n/2-0.1)))
        nf = n+int(rng.integers(0,31))
        rate = float(rng.uniform(1,100))
        complex_input = case_id%2==1
        detrend = case_id%3==0
        method = METHODS[case_id%len(METHODS)]
        edges = case_id%5==0
        alpha = float(rng.uniform(0.001,0.25))
        x = rng.normal(size=n)
        if complex_input:
            x = x+1j*rng.normal(size=n)
        if case_id%9==0:
            x = x + 2.0*np.exp(2j*np.pi*7*np.arange(n)/nf).real
        if case_id%31==0:
            x *= 0.0
        bins_count = (nf if complex_input else nf//2+1)
        selected = None
        if case_id%6==0:
            count = int(rng.integers(0,min(17,bins_count)+1))
            selected = sorted(rng.choice(bins_count,size=count,replace=False).tolist())
        cases.append((x,nw,k,nf,rate,detrend,method,edges,alpha,selected))
        lines.append(f'{n} {k} {nf} {int(complex_input)} {int(detrend)} {METHODS.index(method)} '
                     f'{int(edges)} {len(selected) if selected is not None else -1} {nw} {rate} {alpha}')
        lines.extend(f'{v.real:.17g} {v.imag:.17g}' for v in x)
        if selected:
            lines.append(' '.join(str(i+1) for i in selected))
    stdin = str(len(cases))+'\n'+'\n'.join(lines)+'\n'
    with tempfile.TemporaryDirectory(prefix='scifort-harmonic-scan-') as tmp:
        objects=[]
        for src in SOURCES:
            obj=Path(tmp)/(src.stem+'.o')
            subprocess.run([fc,'-O2','-std=f2018','-I',tmp,'-J',tmp,'-c',str(src),'-o',str(obj)],
                           cwd=tmp,check=True)
            objects.append(str(obj))
        exe=Path(tmp)/'harmonic_scan_ref'
        subprocess.run([fc,'-O2','-std=f2018','-I',tmp,str(DRIVER),*objects,'-o',str(exe)],
                       check=True,cwd=tmp)
        completed=subprocess.run([str(exe)],input=stdin,text=True,capture_output=True,check=True,cwd=tmp)
    rows=completed.stdout.splitlines()
    pos=0
    max_p_error=0.0
    max_adj_error=0.0
    for ic,(x,nw,k,nf,rate,center,method,edges,alpha,selected) in enumerate(cases):
        status,num,m=map(int,rows[pos].split());pos+=1
        assert status==0 and num==(nf if np.iscomplexobj(x) else nf//2+1), (ic,status,num)
        output=np.array([[float(s) for s in row.split()] for row in rows[pos:pos+num]])
        pos+=num
        f,p,fstat,amplitude=harmonic_reference(x,nw,k,nf,rate,center)
        mask=np.zeros(num,dtype=bool)
        if selected is not None:
            mask[selected]=True
        else:
            mask[:]=True
            if not edges:
                mask[0]=False
                if not np.iscomplexobj(x) and nf%2==0:
                    mask[-1]=False
        assert m==int(mask.sum()),(ic,m,mask.sum())
        np.testing.assert_allclose(output[:,0],f,atol=4e-13,rtol=2e-13)
        np.testing.assert_allclose(output[:,1],p,atol=5e-9,rtol=1e-6,err_msg=f'pointwise {ic}')
        expectation=adjust(p,mask,method)
        np.testing.assert_allclose(output[:,2],expectation,atol=5e-9,rtol=2e-6,
                                   err_msg=f'adjusted {ic}, method={method}')
        np.testing.assert_array_equal(output[:,6].astype(bool),mask & (expectation<=alpha))
        np.testing.assert_allclose(output[:,4]+1j*output[:,5],amplitude,atol=1e-9,rtol=2e-8)
        max_p_error=max(max_p_error,float(np.max(abs(output[:,1]-p))))
        max_adj_error=max(max_adj_error,float(np.max(abs(output[:,2]-expectation))))
    assert pos==len(rows)
    print(f'PASS: {len(cases)} SciPy harmonic scans, max raw p error {max_p_error:.3g}, '
          f'correction error {max_adj_error:.3g}')


if __name__=='__main__':
    main()
