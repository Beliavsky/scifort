#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent VAR(p) discrete Lyapunov and finite-time pairwise quadrature.

Compares expected spectra and exact frequency-band integrals, including the
complex off-diagonal phase. Neither reference calculation calls SciFort.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.signal.windows import hann


def reference(a: np.ndarray, sigma: np.ndarray, taper: np.ndarray,
              freqs: np.ndarray, bands: np.ndarray,
              average: bool = False, normalized: bool = False):
    """Independent O(window^2) time-domain spectral definition."""
    n, _, p = a.shape
    length = len(taper)
    companion = np.zeros((n*p, n*p))
    for k in range(p):
        companion[:n, k*n:(k+1)*n] = a[:, :, k]
    companion[n:, :-n] = np.eye(n*(p-1))
    q = np.zeros_like(companion)
    q[:n, :n] = sigma
    state_gamma = solve_discrete_lyapunov(companion, q)
    powers = [(np.linalg.matrix_power(companion, h) @ state_gamma)[:n,:n]
              for h in range(length)]
    spectrum = np.zeros((len(freqs), n, n), dtype=complex)
    integrated = np.zeros((len(bands), n, n), dtype=complex)
    widths = bands[:, 1]-bands[:, 0]
    mids = (bands[:, 0]+bands[:, 1])/2
    for t in range(length):
        for s in range(length):
            h = t-s
            covariance = powers[h] if h >= 0 else powers[-h].T
            scale = taper[t]*taper[s]
            factors_f = scale*np.exp(-2j*np.pi*freqs*h)
            spectrum += factors_f[:, None, None]*covariance[None, :, :]
            # Integrate the finite Fourier kernel exactly, not a sampled grid.
            factors_b = scale*widths*np.sinc(widths*h)*np.exp(-2j*np.pi*mids*h)
            integrated += factors_b[:,None,None]*covariance[None,:,:]
    if normalized:
        spectrum /= np.sum(taper*taper)
        integrated /= np.sum(taper*taper)
    if average:
        integrated /= widths[:,None,None]
    return state_gamma[:n,:n], spectrum, integrated


def run(executable: str, cases: int):
    rng = np.random.default_rng(20261030 + 2718)
    scenarios = []
    for k in range(cases):
        n = 1+k%4
        p = 1+(k//3)%4
        length = 2+k%17
        nf = 1+k%5
        nb = 1+(k//2)%4
        a = rng.normal(size=(n,n,p))
        # Bound the sum of lag operator norms to ensure VAR stability.
        norm_sum = sum(np.linalg.norm(a[:,:,h],2) for h in range(p))
        a *= rng.uniform(0.20,0.88)/norm_sum
        scales = np.geomspace(0.4,3.0,n)
        b = rng.normal(size=(n,n))*scales[:,None]
        sigma = b @ b.T + np.diag(.4*scales**2)
        taper = hann(length,sym=False) if k%3==0 else rng.uniform(-1,1,length)
        freq = rng.uniform(-.5,.5,nf)
        if k%4==0:
            freq[0] = 0.0
        bands = rng.uniform(-.45,.45,size=(nb,2))
        bands.sort(axis=1)
        bands[:,1] += 0.01
        bands[0] = [-.5,.5] if k%5==0 else bands[0]
        average = bool(k%2)
        normalized = bool((k//2)%2)
        scenarios.append((a,sigma,taper,freq,bands,average,normalized))
    lines=[str(cases)]
    for a,sigma,taper,freq,bands,average,normalized in scenarios:
        n,_,p = a.shape
        lines += [f'{n} {p} {len(taper)} {len(freq)} {len(bands)}',
                  ' '.join(f'{v:.17g}' for v in a.flatten(order='F')),
                  ' '.join(f'{v:.17g}' for v in sigma.flatten(order='F')),
                  ' '.join(f'{v:.17g}' for v in taper),
                  ' '.join(f'{v:.17g}' for v in freq),
                  ' '.join(f'{v:.17g}' for v in bands.flatten(order='C')),
                  f'{int(average)} {int(normalized)}']
    process = subprocess.run([executable],input='\n'.join(lines)+'\n',
                             check=True,text=True,capture_output=True)
    rows=process.stdout.splitlines()
    if len(rows)!=cases:
        raise AssertionError((len(rows),cases,process.stderr))
    max_error=0.0
    matrices=0
    for row,(a,sigma,taper,freq,bands,average,normalized) in zip(rows,scenarios):
        n=a.shape[0]
        values=np.fromstring(row,sep=' ')
        if len(values)<2 or not np.array_equal(values[:2],[0,0]):
            raise AssertionError(('Fortran status',values[:2]))
        expected=2+n*n+2*(len(freq)+len(bands))*n*n
        if len(values)!=expected:
            raise AssertionError((len(values),expected))
        gamma=values[2:2+n*n].reshape((n,n),order='F')
        pos=2+n*n
        def decode(count):
            nonlocal pos
            block=values[pos:pos+2*count*n*n].reshape((count,n*n,2))
            pos+=2*count*n*n
            return np.stack([v.reshape((n,n),order='F')
                             for v in block[...,0]+1j*block[...,1]])
        power=decode(len(freq))
        band_power=decode(len(bands))
        gamma_ref,power_ref,band_ref=reference(a,sigma,taper,freq,bands,average,normalized)
        max_case=max(float(np.max(abs(gamma-gamma_ref))),
                     float(np.max(abs(power-power_ref))),
                     float(np.max(abs(band_power-band_ref))))
        scale=max(1.0,float(np.max(abs(power_ref))),float(np.max(abs(band_ref))),
                  float(np.max(abs(gamma_ref))))
        if max_case>=2.e-10*scale:
            raise AssertionError((max_case,scale,a,sigma,bands))
        max_error=max(max_error,max_case)
        matrices+=len(freq)+len(bands)
    print(f'PASS {cases} independent SciPy/NumPy VAR(p) scenarios, '
          f'{matrices} point and exact-band matrices; max error {max_error:.3e}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases',type=int,default=80)
    args=parser.parse_args()
    run(args.executable,args.cases)
