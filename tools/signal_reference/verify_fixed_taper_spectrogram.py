#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent NumPy Fourier sums / Gauss-Legendre band reference.

The band reference integrates the *Fourier transforms* at 96 quadrature nodes;
it does not use SciFort's time-lag sinc kernel. The DPSS reference windows
come from SciPy and the same test includes arbitrary nonorthogonal banks.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.signal.windows import dpss


def reference(x, tapers, freqs, bands, weights, overlap, normalized, average):
    """Independently compute point/window spectra and continuous band integrals."""
    length, n_taper = tapers.shape
    hop = length - overlap
    windows = np.stack([x[start:start+length] for start in
                        range(0, len(x)-length+1, hop)])
    weights = weights / weights.sum()
    energies = np.sum(tapers**2, axis=0) if normalized else np.ones(n_taper)
    weights = weights / energies
    sample_time = np.arange(length)

    def spectral_matrix(f):
        phase = np.exp(-2j*np.pi*f*sample_time)
        z = np.einsum('wlc,lk,l->wkc',windows,tapers,phase,optimize=True)
        return np.einsum('k,wkc,wkd->wcd',weights,z,z.conj(),optimize=True)

    output = [spectral_matrix(f) for f in freqs]
    nodes, quadrature_weights = np.polynomial.legendre.leggauss(96)
    for lower, upper in bands:
        frequencies = 0.5*(upper-lower)*nodes+0.5*(upper+lower)
        factor = 1.0 if average else upper-lower
        integrated = sum(factor*w/2*spectral_matrix(f) for f, w in
                         zip(frequencies, quadrature_weights))
        output.append(integrated)
    return np.stack(output,axis=0), 1.0/np.sum((weights*energies)**2)


def scenarios(cases):
    rng = np.random.default_rng(20261101)
    for i in range(cases):
        length, nchan = 4+i%13, 1+i%5
        n_taper = min(length,1+i%5)
        nw = 2+i%3
        leftover = i%3
        ns = nw*length+leftover
        x = rng.normal(size=(ns,nchan))
        complex_flag = i%3 != 0
        if complex_flag:
            x = x + 1j*rng.normal(size=x.shape)
        if i%2:
            nw_param = min(1.7+.2*(i%4),length*.40)
            tapers = dpss(length,NW=nw_param,Kmax=n_taper,sym=False,norm=2).T
        else:
            tapers = rng.normal(size=(length,n_taper))
        tapers *= rng.uniform(.4,2.,size=(1,n_taper))
        weights = rng.uniform(.1,4.,size=n_taper)
        if i%5==0 and n_taper>1:
            weights[0]=0.0
            tapers[:,0]=0.0
        freqs = rng.uniform(-1.15,1.2,size=i%4)
        bands = np.sort(rng.uniform(-.6,.7,size=((i+1)%4,2)),axis=1)
        if len(bands):
            bands[:,1] += .005
            if i%7==0:
                bands[0]=[-.5,.5]
            if i%9==0:
                bands[0]=[.13,.13001]
        if len(freqs)+len(bands)==0:
            freqs = np.array([.13])
        overlap = 0 if i%3==0 else length//2
        yield (x,tapers,freqs,bands,weights,overlap,bool(i%2),bool((i//2)%2),complex_flag)


def run(executable, cases):
    cases_list = list(scenarios(cases))
    lines = [str(cases)]
    for x,tapers,freqs,bands,weights,overlap,norm,avg,complex_flag in cases_list:
        ns,nc=x.shape
        length,nk=tapers.shape
        lines.append(f'{ns} {nc} {length} {nk} {len(freqs)} {len(bands)} '
                     f'{overlap} {int(norm)} {int(avg)} {int(complex_flag)}')
        z=x.flatten(order='F')
        ri=np.column_stack((z.real,z.imag)).ravel()
        lines.append(' '.join(f'{v:.17g}' for v in ri))
        lines.append(' '.join(f'{v:.17g}' for v in tapers.flatten(order='F')))
        lines.append(' '.join(f'{v:.17g}' for v in weights))
        if len(freqs):
            lines.append(' '.join(f'{v:.17g}' for v in freqs))
        for left,right in bands:
            lines.append(f'{left:.17g} {right:.17g}')
    proc = subprocess.run([str(executable)],input='\n'.join(lines)+'\n',
                          text=True,capture_output=True,check=True)
    rows=proc.stdout.splitlines()
    if len(rows)!=cases:
        raise AssertionError((len(rows),cases,proc.stderr[:2000]))
    max_err = 0.0
    total_matrices = 0
    for i,(row,parameters) in enumerate(zip(rows,cases_list)):
        x,tapers,freqs,bands,weights,overlap,norm,avg,_ = parameters
        expected, effective = reference(x,tapers,freqs,bands,weights,overlap,norm,avg)
        data=np.fromstring(row,sep=' ')
        if len(data)<3 or data[0]!=0:
            raise AssertionError(('invalid Fortran status',i,data[:3]))
        nf,nw,nc,_=expected.shape
        if data[1]!=nw or len(data)!=3+2*nf*nw*nc*nc:
            raise AssertionError(('shape mismatch',i,data[:3],expected.shape))
        measured=data[3:].reshape(nw,nf,nc,nc,2).transpose(1,0,3,2,4)
        measured=measured[...,0]+1j*measured[...,1]
        error=np.max(abs(measured-expected))
        max_err=max(max_err,float(error))
        tolerance=5e-11*(1+np.max(abs(expected)))
        if error>tolerance:
            raise AssertionError(('matrix mismatch',i,error,tolerance,measured,expected))
        if not np.isclose(data[2],effective,atol=1e-12):
            raise AssertionError(('effective taper count',i,data[2],effective))
        if not np.allclose(measured,measured.swapaxes(-2,-1).conj(),atol=2e-11):
            raise AssertionError(('non-Hermitian output',i))
        # Every taper-outer-product and positive frequency-band integral is PSD.
        if np.linalg.eigvalsh(measured).min() < -2e-10*max(1,np.max(abs(measured))):
            raise AssertionError(('spectrum not PSD',i))
        total_matrices+=nf*nw
    print(f'PASS {cases} independent fixed-taper spectrogram scenarios, '
          f'{total_matrices} matrices; max complex entry error {max_err:.3e}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases',type=int,default=50)
    args=parser.parse_args()
    run(args.executable,args.cases)
