#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check VAR(p) conditional canonical mode/population band targets independently."""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.signal.windows import hann
from verify_varp_taper_spectrum import reference


def canonical(spectra):
    out=[]
    for s in spectra:
        a=np.array([0,1]); b=np.array([2,3]); c=np.array([4])
        aa=s[np.ix_(a,a)]-s[np.ix_(a,c)]@np.linalg.solve(s[np.ix_(c,c)],s[np.ix_(c,a)])
        bb=s[np.ix_(b,b)]-s[np.ix_(b,c)]@np.linalg.solve(s[np.ix_(c,c)],s[np.ix_(c,b)])
        ab=s[np.ix_(a,b)]-s[np.ix_(a,c)]@np.linalg.solve(s[np.ix_(c,c)],s[np.ix_(c,b)])
        la=cholesky(aa,lower=True)
        lb=cholesky(bb,lower=True)
        whitened=solve_triangular(la,ab,lower=True)
        whitened=solve_triangular(lb,whitened.conj().T,lower=True).conj().T
        modes=np.linalg.svd(whitened,compute_uv=False)**2
        aggregate=-np.expm1(np.log1p(-modes).sum())
        out.append((aggregate,*modes))
    return np.asarray(out)


def run(executable,cases):
    rng=np.random.default_rng(20261030+771)
    scenarios=[]
    for k in range(cases):
        n=5; p=1+k%3; l=4+k%14; nf=1+k%4; nb=1+(k//2)%3
        a=rng.normal(size=(n,n,p))
        a*=rng.uniform(.22,.78)/sum(np.linalg.norm(a[:,:,h],2) for h in range(p))
        z=rng.normal(size=(n,n))
        sigma=z@z.T+np.eye(n)*.8
        taper=hann(l,sym=False) if k%2==0 else rng.uniform(.1,1,l)
        freq=rng.uniform(-.4,.4,nf)
        bands=np.column_stack([rng.uniform(-.4,-.02,nb),rng.uniform(.03,.45,nb)])
        scenarios.append((a,sigma,taper,freq,bands))
    lines=[str(cases)]
    for a,sigma,taper,freq,bands in scenarios:
        lines += [f'{a.shape[-1]} {len(taper)} {len(freq)} {len(bands)}',
                  ' '.join(f'{v:.17g}' for v in a.flatten(order='F')),
                  ' '.join(f'{v:.17g}' for v in sigma.flatten(order='F')),
                  ' '.join(f'{v:.17g}' for v in taper),
                  ' '.join(f'{v:.17g}' for v in freq),
                  ' '.join(f'{v:.17g}' for v in bands.flatten())]
    process=subprocess.run([executable],input='\n'.join(lines)+'\n',
                           text=True,capture_output=True,check=True)
    rows=process.stdout.splitlines()
    if len(rows)!=cases:
        raise AssertionError((len(rows),cases,process.stderr))
    max_error=0.0
    compared=0
    for row,(a,sigma,taper,freq,bands) in zip(rows,scenarios):
        vals=np.fromstring(row,sep=' ')
        if vals.size<2 or not np.array_equal(vals[:2],[0,0]):
            raise AssertionError(('bad Fortran status', vals[:2]))
        nf,nb=len(freq),len(bands)
        if len(vals)!=2+3*(nf+nb):
            raise AssertionError((len(vals),nf,nb))
        _,point,integrated=reference(a,sigma,taper,freq,bands)
        expected=np.vstack((canonical(point),canonical(integrated)))
        measured=vals[2:].reshape((-1,3))
        error=float(np.max(abs(measured-expected)))
        if error>3e-9:
            raise AssertionError((error,measured,expected))
        max_error=max(max_error,error)
        compared+=(nf+nb)*3
    print(f'PASS {cases} independent conditional canonical VAR(p) cases, '
          f'{compared} scalar targets; max error {max_error:.3e}')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases',type=int,default=50)
    args=parser.parse_args()
    run(args.executable,args.cases)
