#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Monte Carlo confirmation of exact finite-window VAR(2) point/band spectra.

The independent reference is from solve_discrete_lyapunov; the sampling check
simulates the recursion directly and does not use any SciFort implementation.
"""
from __future__ import annotations

import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.signal.windows import hann
from verify_varp_taper_spectrum import reference


def run(trials: int = 20000):
    rng=np.random.default_rng(20261030+4242)
    n=3; length=16; burn=130
    a=np.zeros((n,n,2))
    a[:,:,0]=np.array([[.57,.08,-.04],[.06,.42,.12],[.04,-.08,.55]])
    a[:,:,1]=np.array([[-.18,.02,.0],[.01,-.12,.02],[.01,.01,-.21]])
    sigma=np.array([[1.0,.25,.18],[.25,1.25,-.2],[.18,-.2,.9]])
    taper=hann(length,sym=False)
    freq=np.array([-.24,0.,.16,.43])
    bands=np.array([[-.3,-.02],[.0,.25]])
    gamma_ref,point_ref,band_ref=reference(a,sigma,taper,freq,bands)
    # Independently simulate many simultaneously evolving VAR(2) paths.
    cholesky=np.linalg.cholesky(sigma)
    prev=np.zeros((trials,n)); prev2=np.zeros_like(prev)
    sample=np.empty((trials,length,n))
    for t in range(burn+length):
        noise=rng.standard_normal((trials,n))@cholesky.T
        value=prev@a[:,:,0].T+prev2@a[:,:,1].T+noise
        if t>=burn:
            sample[:,t-burn,:]=value
        prev2,prev=prev,value
    def observed(frequencies):
        factors=np.exp(-2j*np.pi*np.outer(frequencies,np.arange(length)))*taper
        z=np.einsum('ft,ktc->kfc',factors,sample,optimize=True)
        return np.einsum('kfi,kfj->fij',z,z.conj(),optimize=True)/trials
    point_obs=observed(freq)
    band_obs=np.zeros_like(band_ref)
    nodes,quad_weights=leggauss(24)
    for b,(lo,hi) in enumerate(bands):
        f=.5*(hi-lo)*nodes+.5*(hi+lo)
        obs=observed(f)
        band_obs[b]=np.tensordot(.5*(hi-lo)*quad_weights,obs,axes=(0,0))
    error=max(float(np.max(abs(point_ref-point_obs))),
              float(np.max(abs(band_ref-band_obs))))
    scale=max(float(np.max(abs(point_ref))),float(np.max(abs(band_ref))))
    if error > .06*scale:
        raise AssertionError((error,scale,'tapered VAR(2) Monte Carlo mismatch'))
    print(f'PASS {trials:,} independent simulated VAR(2) windows; '
          f'point and band max absolute error {error:.4g}, '
          f'relative global scale {error/scale:.3%}')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trials',type=int,default=20000)
    args=parser.parse_args()
    run(args.trials)
