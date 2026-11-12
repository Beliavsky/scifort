#!/usr/bin/env python3
"""Reproducible exploratory VAR(1) pointwise bootstrap coverage comparison.

Three-channel VAR(1) simulations compare matrix-Frobenius and target-influence
block selection against unit blocks. An independent high-variance third channel
is deliberately not among the coherence groups. Results are simulation-specific.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.signal.windows import hann


def make_var_matrices(rng, nw: int, nperseg: int, a: np.ndarray,
                      sigma: np.ndarray) -> np.ndarray:
    n = nw*nperseg
    innovations = rng.standard_normal((n+400,3)) @ np.linalg.cholesky(sigma).T
    x = np.zeros_like(innovations)
    for t in range(1,len(x)):
        x[t] = a @ x[t-1] + innovations[t]
    taper = hann(nperseg,sym=False)
    fourier = np.sum(x[400:].reshape(nw,nperseg,3)*taper[None,:,None],axis=1)
    return (fourier[:,:,None] * fourier[:,None,:]).astype(complex)


def true_zero_frequency_coherence(a: np.ndarray,sigma: np.ndarray) -> float:
    transfer = np.linalg.inv(np.eye(3)-a)
    spectral = transfer @ sigma @ transfer.T
    return float(abs(spectral[0,1])**2/(spectral[0,0]*spectral[1,1]))


def run(driver: str, trials: int) -> None:
    rng = np.random.default_rng(20261030)
    nw, nperseg = 48, 12
    scenarios = (
        ('persistent target, white nuisance',
         np.array([[.95,.035,0.],[0.,.94,0.],[0.,0.,.20]])),
        ('weaker target, persistent nuisance',
         np.array([[.60,.08,0.],[0.,.55,0.],[0.,0.,.96]])),
    )
    sigma = np.array([[1.,.5,0.],[.5,1.,0.],[0.,0.,900.]])
    lines = [str(trials*len(scenarios))]
    truth = []
    for _, a in scenarios:
        target = true_zero_frequency_coherence(a,sigma)
        truth.append(target)
        for trial in range(trials):
            windows = make_var_matrices(rng,nw,nperseg,a,sigma)
            lines.append(f'{nw} 3 {1000+trial}')
            for matrix in windows:
                values = matrix.flatten(order='F')
                joined = np.column_stack((values.real,values.imag)).ravel()
                lines.append(' '.join(f'{x:.17g}' for x in joined))
    process = subprocess.run([driver],input='\n'.join(lines)+'\n',text=True,
                             capture_output=True,check=True)
    observed = np.array([[float(x) for x in line.split()]
                         for line in process.stdout.strip().splitlines()])
    assert observed.shape == (trials*len(scenarios),9),observed.shape
    for i,(label,_) in enumerate(scenarios):
        rows = observed[i*trials:(i+1)*trials]
        target = truth[i]
        coverage = [np.mean((rows[:,j] <= target)&(target <= rows[:,j+1]))
                    for j in (3,5,7)]
        width = [np.mean(rows[:,j+1]-rows[:,j]) for j in (3,5,7)]
        print(f'VAR(1) {label}; true squared coherence {target:.4f}; '
              f'{trials} trials, {nw} nonoverlapping {nperseg}-sample windows')
        print(f'  length median: spectral {np.median(rows[:,0]):.1f}, '
              f'influence {np.median(rows[:,1]):.1f} (range spectral '
              f'{int(rows[:,0].min())}..{int(rows[:,0].max())}, influence '
              f'{int(rows[:,1].min())}..{int(rows[:,1].max())})')
        print(f'  empirical nominal-95% coverage: spectral {coverage[0]:.3f}, '
              f'influence {coverage[1]:.3f}, fixed-1 {coverage[2]:.3f}')
        print(f'  mean interval width: spectral {width[0]:.3f}, '
              f'influence {width[1]:.3f}, fixed-1 {width[2]:.3f}')
    print('Exploratory simulation only: percentile intervals are not calibrated.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    parser.add_argument('--trials',type=int,default=48)
    args = parser.parse_args()
    run(args.driver,args.trials)
