#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Monte Carlo canonical-bootstrap coverage for the matching Hann-smoothed VAR(1) target.

The target is C(E[Z Z*]) for each finite tapered Fourier window Z,
NOT C(S(0)), C(E[sample covariance]), or E[C(sample covariance)].
Intervals remain exploratory because short-window and dependent-block bias
can cause substantial undercoverage regardless of the target correction.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.signal.windows import hann

from study_canonical_influence_var import make_var_matrices, true_zero_frequency_coherence
from verify_var1_taper_spectrum import direct_reference


METHODS = ('percentile', 'basic', 'bc', 'bca', 'studentized')


def wilson(k: int, n: int) -> tuple[float, float]:
    """Wilson score interval for a binomial Monte Carlo coverage proportion."""
    z = 1.959963984540054
    if n == 0:
        return float('nan'), float('nan')
    p = k / n
    denominator = 1.0 + z*z/n
    center = (p + z*z/(2*n))/denominator
    halfwidth = z/denominator * np.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return float(center-halfwidth), float(center+halfwidth)


def scenarios():
    """Fixed VAR coefficients and innovations; third channel is a nuisance."""
    sigma = np.array([[1.0,.5,0.0],[.5,1.0,0.0],[0.0,0.0,900.0]])
    return (
        ('persistent target / white nuisance',
         np.array([[.95,.035,0.0],[0.0,.94,0.0],[0.0,0.0,.20]]),sigma),
        ('moderate target / persistent nuisance',
         np.array([[.60,.08,0.0],[0.0,.55,0.0],[0.0,0.0,.96]]),sigma),
        ('lag-induced target (independent innovations)',
         np.array([[.77,.32,0.0],[0.0,.48,0.0],[0.0,0.0,.3]]),
         np.diag([1.0,1.0,900.0])),
        ('opposite-sign autoregression',
         np.array([[.85,-.16,0.0],[0.0,-.62,0.0],[0.0,0.0,.7]]),sigma),
    )


def true_window_coherence(a: np.ndarray, sigma: np.ndarray, nperseg: int) -> float:
    taper = hann(nperseg,sym=False)
    _,spectrum = direct_reference(a,sigma,taper,np.array([0.0]))
    s = spectrum[0]
    return float(abs(s[0,1])**2/(s[0,0].real*s[1,1].real))


def run(driver: str, trials: int = 64, nw: int = 48,
        nperseg: int = 12, seed: int = 20261030) -> None:
    if trials < 1 or nw < 16 or nperseg < 3:
        raise ValueError('require trials>=1, windows>=16 and samples/window>=3')
    rng = np.random.default_rng(seed)
    models = scenarios()
    lines = [str(trials*len(models))]
    for _,a,sigma in models:
        for trial in range(trials):
            windows = make_var_matrices(rng,nw,nperseg,a,sigma)
            lines.append(f'{nw} 3 {1000+trial}')
            for matrix in windows:
                values = matrix.flatten(order='F')
                pairs = np.column_stack((values.real,values.imag)).ravel()
                lines.append(' '.join(f'{x:.17g}' for x in pairs))
    process = subprocess.run([driver],input='\n'.join(lines)+'\n',text=True,
                             capture_output=True,check=True)
    observations = [np.fromstring(row,sep=' ') for row in process.stdout.splitlines()]
    if len(observations) != trials*len(models) or any(len(row) != 13 for row in observations):
        raise AssertionError((len(observations),[len(r) for r in observations[:3]],
                              process.stderr))
    print(f'VAR(1) finite-window target study: {len(models)} models x {trials} datasets; '
          f'{nw} nonoverlapping windows, {nperseg} Hann samples/window, '
          '199 circular moving-block draws, block length 4')
    print('Coverage percentages are conditional on the interval being valid; '
          'Wilson intervals describe Monte Carlo uncertainty only.')
    for index,(label,a,sigma) in enumerate(models):
        target = true_window_coherence(a,sigma,nperseg)
        zero_target = true_zero_frequency_coherence(a,sigma)
        rows = np.vstack(observations[index*trials:(index+1)*trials])
        print(f'\n{label}: C(E[window spectrum])={target:.5f}, '
              f'C(S(0))={zero_target:.5f}, '
              f'difference={target-zero_target:+.5f}')
        # 10 endpoint columns, 3 validity columns; all methods use the same
        # point estimate and seeded bootstrap replicates.
        for j, method in enumerate(METHODS):
            lo,hi = rows[:,2*j],rows[:,2*j+1]
            finite = np.isfinite(lo)&np.isfinite(hi)&(lo<=hi)
            valid = (rows[:,10+j-2] == 1.0) if j >= 2 else np.ones(trials,dtype=bool)
            valid &= finite
            n = int(valid.sum())
            if not n:
                print(f'  {method:11s} no valid intervals')
                continue
            covers = int(np.sum((lo[valid]<=target)&(target<=hi[valid])))
            low_ci,up_ci = wilson(covers,n)
            width = float(np.mean(hi[valid]-lo[valid]))
            print(f'  {method:11s} coverage {100*covers/n:5.1f}% '
                  f'(MC 95% Wilson {100*low_ci:4.1f}..{100*up_ci:4.1f}%), '
                  f'valid {n}/{trials}, mean width {width:.4f}')
    print('\nStatistical note: pointwise resampling intervals can still undercover '
          'a correctly defined population functional. No calibrated-95% claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    parser.add_argument('--trials',type=int,default=64)
    parser.add_argument('--windows',type=int,default=48)
    parser.add_argument('--samples-per-window',type=int,default=12)
    parser.add_argument('--seed',type=int,default=20261030)
    args = parser.parse_args()
    run(args.driver,args.trials,args.windows,args.samples_per_window,args.seed)
