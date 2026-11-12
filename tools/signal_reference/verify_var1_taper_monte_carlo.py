#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Monte Carlo check of finite-taper VAR(1) matrix expectation using stationary starts."""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.signal.windows import hann

from verify_var1_taper_spectrum import direct_reference


def run() -> None:
    rng = np.random.default_rng(20261030)
    a = np.array([[.74,.23,0.0],[-.08,.57,.0],[0.0,0.0,-.35]])
    b = np.array([[1.0,0.0,0.0],[.3,1.1,0.0],[0.2,-.1,.7]])
    sigma = b@b.T
    covariance = solve_discrete_lyapunov(a,sigma)
    taper = hann(12,sym=False)
    frequencies = np.array([0.0,1/6,0.37])
    _, exact = direct_reference(a,sigma,taper,frequencies)
    n = 60000
    x = rng.standard_normal((n,3))@np.linalg.cholesky(covariance).T
    samples = np.empty((12,n,3))
    for t in range(12):
        x = x@a.T + rng.standard_normal((n,3))@b.T
        samples[t] = x
    fourier = np.einsum('tnc,t,ft->fnc',samples,taper,
                        np.exp(-2j*np.pi*frequencies[:,None]*np.arange(12)))
    estimated = np.einsum('fna,fnb->fab',fourier,fourier.conj())/n
    relative = float(np.max(np.linalg.norm(estimated-exact,axis=(1,2)) /
                                  np.linalg.norm(exact,axis=(1,2))))
    assert relative < .035,relative
    print(f'PASS stationary-start VAR(1) Monte Carlo: 60000 windows, '
          f'maximum relative spectral-matrix difference {relative:.3%}')


if __name__ == '__main__':
    run()
