#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent SciPy Lyapunov and direct (t,s) window covariance reference."""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.signal.windows import hann


def direct_reference(a: np.ndarray, sigma: np.ndarray,
                     taper: np.ndarray, frequencies: np.ndarray):
    """Sum E[x_t x_s'] over the individual taper sample pairs."""
    gamma = solve_discrete_lyapunov(a, sigma)
    l = taper.size
    n = a.shape[0]
    spectrum = np.zeros((len(frequencies), n, n), complex)
    powers = [np.linalg.matrix_power(a, k) @ gamma for k in range(l)]
    for t in range(l):
        for s in range(l):
            covariance = powers[t-s] if t >= s else powers[s-t].T
            factors = taper[t] * taper[s] * np.exp(-2j*np.pi*frequencies*(t-s))
            spectrum += factors[:, None, None]*covariance[None, :, :]
    return gamma, spectrum


def run(executable: str, cases: int) -> None:
    rng = np.random.default_rng(20261030)
    scenarios = []
    for k in range(cases):
        n = 1 + (k % 4)
        l = 2 + (k % 16)
        nf = 1 + (k % 5)
        a = rng.normal(size=(n,n))
        radius = max(abs(np.linalg.eigvals(a)))
        a *= (0.15 + 0.8*rng.uniform())/radius
        scale = np.geomspace(0.5, 4.0, num=n)
        b = rng.normal(size=(n,n))*scale[:,None]
        sigma = b @ b.T + np.diag(scale*scale*.4)
        taper = hann(l,sym=False) if k % 3 == 0 else rng.uniform(-1.0,1.0,l)
        freq = rng.uniform(-.5,.5,size=nf)
        if k % 5 == 0:
            freq[0] = 0.0
        scenarios.append((a,sigma,taper,freq))
    lines = [str(cases)]
    for a,sigma,taper,freq in scenarios:
        n = a.shape[0]
        lines.append(f'{n} {len(taper)} {len(freq)}')
        lines.append(' '.join(f'{v:.17g}' for v in a.flatten(order='F')))
        lines.append(' '.join(f'{v:.17g}' for v in sigma.flatten(order='F')))
        lines.append(' '.join(f'{v:.17g}' for v in taper))
        lines.append(' '.join(f'{v:.17g}' for v in freq))
    process = subprocess.run([executable],input='\n'.join(lines)+'\n',
                             text=True,capture_output=True,check=True)
    rows = process.stdout.splitlines()
    if len(rows) != cases:
        raise AssertionError((len(rows),cases,process.stderr))
    max_abs = 0.0
    matrices = 0
    for row,(a,sigma,taper,freq) in zip(rows,scenarios):
        n = a.shape[0]
        values = np.fromstring(row,sep=' ')
        assert values[0] == 0, values[:5]
        expected_count = 1 + n*n + 2*len(freq)*n*n
        assert len(values) == expected_count,(len(values),expected_count)
        gamma_f = values[1:1+n*n].reshape((n,n),order='F')
        real_imag = values[1+n*n:].reshape(len(freq),n*n,2)
        # Each spectral matrix was written in Fortran column order.
        power_f = np.stack([v.reshape((n,n),order='F')
                            for v in real_imag[...,0]+1j*real_imag[...,1]])
        gamma_ref,power_ref = direct_reference(a,sigma,taper,freq)
        abs_diff = max(float(np.max(abs(gamma_f-gamma_ref))),
                       float(np.max(abs(power_f-power_ref))))
        scale = max(1.0,float(np.max(abs(power_ref))),float(np.max(abs(gamma_ref))))
        assert abs_diff < 2.e-10*scale, (abs_diff,scale,a,sigma)
        max_abs = max(max_abs,abs_diff)
        matrices += len(freq)
    print(f'PASS {cases} independent SciPy/NumPy VAR(1) cases, '
          f'{matrices} finite-tapered spectral matrices; '
          f'max abs discrepancy {max_abs:.3e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases',type=int,default=80)
    args = parser.parse_args()
    run(args.executable,args.cases)
