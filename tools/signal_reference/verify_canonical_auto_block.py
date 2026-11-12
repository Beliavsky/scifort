#!/usr/bin/env python3
"""Independent NumPy pilot reference and small correlated-AR(1) coverage study.

Inputs to the Fortran driver are complete complex spectral matrices. The
coverage study constructs one-sided window periodograms of a *bivariate*
Gaussian AR(1) with the same AR coefficient in each channel; its population
squared coherence at any frequency is rho**2. Study results are descriptive,
not a claim of calibrated coverage or universally optimal block selection.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.signal.windows import hann


def pilot_reference(spectra: np.ndarray) -> tuple[int, int, int, int, float, float, float]:
    """Compute a separate scalar Frobenius autocovariance plug-in calculation."""
    _, nw, _, _ = spectra.shape
    centered = spectra - np.mean(spectra, axis=1, keepdims=True)
    scale = np.max(np.abs(centered))
    cap = min(nw - 2, max(1, math.isqrt(nw)))
    maximum = max(1, min(nw - 1, nw // 3))
    if scale == 0:
        return 1, 0, cap, maximum, 0.0, 1.0, 0.0
    centered = centered / scale
    variance = np.vdot(centered, centered).real / nw
    if variance <= 0:
        return 1, 0, cap, maximum, 0.0, 1.0, 0.0
    acf = np.empty(cap)
    for lag in range(1, cap + 1):
        lhs = centered[:, :-lag, :, :]
        rhs = centered[:, lag:, :, :]
        acf[lag - 1] = (np.vdot(lhs, rhs).real / (nw - lag)) / variance
    significant = 0
    for lag, corr in enumerate(acf, 1):
        if abs(corr) <= 2.0 / math.sqrt(nw):
            break
        significant = lag
    if not significant:
        return 1, 0, cap, maximum, float(acf[0]), 1.0, 0.0
    lags = np.arange(1, significant + 1)
    weights = 1.0 - lags / (significant + 1.0)
    long_run = float(1.0 + 2.0 * np.dot(weights, acf[:significant]))
    first_moment = float(2.0 * np.dot(weights * lags, acf[:significant]))
    if long_run <= math.sqrt(np.finfo(float).eps):
        return 1, significant, cap, maximum, float(acf[0]), 1.0, 0.0
    ratio = first_moment / long_run
    proposed = (1.5 * nw * ratio**2) ** (1.0 / 3.0)
    length = max(1, min(maximum, int(math.floor(min(maximum, proposed) + .5))))
    return length, significant, cap, maximum, float(acf[0]), long_run, ratio


def synthetic_matrices(rng: np.random.Generator, nw: int, nf: int, nc: int,
                       phi: float) -> np.ndarray:
    """Generate positive-definite Hermitian window spectra with AR(1) variation."""
    z = np.zeros(nw)
    innov = rng.normal(size=nw)
    for t in range(nw):
        z[t] = phi * (z[t - 1] if t else 0.0) + np.sqrt(1 - phi**2) * innov[t]
    z = np.clip(z, -3.0, 3.0)
    result = np.zeros((nf, nw, nc, nc), dtype=np.complex128)
    for f in range(nf):
        a = rng.normal(size=(nc,nc)) + 1j*rng.normal(size=(nc,nc))
        baseline = a @ a.conj().T + 3.0 * np.eye(nc)
        b = rng.normal(size=(nc,nc)) + 1j*rng.normal(size=(nc,nc))
        perturbation = (b + b.conj().T) * 0.12
        result[f] = baseline + z[:,None,None] * perturbation
    return result


def ar1_window_periodograms(rng: np.random.Generator, *, nw: int, nperseg: int,
                              phi: float, rho: float) -> np.ndarray:
    """Spectral matrices at frequency zero from a stationary bivariate AR(1)."""
    n = nw*nperseg
    normal = rng.normal(size=(n+200,2))
    eps = np.column_stack((normal[:,0], rho*normal[:,0] +
                           np.sqrt(1-rho**2)*normal[:,1]))
    x = np.zeros_like(eps)
    for t in range(1, len(x)):
        x[t] = phi*x[t-1] + eps[t]
    x = x[200:]
    taper = hann(nperseg, sym=False)
    fourier = np.sum(x.reshape(nw,nperseg,2) * taper[None,:,None], axis=1)
    return (fourier[:, :, None] * fourier[:, None, :]).astype(complex)[None, :, :, :]


def invoke(executable: Path, cases: list[tuple[int,np.ndarray,int]]) -> np.ndarray:
    text = [str(len(cases))]
    for mode, power, seed in cases:
        nf, nw, nc, _ = power.shape
        text.append(f'{mode} {nf} {nw} {nc} {seed}')
        for f in range(nf):
            for w in range(nw):
                for i in range(nc):
                    for j in range(nc):
                        x = power[f,w,i,j]
                        text.append(f'{x.real:.17g} {x.imag:.17g}')
    completed = subprocess.run([str(executable)], input='\n'.join(text)+'\n',
                               text=True,capture_output=True,check=True)
    results = np.asarray([[float(v) for v in line.split()]
                          for line in completed.stdout.strip().splitlines()])
    if results.shape[0] != len(cases):
        raise AssertionError(f'expected {len(cases)} answers, got {results.shape}')
    return results


def run_reference(executable: Path, cases: int = 85) -> None:
    rng = np.random.default_rng(20261030)
    inputs = []
    samples = []
    phi_choices = [-.90, -.70, -.3, 0.0, .20, .55, .80, .95]
    for rep in range(cases):
        nw = int(rng.integers(18,100))
        nf = int(rng.integers(1,4))
        nc = int(rng.integers(2,5))
        phi = phi_choices[rep % len(phi_choices)]
        spectra = synthetic_matrices(rng,nw,nf,nc,phi)
        if rep == 0:
            spectra[:] = spectra.mean(axis=1,keepdims=True)
        if rep == 1:
            spectra *= 1.e100
        if rep == 2:
            spectra *= 1.e-100
        samples.append(pilot_reference(spectra))
        inputs.append((0,spectra,rep+1))
    actual = invoke(executable,inputs)
    for rep, (expected, observed) in enumerate(zip(samples,actual)):
        np.testing.assert_array_equal(observed[:4],expected[:4],err_msg=f'case {rep}')
        np.testing.assert_allclose(observed[4:],expected[4:],rtol=1e-11,atol=2e-13,
                                   err_msg=f'case {rep}')
    distinct = np.unique(actual[:,0].astype(int))
    print(f'PASS independent NumPy block-selection reference: {cases} Hermitian scenarios; '
          f'selected lengths {distinct.tolist()}',flush=True)


def run_coverage(executable: Path, replications: int = 64) -> None:
    rng = np.random.default_rng(1179)
    rho = .60
    phi = .94
    nw = 48
    nperseg = 12
    cases = [(1,ar1_window_periodograms(rng,nw=nw,nperseg=nperseg,phi=phi,rho=rho),
              1300+r) for r in range(replications)]
    observations = invoke(executable,cases)
    truth = rho**2
    auto = (observations[:,2] <= truth) & (truth <= observations[:,3])
    independent = (observations[:,7] <= truth) & (truth <= observations[:,8])
    if not np.all(observations[:,6] == 1) or not np.all(observations[:,9] == 1):
        raise AssertionError('unexpected invalid bootstrap aggregate')
    if not np.all((observations[:,1] >= 0) & (observations[:,1] <= 1)):
        raise AssertionError('coherence outside [0,1]')
    if not np.all(observations[:,2] <= observations[:,3]):
        raise AssertionError('reversed percentile interval')
    print(f'Exploratory AR(1) study: rho={rho}, phi={phi}, {replications} trials, '
          f'{nw} nonoverlapping {nperseg}-sample windows, 99 draws each',flush=True)
    print(f'  selected blocks: median {np.median(observations[:,0]):.1f}, '
          f'range {int(observations[:,0].min())}..{int(observations[:,0].max())}',flush=True)
    print(f'  empirical 95% percentile coverage: automatic {auto.mean():.3f}, '
          f'length-1 {independent.mean():.3f} (simulation-specific, not calibrated)',flush=True)
    print(f'  average interval widths: automatic '
          f'{np.mean(observations[:,3]-observations[:,2]):.4f}, length-1 '
          f'{np.mean(observations[:,8]-observations[:,7]):.4f}',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver',type=Path)
    parser.add_argument('--reference-cases',type=int,default=85)
    parser.add_argument('--coverage-trials',type=int,default=64)
    args = parser.parse_args()
    run_reference(args.driver.resolve(),args.reference_cases)
    run_coverage(args.driver.resolve(),args.coverage_trials)
