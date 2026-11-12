#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reproducible conditional, two-mode VAR(2)/VAR(3) bootstrap coverage study.

This is a simulation diagnostic, not evidence that bootstrap intervals are
calibrated.  The estimand is canonical(E[window periodogram or integrated
periodogram]), with the same taper and band normalization as the experiment.
The expected nonlinear finite-sample estimator is a distinct quantity.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import subprocess

import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.signal.windows import hann
from verify_varp_taper_spectrum import reference
from verify_varp_canonical_target import canonical

METHODS = ('percentile', 'basic', 'bc', 'bca', 'studentized')
LABELS = ('f=0.1375', 'band [0,0.125]', 'band [0.125,0.5]')
STATISTICS = ('mode 1', 'mode 2', 'group')
FREQUENCIES = np.array([0.1375])
BANDS = np.array([[0.0, 0.125], [0.125, 0.5]])


@dataclass(frozen=True)
class Model:
    name: str
    a: np.ndarray  # (5,5,p), first lag first
    sigma: np.ndarray


def _innovation_covariance(weak: bool = False, tied: bool = False) -> np.ndarray:
    """Construct a positive-definite covariance from independent latent shocks."""
    mix = np.eye(5)
    mix[0, 4] = 0.35
    mix[1, 4] = 0.20
    mix[2, 4] = 0.35
    mix[3, 4] = 0.20
    mix[2, 0] = 0.72
    mix[2, 2] = 0.78
    mix[3, 1] = (0.16 if weak else (0.72 if tied else 0.43))
    mix[3, 3] = 0.90 if weak else (0.78 if tied else 0.84)
    return mix @ mix.T


def _diagonal_var(*diagonals: list[float]) -> np.ndarray:
    return np.stack([np.diag(d) for d in diagonals], axis=-1)


def models() -> tuple[Model, ...]:
    """Stable models with two associations, one nuisance channel, and lag coupling."""
    moderate = _diagonal_var([.60, .45, .58, .46, .75],
                             [-.12, .06, -.10, .04, -.06])
    moderate[0, 4, 0] = .10
    moderate[2, 4, 1] = -.07
    lagged = _diagonal_var([.65, .57, .68, .47, .78],
                           [-.18, -.08, -.21, -.06, -.09],
                           [.07, .04, .08, .03, .02])
    lagged[0, 2, 0] = .10
    lagged[1, 3, 1] = -.14
    lagged[2, 4, 2] = .08
    weak = _diagonal_var([.76, .59, .73, .57, .64],
                         [-.17, -.06, -.11, -.07, .04])
    weak[0, 2, 0] = .07
    tied = _diagonal_var([.49, .49, .49, .49, .67],
                         [-.04, -.04, -.04, -.04, -.02],
                         [.02, .02, .02, .02, .01])
    return (Model('VAR(2) separated modes', moderate, _innovation_covariance()),
            Model('VAR(3) cross-lag coupling', lagged, _innovation_covariance()),
            Model('VAR(2) weak second mode', weak, _innovation_covariance(weak=True)),
            Model('VAR(3) near-tied modes', tied, _innovation_covariance(tied=True)))


def companion_and_covariance(model: Model) -> tuple[np.ndarray, np.ndarray]:
    """Solve stationary companion-state covariance, without a burn-in approximation."""
    n, _, p = model.a.shape
    transition = np.zeros((n*p, n*p))
    transition[:n] = np.concatenate([model.a[:, :, lag] for lag in range(p)], axis=1)
    if p > 1:
        transition[n:, :-n] = np.eye(n*(p-1))
    if max(abs(np.linalg.eigvals(transition))) >= 1.0:
        raise ValueError('VAR companion matrix is not stationary')
    noise = np.zeros_like(transition)
    noise[:n, :n] = model.sigma
    state_covariance = solve_discrete_lyapunov(transition, noise)
    if not np.allclose(state_covariance, transition @ state_covariance @
                       transition.T + noise, atol=1e-10, rtol=1e-10):
        raise AssertionError('companion Lyapunov residual')
    return transition, state_covariance


def simulate_windows(rng: np.random.Generator, model: Model, nw: int,
                     length: int) -> np.ndarray:
    """Simulate exactly stationary Gaussian VAR(p) windows (no warmup)."""
    transition, covariance = companion_and_covariance(model)
    n = model.a.shape[0]
    state = np.linalg.cholesky(covariance) @ rng.standard_normal(len(transition))
    innovation_factor = np.linalg.cholesky(model.sigma)
    data = np.empty((nw*length, n))
    for t in range(len(data)):
        state = transition @ state
        state[:n] += innovation_factor @ rng.standard_normal(n)
        data[t] = state[:n]
    return data.reshape(nw, length, n)


def window_spectra(windows: np.ndarray, taper: np.ndarray) -> np.ndarray:
    """Exact per-window Fourier periodogram and analytic band-averaged spectra."""
    nw, length, n = windows.shape
    if taper.shape != (length,):
        raise ValueError('taper shape mismatch')
    energy = np.dot(taper, taper)
    if energy <= 0:
        raise ValueError('zero taper')
    samples = windows * taper[None, :, None]
    t = np.arange(length)
    phase = np.exp(-2j*np.pi*FREQUENCIES[0]*t)
    fourier = np.einsum('wti,t->wi', samples, phase)
    point = np.einsum('wi,wj->wij', fourier, fourier.conj()) / energy
    bins = [point]
    lag = t[:, None]-t[None, :]
    for lower, upper in BANDS:
        width, midpoint = upper-lower, (lower+upper)/2
        kernel = np.sinc(width*lag) * np.exp(-2j*np.pi*midpoint*lag)
        average = np.einsum('wti,ts,wsj->wij', samples, kernel, samples,
                            optimize=True) / energy
        bins.append((average+average.swapaxes(1,2).conj())/2)
    output = np.stack(bins)
    if output.shape != (1+len(BANDS), nw, n, n):
        raise AssertionError('unexpected spectrum dimensions')
    return output


def population_targets(model: Model, taper: np.ndarray) -> np.ndarray:
    """Independent SciPy Lyapunov, time-pair sums and SVD canonical targets."""
    _, frequencies, integrated = reference(model.a, model.sigma, taper,
                                            FREQUENCIES, BANDS,
                                            normalized=True, average=True)
    return canonical(np.concatenate((frequencies, integrated)))[:, [1,2,0]]


def wilson(successes: int, count: int) -> tuple[float, float]:
    """Two-sided 95% Monte Carlo uncertainty, not a bootstrap confidence interval."""
    if not count:
        return float('nan'), float('nan')
    z = 1.959963984540054
    p = successes/count
    divisor = 1+z*z/count
    center = (p+z*z/(2*count))/divisor
    half = z*np.sqrt(p*(1-p)/count+z*z/(4*count*count))/divisor
    return float(center-half), float(center+half)


def encode_dataset(spectra: np.ndarray, seed: int, block_length: int,
                   resamples: int) -> list[str]:
    """Send the full complex Hermitian cross-spectrum, Fortran column order."""
    nf, nw, nc, _ = spectra.shape
    lines = [f'{nw} {nc} {nf} {block_length} {resamples} {seed}',
             " ".join(f'{f:.17g}' for f in [FREQUENCIES[0],
                          *((BANDS[:,0]+BANDS[:,1])/2)])]
    for t in range(nw):
        for f in range(nf):
            values = spectra[f, t].flatten(order='F')
            joined = np.column_stack((values.real, values.imag)).ravel()
            lines.append(' '.join(f'{x:.17g}' for x in joined))
    return lines


def parse_output(stdout: str, expected: int, nf: int = 3,
                 nm: int = 2) -> np.ndarray:
    """Validate every reported interval and its validity flag."""
    records = [np.fromstring(line, sep=' ') for line in stdout.splitlines()]
    ncols = nf*(nm+1)*(1+3*len(METHODS))
    if len(records) != expected or any(len(record) != ncols for record in records):
        raise AssertionError(f'expected {expected} rows x {ncols}, got '
                             f'{len(records)} rows and {[len(r) for r in records[:3]]}')
    data = np.stack(records).reshape((expected, nf, nm+1, 1+3*len(METHODS)))
    point = data[..., 0]
    if not np.all(np.isfinite(point) & (point >= -1e-10) & (point <= 1+1e-10)):
        raise AssertionError('invalid full-sample canonical estimate')
    for m in range(len(METHODS)):
        lower, upper, valid = (data[...,1+3*m+i] for i in range(3))
        if not np.all(np.isin(valid, (0, 1))):
            raise AssertionError('invalid Boolean result')
        mask = valid == 1
        if not np.all(np.isfinite(lower[mask]) & np.isfinite(upper[mask])):
            raise AssertionError('nonfinite valid interval')
        if not np.all((lower[mask] >= -1e-10) & (upper[mask] <= 1+1e-10) &
                      (lower[mask] <= upper[mask])):
            raise AssertionError('interval outside parameter support')
    return data


def conduct_study(executable: str | Path, *, trials: int = 8, nw: int = 40,
                  length: int = 16, resamples: int = 79, block_length: int = 4,
                  seed: int = 20261030, report: bool = True,
                  csv_path: str | Path | None = None) -> dict:
    """Run the Fortran bootstrap and summarize valid-case coverage by scalar estimand."""
    if trials < 1 or nw < 8 or length < 3 or resamples < 20 or not 1 <= block_length < nw:
        raise ValueError('invalid simulation sizes or resampling configuration')
    rng = np.random.default_rng(seed)
    taper = hann(length, sym=False)
    all_models = models()
    targets = [population_targets(model, taper) for model in all_models]
    spec_list = []
    lines = [str(trials*len(all_models))]
    for k, model in enumerate(all_models):
        for trial in range(trials):
            spectra = window_spectra(simulate_windows(rng, model, nw, length), taper)
            spec_list.append(spectra)
            lines.extend(encode_dataset(spectra, 1000 + k*trials + trial,
                                        block_length, resamples))
    process = subprocess.run([str(executable)], input='\n'.join(lines)+'\n',
                             capture_output=True, text=True, check=True)
    observations = parse_output(process.stdout, trials*len(all_models))
    # Independently refit the full mean cross-spectrum for every dataset.
    independent = np.stack([canonical(np.mean(spectrum, axis=1))[:, [1,2,0]]
                            for spectrum in spec_list])
    max_error = float(np.max(abs(observations[...,0]-independent)))
    if max_error > 5e-9:
        raise AssertionError(f'Fortran point canonical modes disagree with SVD: {max_error}')
    summaries = {}
    csv_rows = []
    if report:
        print(f'VAR(p) conditional canonical coherence study: {len(all_models)} models x '
              f'{trials} datasets, {nw} stationary nonoverlapping '
              f'{length}-sample Hann windows, {resamples} circular block draws '
              f'of length {block_length}')
        print(f'Independent point-estimate SVD comparison: max error {max_error:.3e}')
        print('Coverage is conditional on valid intervals; Wilson ranges quantify '
              'only finite Monte Carlo simulation uncertainty.')
    for k, model in enumerate(all_models):
        obs = observations[k*trials:(k+1)*trials]
        target = targets[k]
        rows = []
        if report:
            print(f'\n{model.name}; target rows are frequency, low band, high band:')
            for label, items in zip(LABELS, target):
                print(f'  {label:19s} mode1 {items[0]:.4f}, mode2 {items[1]:.4f}, '
                      f'aggregate {items[2]:.4f}')
        for f, name in enumerate(LABELS):
            for s, statistic in enumerate(STATISTICS):
                point_bias = float(np.mean(obs[:,f,s,0]-target[f,s]))
                point_rmse = float(np.sqrt(np.mean((obs[:,f,s,0]-target[f,s])**2)))
                for m, method in enumerate(METHODS):
                    lower, upper, flags = (obs[:,f,s,1+3*m+i] for i in range(3))
                    valid = (flags == 1) & np.isfinite(lower) & np.isfinite(upper)
                    count = int(valid.sum())
                    hits = int(np.sum((lower[valid] <= target[f,s]) &
                                      (target[f,s] <= upper[valid])))
                    ci = wilson(hits, count)
                    avg_width = float(np.mean(upper[valid]-lower[valid])) if count else float('nan')
                    rows.append((name, statistic, method, hits, count, ci, avg_width))
                    csv_rows.append(dict(model=model.name, frequency=name,
                        statistic=statistic, method=method,
                        population_target=float(target[f,s]),
                        population_root_gap=float(target[f,0]-target[f,1]),
                        datasets=trials, valid_intervals=count, covered=hits,
                        coverage=(hits/count if count else float('nan')),
                        wilson_lower=ci[0], wilson_upper=ci[1],
                        mean_interval_width=avg_width,
                        full_estimator_bias=point_bias,
                        full_estimator_rmse=point_rmse))
                    if report:
                        display = f'{100*hits/count:5.1f}%' if count else '  n/a'
                        print(f'  {name:19s} {statistic:7s} {method:11s} '
                              f'coverage {display}, valid {count:3d}/{trials}, '
                              f'Wilson 95% {100*ci[0]:5.1f}..{100*ci[1]:5.1f}%, '
                              f'mean width {avg_width:.4f}')
        summaries[model.name] = (target, rows)
    if csv_path is not None:
        destination = Path(csv_path)
        with destination.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
            writer.writeheader()
            writer.writerows(csv_rows)
        if report:
            print(f'\nWrote {len(csv_rows)} coverage/statistical-bias rows to {destination}')
    return {'targets': targets, 'summaries': summaries, 'data': observations,
            'max_point_error': max_error, 'csv_rows': csv_rows}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver', type=Path)
    parser.add_argument('--trials', type=int, default=8)
    parser.add_argument('--windows', type=int, default=40)
    parser.add_argument('--samples-per-window', type=int, default=16)
    parser.add_argument('--resamples', type=int, default=79)
    parser.add_argument('--block-length', type=int, default=4)
    parser.add_argument('--seed', type=int, default=20261030)
    parser.add_argument('--csv', type=Path, help='save 180 summary rows as CSV')
    args = parser.parse_args()
    conduct_study(args.driver, trials=args.trials, nw=args.windows,
                  length=args.samples_per_window, resamples=args.resamples,
                  block_length=args.block_length, seed=args.seed,
                  csv_path=args.csv)
