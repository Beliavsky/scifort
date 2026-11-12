#!/usr/bin/env python3
"""Independent NumPy test of leave-one canonical influence and pilot autocovariance.

The Fortran driver only reads CSD matrices and emits diagnostics; all reference
coherence, jackknife, and block-length calculations are implemented here.
"""
from __future__ import annotations

import math
import subprocess
import sys

import numpy as np


def canonical_values(matrix: np.ndarray, conditioned: bool, ridge: float,
                     shrinkage: float) -> np.ndarray:
    """Return [determinant group coherence, root1, root2] independently."""
    nc = matrix.shape[0]
    regular = (matrix + matrix.conj().T) / 2
    regular = (1 - shrinkage) * regular
    regular[np.diag_indices(nc)] = np.diag(matrix).real + ridge * np.trace(matrix).real / nc
    aa = regular[np.ix_([0, 1], [0, 1])]
    bb = regular[np.ix_([2, 3], [2, 3])]
    ab = regular[np.ix_([0, 1], [2, 3])]
    if conditioned:
        ac = regular[np.ix_([0, 1], [4])]
        bc = regular[np.ix_([2, 3], [4])]
        cc = regular[4, 4]
        aa = aa - ac @ ac.conj().T / cc
        bb = bb - bc @ bc.conj().T / cc
        ab = ab - ac @ bc.conj().T / cc
    la = np.linalg.cholesky(aa)
    lb = np.linalg.cholesky(bb)
    whitened = np.linalg.solve(la, ab) @ np.linalg.inv(lb.conj().T)
    roots = np.linalg.svd(whitened, compute_uv=False) ** 2
    aggregate = 1 - np.prod(1 - roots)
    return np.r_[aggregate, roots]


def influence_reference(spectra: np.ndarray, mode: int, frequency: int,
                        conditioned: bool, ridge: float, shrinkage: float) -> tuple[float, np.ndarray]:
    # Spectra shape (window, frequency, channel, channel); recompute every
    # leave-one conditional canonical fit, as opposed to linearizing the CSD.
    n = spectra.shape[0]
    total = spectra.sum(axis=0)
    freq_indices = list(range(spectra.shape[1])) if frequency == 0 else [frequency-1]

    def score(means: np.ndarray) -> float:
        return float(np.mean([canonical_values(means[f], conditioned, ridge, shrinkage)[mode]
                              for f in freq_indices]))

    full = score(total / n)
    influences = np.array([(n-1)*(full-score((total-spectra[t])/(n-1)))
                           for t in range(n)])
    return full, influences


def pilot_reference(scores: np.ndarray, lag_cap: int, length_cap: int) -> tuple[int, int, np.ndarray]:
    n = len(scores)
    centered = scores - scores.mean()
    # The Frobenius pilot is invariant to common normalization. The target
    # score is one-dimensional, so this is exactly its scalar autocovariance.
    gamma0 = np.sum(centered ** 2) / n
    rho = np.zeros(lag_cap)
    if gamma0 == 0:
        return 1, 0, rho
    cutoff = 2 / np.sqrt(n)
    used = 0
    for k in range(1, lag_cap+1):
        rho[k-1] = np.dot(centered[:n-k], centered[k:]) / ((n-k)*gamma0)
        if abs(rho[k-1]) <= cutoff:
            break
        used = k
    if used == 0:
        return 1, 0, rho
    weights = 1 - np.arange(1, used+1)/(used+1)
    longrun = 1 + 2 * np.dot(weights, rho[:used])
    moment = 2 * np.dot(np.arange(1, used+1)*weights, rho[:used])
    if longrun <= np.sqrt(np.finfo(float).eps):
        return 1, used, rho
    proposed = (1.5 * n * (moment/longrun)**2) ** (1 / 3)
    return max(1, min(length_cap, int(math.floor(min(length_cap, proposed)+0.5)))), used, rho


def band_integral(spec: np.ndarray) -> np.ndarray:
    """Integrate linearly interpolated matrices between 0.05 and 0.42 Hz."""
    a, b, mid = 0.05, 0.42, 0.25
    s0, s1, s2 = spec[:, 0], spec[:, 1], spec[:, 2]
    left = s0 * ((mid-a)/mid) + s1 * (a/mid)
    right = s1 * ((0.5-b)/mid) + s2 * ((b-mid)/mid)
    area = 0.5*(mid-a)*(left+s1) + 0.5*(b-mid)*(s1+right)
    return (area/(b-a))[:, None]


def cases():
    rng = np.random.default_rng(20261030)
    for caseid in range(40):
        nw = 32 + caseid % 15
        spec = np.zeros((nw, 3, 5, 5), dtype=np.complex128)
        dependent = 0.0
        for t in range(nw):
            dependent = (0.92 if caseid % 4 == 0 else 0.43)*dependent + rng.normal()
            for f in range(3):
                scale = 1 + 0.04*f
                diagonals = np.array([2.5, 2.2, 2.8, 2.4, 20 + 0.5*dependent]) * scale
                m = np.diag(diagonals.astype(complex))
                c13 = 0.68 + 0.08*np.sin(.15*t+.13*f+.02*caseid) + .06j
                c24 = 0.43 + 0.05*np.cos(.21*t+.25*f) - .03j
                connections = {(0,2): c13, (1,3): c24, (0,4): .15+.03j,
                               (1,4): .06-.04j, (2,4): .11-.02j, (3,4): .05+.08j}
                for (a,b), val in connections.items():
                    m[a,b] = val
                    m[b,a] = np.conj(val)
                spec[t,f] = m
        mode = caseid % 3  # 0: aggregate, 1 or 2: canonical mode
        conditioned = caseid % 2 == 0
        band = caseid % 5 == 0
        frequency = 0 if band else (caseid // 3) % 4
        ridge = 0.02 if caseid % 3 == 0 else 0.0
        shrinkage = 0.12 if caseid % 4 == 1 else 0.0
        lag_cap = min(nw-2, 4 + caseid % 6)
        cap = min(nw-1, 6 + caseid % 12)
        yield spec, mode, frequency, conditioned, band, ridge, shrinkage, lag_cap, cap


def main(driver: str) -> None:
    samples = list(cases())
    lines = [str(len(samples))]
    references = []
    for spec, mode, frequency, cond, band, ridge, shrinkage, lag_cap, cap in samples:
        nw, nf, nc, _ = spec.shape
        lines.append(f'{nw} {nf} {nc} {mode} {frequency} {int(cond)} {int(band)} '
                     f'{lag_cap} {cap} {ridge:.17g} {shrinkage:.17g}')
        lines.append(' '.join(f'{v:.17g}' for v in [0., 0.25, 0.5]))
        lines.append(' '.join(str(i+1) for i in range(nw)))
        for t in range(nw):
            for f in range(nf):
                column_major = spec[t,f].flatten(order='F')
                elements = np.column_stack((column_major.real,column_major.imag)).ravel()
                lines.append(' '.join(f'{v:.17g}' for v in elements))
        target = band_integral(spec) if band else spec
        full, infl = influence_reference(target,mode,frequency,cond,ridge,shrinkage)
        selected, used, acf = pilot_reference(infl,lag_cap,cap)
        references.append((full,infl,selected,used,acf))
    completed = subprocess.run([driver], input='\n'.join(lines)+'\n', text=True,
                               capture_output=True, check=True)
    received = completed.stdout.strip().splitlines()
    assert len(received) == 3*len(samples), (len(received), len(samples), completed.stderr)
    max_err = 0.0
    for i, (full, scores, length, used, acf) in enumerate(references):
        status, got_length, got_used, got_full, first_lag = received[i*3].split()
        actual = np.fromstring(received[i*3+1],sep=' ')
        pilot = np.fromstring(received[i*3+2],sep=' ')
        assert int(status) == 0 and int(got_length) == length and int(got_used) == used, (i,received[i*3],length,used)
        assert actual.shape == scores.shape and pilot.shape == acf.shape
        err = max(np.max(np.abs(actual-scores)),abs(float(got_full)-full),
                  np.max(np.abs(pilot-acf)),abs(float(first_lag)-acf[0]))
        max_err = max(max_err,err)
        assert err < 1e-9, (i,err)
    print(f'PASS {len(samples)} independent complex/conditioned/root/band NumPy '
          f'influence pilots, {sum(len(x[1]) for x in references)} leave-one scores; '
          f'max absolute error {max_err:.3e}')


if __name__ == '__main__':
    main(sys.argv[1])
