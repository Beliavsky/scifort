#!/usr/bin/env python3
"""Independent NumPy/SciPy jackknife of conditional canonical spectral modes.

The reference uses NumPy Cholesky/SVD, SciPy interpolation/integration, and
independent explicit leave-one-window calculations, not the SciFort solver.
"""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid


@dataclass
class Case:
    frequencies: np.ndarray
    samples: np.ndarray  # frequency, window, channel, channel
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    bands: np.ndarray | None
    ridge: float
    shrink: float
    tol: float = 1e-6


def canonical(csd: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray,
              ridge: float, shrink: float) -> tuple[np.ndarray, float]:
    csd = (csd + csd.conj().T) / 2
    s = (1 - shrink) * csd
    np.fill_diagonal(s, np.real(np.diag(csd)) + ridge * np.real(np.diag(csd)).mean())
    def block(i, j):
        return s[np.ix_(i, j)]
    ca = block(a, a)
    cb = block(b, b)
    cab = block(a, b)
    if len(c):
        cc = block(c, c)
        ca -= block(a, c) @ np.linalg.solve(cc, block(c, a))
        cb -= block(b, c) @ np.linalg.solve(cc, block(c, b))
        cab -= block(a, c) @ np.linalg.solve(cc, block(c, b))
    la = np.linalg.cholesky(ca)
    lb = np.linalg.cholesky(cb)
    white = np.linalg.solve(la, cab) @ np.linalg.inv(lb.conj().T)
    roots = np.clip(np.linalg.svd(white, compute_uv=False)**2, 0, 1)
    group = 1.0 - np.prod(1.0 - roots)
    return roots, group


def integrated(s: Case) -> np.ndarray:
    if s.bands is None:
        return s.samples
    output = []
    # Frequency interpolation is performed on each complex matrix component;
    # trapezoid of the piecewise-linear interpolant integrates exact endpoints.
    for lo, hi in s.bands:
        interior = s.frequencies[(s.frequencies > lo) & (s.frequencies < hi)]
        grid = np.r_[lo, interior, hi]
        interpolator = interp1d(s.frequencies, s.samples, axis=0, kind='linear')
        output.append(trapezoid(interpolator(grid), x=grid, axis=0))
    return np.stack(output, axis=0)


def regular(values: np.ndarray, k: int, tol: float) -> bool:
    if values[k] <= tol or values[k] >= 1 - tol:
        return False
    if k and values[k - 1] - values[k] <= tol:
        return False
    if k + 1 < len(values) and values[k] - values[k + 1] <= tol:
        return False
    return True


def reference(s: Case):
    x = integrated(s)
    nf, nt = x.shape[:2]
    nm = min(len(s.a), len(s.b))
    estimates = np.zeros((nf, nm))
    se = estimates.copy()
    se_z = estimates.copy()
    lows = estimates.copy()
    highs = estimates.copy()
    bias = estimates.copy()
    valid = np.ones((nf, nm), dtype=bool)
    group = np.zeros(nf)
    se_group = np.zeros(nf)
    for f in range(nf):
        estimates[f], group[f] = canonical(x[f].mean(axis=0), s.a, s.b, s.c,
                                            s.ridge, s.shrink)
        loo_roots = np.zeros((nt, nm))
        loo_groups = np.zeros(nt)
        for i in range(nt):
            loo_roots[i], loo_groups[i] = canonical(
                (x[f].sum(axis=0) - x[f, i]) / (nt - 1),
                s.a, s.b, s.c, s.ridge, s.shrink)
        se_group[f] = np.sqrt((nt - 1) / nt * np.square(
            loo_groups - loo_groups.mean()).sum())
        for k in range(nm):
            valid[f, k] = (regular(estimates[f], k, s.tol)
                           and all(regular(loo_roots[i], k, s.tol) for i in range(nt)))
            if not valid[f, k]:
                continue
            roots = loo_roots[:, k]
            se[f, k] = np.sqrt((nt - 1) / nt * np.square(roots - roots.mean()).sum())
            bias[f, k] = (nt - 1) * (roots.mean() - estimates[f, k])
            zz = np.arctanh(np.sqrt(np.clip(roots, 0, 1 - 1e-12)))
            se_z[f, k] = np.sqrt((nt - 1) / nt * np.square(zz - zz.mean()).sum())
            center = np.arctanh(np.sqrt(min(1 - 1e-12, estimates[f, k])))
            lows[f, k] = np.tanh(max(0.0, center - 1.959963984540054 * se_z[f, k]))**2
            highs[f, k] = np.tanh(center + 1.959963984540054 * se_z[f, k])**2
    return estimates, se, se_z, lows, highs, bias, valid, group, se_group


def cases() -> list[Case]:
    rng = np.random.default_rng(20311030)
    all_cases = []
    for j in range(100):
        nc = int(rng.integers(3, 7))
        na = int(rng.integers(1, nc))
        nb = int(rng.integers(1, nc - na + 1))
        c = np.arange(na + nb, nc, dtype=int)
        nf = int(rng.integers(2, 6))
        nt = int(rng.integers(5, 11))
        ff = np.cumsum(rng.uniform(0.12, 0.35, size=nf))
        xx = rng.normal(size=(nf, nt, nc, nc)) + 1j * rng.normal(size=(nf, nt, nc, nc))
        # Give the windows a common stationary covariance component, with
        # window-specific variation, all strictly positive definite.
        shared = rng.normal(size=(nc, nc)) + 1j * rng.normal(size=(nc, nc))
        shared = shared @ shared.conj().T + 0.8 * np.eye(nc)
        sample = np.empty((nf, nt, nc, nc), dtype=np.complex128)
        for f in range(nf):
            for t in range(nt):
                cov = xx[f, t] @ xx[f, t].conj().T / nc
                sample[f, t] = (0.8 * shared + 0.2 * cov) * (1 + 0.07 * f)
        bands = None
        if j % 3 == 0:
            bands = np.array([[0.7 * ff[0] + 0.3 * ff[1],
                               0.6 * ff[-2] + 0.4 * ff[-1]]]) if nf > 2 else np.array(
                                   [[ff[0] + .05 * (ff[1] - ff[0]),
                                     ff[0] + .9 * (ff[1] - ff[0])]])
            if nf > 3:
                bands = np.vstack([bands, [ff[1] + .1 * (ff[2] - ff[1]),
                                                 ff[-2] + .8 * (ff[-1] - ff[-2])]])
        all_cases.append(Case(ff, sample, np.arange(na), np.arange(na, na + nb),
                              c, bands, float(rng.uniform(0.005, 0.07)),
                              float(rng.uniform(0.03, 0.4))))
    return all_cases


def stdin_payload(test_cases: list[Case]) -> str:
    lines = [str(len(test_cases))]
    for s in test_cases:
        nf, nt, nc = s.samples.shape[:3]
        nbands = 0 if s.bands is None else len(s.bands)
        lines.append(f'{nf} {nt} {nc} {len(s.a)} {len(s.b)} {len(s.c)} {nbands} '
                     f'{s.ridge:.17g} {s.shrink:.17g} {s.tol:.17g}')
        for arr in (s.a + 1, s.b + 1):
            lines.append(' '.join(map(str, arr)))
        if len(s.c):
            lines.append(' '.join(map(str, s.c + 1)))
        if nbands:
            lines.append(' '.join(map(str, s.bands.ravel())))
        lines.append(' '.join(f'{x:.17g}' for x in s.frequencies))
        lines.append(' '.join(map(str, np.arange(nt, dtype=float) + 1)))
        for f in range(nf):
            for t in range(nt):
                for row in s.samples[f, t]:
                    lines.append(' '.join(f'({v.real:.17e},{v.imag:.17e})' for v in row))
    return '\n'.join(lines) + '\n'


def main(binary: str):
    test_cases = cases()
    run = subprocess.run([binary], input=stdin_payload(test_cases), text=True,
                         capture_output=True, check=True)
    lines = iter(run.stdout.splitlines())
    counts = 0
    maximum_error = 0.0
    for i, s in enumerate(test_cases):
        nf, nm, nw = [int(x) for x in next(lines).split()]
        ref = reference(s)
        assert nf == ref[0].shape[0] and nm == ref[0].shape[1] and nw == len(s.samples[0])
        for f in range(nf):
            row = next(lines).split()
            assert len(row) == 4 and int(row[0]) == f + 1 and int(row[3]) == 1
            maximum_error = max(maximum_error, abs(float(row[1]) - ref[7][f]),
                                abs(float(row[2]) - ref[8][f]))
            for k in range(nm):
                row = next(lines).split()
                assert int(row[0]) == f + 1 and int(row[1]) == k + 1
                got_valid = int(row[-1]) != 0
                assert got_valid == ref[6][f, k], (i, f, k, got_valid, ref[6][f, k])
                expected = [r[f, k] for r in ref[:6]]
                got = [float(value) for value in row[2:-1]]
                if got_valid:
                    maximum_error = max(maximum_error, float(np.max(np.abs(
                        np.asarray(got) - np.asarray(expected)))))
                    counts += 1
    assert maximum_error < 1e-9, maximum_error
    print(f'PASS {len(test_cases)} independently verified NumPy/SciPy cases, '
          f'{counts} canonical modes, max error {maximum_error:.3g}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('driver')
    args = parser.parse_args()
    main(args.driver)
