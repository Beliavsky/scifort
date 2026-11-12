#!/usr/bin/env python3
"""Independent NumPy metric-weighted principal angles and Procrustes reference."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
import subprocess

import numpy as np
from scipy.optimize import linear_sum_assignment

from verify_canonical_subspace_tracking import (
    blocks, cases as legacy_cases, comparison, unitary,
)


@dataclass
class MetricCase:
    freq: np.ndarray
    values: np.ndarray
    a: np.ndarray
    b: np.ndarray
    valid: np.ndarray
    gap: float
    minimum: float
    sa: np.ndarray
    sb: np.ndarray


def metric_compare(pa, pb, ca, cb, psa, psb, csa, csb):
    la = np.linalg.cholesky(0.5 * (psa + csa))
    lb = np.linalg.cholesky(0.5 * (psb + csb))
    return comparison(la.conj().T @ pa, lb.conj().T @ pb,
                      la.conj().T @ ca, lb.conj().T @ cb)


def reference(s: MetricCase):
    nf, nm = s.values.shape
    a, b = s.a.copy(), s.b.copy()
    ids = np.zeros((nf, nm), dtype=int)
    first, size = ids.copy(), ids.copy()
    matched, new = ids.astype(bool), ids.astype(bool)
    cosines = np.full((nf, nm), np.nan)
    similarity = np.full((nf, nm), np.nan)
    prev = []
    next_id = 0
    for f in range(nf):
        if not s.valid[f]:
            prev = []
            continue
        curr = blocks(s.values[f], s.gap)
        curr_good = {}
        for j, (lo, hi) in enumerate(curr):
            first[f, lo:hi] = lo + 1
            size[f, lo:hi] = hi - lo
            curr_good[j] = (np.all(np.isfinite(s.values[f, lo:hi]))
                            and s.values[f, lo:hi].min() > 256 * np.finfo(float).eps
                            and np.all(np.isfinite(a[f, :, lo:hi]))
                            and np.all(np.isfinite(b[f, :, lo:hi])))
        links = {}
        if prev and np.isfinite(s.freq[f]) and np.isfinite(s.freq[f - 1]) \
                and s.freq[f] > s.freq[f - 1]:
            nmax = max(len(prev), len(curr))
            scores = np.zeros((nmax, nmax))
            for i, (plo, phi, good) in enumerate(prev):
                if not good:
                    continue
                for j, (clo, chi) in enumerate(curr):
                    if not curr_good[j] or (phi - plo) != (chi - clo):
                        continue
                    candidate = metric_compare(a[f - 1, :, plo:phi], b[f - 1, :, plo:phi],
                                               a[f, :, clo:chi], b[f, :, clo:chi],
                                               s.sa[f - 1], s.sb[f - 1], s.sa[f], s.sb[f])
                    if candidate is None:
                        continue
                    vmin = candidate[0].min()
                    if vmin > 64 * np.finfo(float).eps and vmin >= s.minimum:
                        scores[i, j] = 1 + vmin / (nmax + 1)
            ii, jj = linear_sum_assignment(-scores)
            for i, j in zip(ii, jj):
                if i < len(prev) and j < len(curr) and scores[i, j] > 0:
                    links[j] = i
        for j, (lo, hi) in enumerate(curr):
            if j in links:
                plo, phi, _ = prev[links[j]]
                val, rot = metric_compare(a[f - 1, :, plo:phi], b[f - 1, :, plo:phi],
                                          a[f, :, lo:hi], b[f, :, lo:hi],
                                          s.sa[f - 1], s.sb[f - 1], s.sa[f], s.sb[f])
                ids[f, lo:hi] = ids[f - 1, plo]
                aa = a[f, :, lo:hi].copy()
                bb = b[f, :, lo:hi].copy()
                a[f, :, lo:hi] = aa @ rot
                b[f, :, lo:hi] = bb @ rot
                cosines[f, lo:hi] = val
                similarity[f, lo:hi] = val[-1]
                matched[f, lo:hi] = True
            else:
                next_id += 1
                ids[f, lo:hi] = next_id
                new[f, lo:hi] = True
        prev = [(lo, hi, curr_good[j]) for j, (lo, hi) in enumerate(curr)]
    return ids, first, size, matched, new, cosines, similarity, a, b


def cases(seed: int = 20261030, count: int = 80) -> list[MetricCase]:
    rng = np.random.default_rng(seed + 19)
    scenarios = []
    for k, s in enumerate(legacy_cases(seed=seed, ncases=count)):
        nf, na, nm = s.a.shape
        nb = s.b.shape[1]
        sa = np.zeros((nf, na, na), dtype=np.complex128)
        sb = np.zeros((nf, nb, nb), dtype=np.complex128)
        a, b = s.a.copy(), s.b.copy()
        for f in range(nf):
            for metric, coefficients, n in ((sa, a, na), (sb, b, nb)):
                u = unitary(rng, n)
                spectrum = (u * np.exp(rng.uniform(-0.9, 0.9, n))) @ u.conj().T
                metric[f] = spectrum
                # Canonical coefficients are normalized in their own metric.
                chol = np.linalg.cholesky(spectrum)
                coefficients[f] = np.linalg.solve(chol.conj().T, coefficients[f])
        if k % 4 == 0:
            # Exact change of coordinates, common to all bins in one case.
            for metric, coeff, n in ((sa, a, na), (sb, b, nb)):
                u = unitary(rng, n)
                scale = np.exp(rng.uniform(-2, 2, n))
                transform = u * scale[None, :]
                for f in range(nf):
                    metric[f] = transform @ metric[f] @ transform.conj().T
                    coeff[f] = np.linalg.solve(transform.conj().T, coeff[f])
        scenarios.append(MetricCase(s.freq, s.values, a, b, s.valid,
                                    s.gap, s.minimum, sa, sb))
    return scenarios


def check_congruence_invariance(seed: int = 20261031, count: int = 120) -> None:
    """Non-diagonal complex channel mixing must not change overlap statistics."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(count):
        dimensions = 2 + int(rng.integers(0, 3))
        n = dimensions + 1
        covs = []
        vectors = []
        transforms = []
        for _group in range(2):
            w = []
            for _frequency in range(2):
                mixing = unitary(rng, n)
                cov = (mixing * np.exp(rng.uniform(-0.8, 0.8, n))) @ mixing.conj().T
                w.append(cov)
            covs.append(w)
            v = []
            for f in range(2):
                raw = rng.normal(size=(n, dimensions)) + 1j*rng.normal(size=(n, dimensions))
                raw, _ = np.linalg.qr(raw)
                v.append(np.linalg.solve(np.linalg.cholesky(w[f]).conj().T, raw))
            vectors.append(v)
            transform = (unitary(rng, n) * np.exp(rng.uniform(-1.1, 1.1, n))) \
                @ unitary(rng, n)
            transforms.append(transform)
        ref = metric_compare(vectors[0][0], vectors[1][0],
                             vectors[0][1], vectors[1][1],
                             covs[0][0], covs[1][0], covs[0][1], covs[1][1])
        assert ref is not None
        moved_c = []
        moved_v = []
        for group in range(2):
            transform = transforms[group]
            moved_c.append([transform @ w @ transform.conj().T for w in covs[group]])
            moved_v.append([np.linalg.solve(transform.conj().T, v)
                            for v in vectors[group]])
        moved = metric_compare(moved_v[0][0], moved_v[1][0],
                               moved_v[0][1], moved_v[1][1],
                               moved_c[0][0], moved_c[1][0],
                               moved_c[0][1], moved_c[1][1])
        assert moved is not None
        err = max(np.max(np.abs(ref[0]-moved[0])),
                  np.max(np.abs(ref[1]-moved[1])))
        worst = max(worst, err)
        assert err < 3.e-11, ('congruence invariance', err)
    print(f'PASS {count} complex non-diagonal congruence invariance cases, '
          f'max error {worst:.3g}', flush=True)


def protocol(scenarios: list[MetricCase]) -> str:
    chunks = [str(len(scenarios))]
    for s in scenarios:
        nf, nm = s.values.shape
        na, nb = s.a.shape[1], s.b.shape[1]
        chunks.append(f'{nf} {nm} {na} {nb} {s.gap:.17g} {s.minimum:.17g}')
        chunks.append(' '.join(f'{x:.17g}' for x in s.freq))
        chunks.append(' '.join(str(int(x)) for x in s.valid))
        for f in range(nf):
            chunks.append(' '.join(f'{x:.17g}' for x in s.values[f]))
            for k in range(nm):
                for matrix in (s.a, s.b):
                    for value in matrix[f, :, k]:
                        chunks.append(f'({value.real:.17g},{value.imag:.17g})')
            for matrix in (s.sa, s.sb):
                for k in range(matrix.shape[1]):
                    for value in matrix[f, :, k]:
                        chunks.append(f'({value.real:.17g},{value.imag:.17g})')
    return '\n'.join(chunks) + '\n'


def run(driver: str, scenarios: list[MetricCase]) -> None:
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
    proc = subprocess.run([driver], input=protocol(scenarios), text=True,
                          check=True, capture_output=True, env=env)
    lines = iter(proc.stdout.splitlines())
    worst = 0.0
    observed = 0
    linked = 0
    for case, s in enumerate(scenarios):
        nf, nm = map(int, next(lines).split())
        assert (nf, nm) == s.values.shape
        truth = reference(s)
        na, nb = s.a.shape[1], s.b.shape[1]
        for f in range(nf):
            for k in range(nm):
                fields = next(lines).split()
                integers = np.asarray(fields[:7], dtype=int)
                assert (integers[0], integers[1]) == (f + 1, k + 1)
                expected = [truth[0][f, k], truth[1][f, k], truth[2][f, k],
                            int(truth[3][f, k]), int(truth[4][f, k])]
                if integers[2:].tolist() != expected:
                    raise AssertionError(('assignment', case, f, k, integers[2:].tolist(), expected))
                vals = np.asarray(fields[7:], dtype=float)
                for got, wanted in zip(vals[:2], (truth[5][f, k], truth[6][f, k])):
                    if np.isnan(wanted):
                        assert np.isnan(got), (case, f, k, got, wanted)
                    else:
                        err = abs(got - wanted)
                        worst = max(worst, err)
                        assert err < 3.e-10, ('cosines', case, f, k, err)
                actual = vals[2::2] + 1j * vals[3::2]
                for got, wanted in ((actual[:na], truth[7][f, :, k]),
                                    (actual[na:], truth[8][f, :, k])):
                    err = np.max(np.abs(got - wanted))
                    worst = max(worst, err)
                    assert err < 2.e-9, ('Procrustes', case, f, k, err)
                observed += 1
                linked += bool(truth[3][f, k])
    assert next(lines, None) is None
    assert linked > 100
    print(f'PASS {len(scenarios)} metric-aware complex scenarios, {observed} modes, '
          f'{linked} linked, worst abs error {worst:.3g}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    args = parser.parse_args()
    check_congruence_invariance()
    run(args.driver, cases())
