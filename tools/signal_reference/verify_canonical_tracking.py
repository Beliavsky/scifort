#!/usr/bin/env python3
"""Independent randomized SciPy Hungarian reference for canonical-mode tracking.

This deliberately constructs valid coupled A/B singular directions as orthonormal
unitary columns, changes their phases and crossing strengths, and validates the
Fortran assignment and continuity gauge against scipy.optimize and NumPy.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess

import numpy as np
from scipy.linalg import expm
from scipy.optimize import linear_sum_assignment


def random_unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    z = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    q, r = np.linalg.qr(z)
    return q * np.conj(np.diag(r) / np.abs(np.diag(r)))


def generate(rng: np.random.Generator, case_id: int) -> dict:
    na = int(rng.integers(2, 7))
    nb = int(rng.integers(2, 7))
    nm = min(na, nb)
    nf = int(rng.integers(6, 13))
    f = np.arange(nf, dtype=float)
    frequencies = np.copy(f)
    if case_id % 4 == 0:
        frequencies[nf // 2 :] -= 2.0 * nf
    valid = np.ones(nf, dtype=bool)
    if case_id % 6 == 0:
        valid[nf // 2] = False
    stable = np.ones((nf, nm), dtype=bool)
    if case_id % 5 == 0:
        stable[nf // 2, 0] = False
    minimum = [0.4, 0.6, 0.85, 0.95][case_id % 4]
    qa = random_unitary(rng, na)
    qb = random_unitary(rng, nb)
    za = rng.standard_normal((na, na)) + 1j * rng.standard_normal((na, na))
    zb = rng.standard_normal((nb, nb)) + 1j * rng.standard_normal((nb, nb))
    ha = 0.08 * (za - za.conj().T) / np.sqrt(na)
    hb = 0.08 * (zb - zb.conj().T) / np.sqrt(nb)
    a = np.zeros((nf, na, nm), dtype=complex)
    b = np.zeros((nf, nb, nm), dtype=complex)
    values = np.zeros((nf, nm))
    for j in range(nf):
        physical = np.array([
            0.22 + 0.28 * (1 + np.sin(0.67 * j + 1.12 * k))
            for k in range(nm)
        ])
        order = np.argsort(-physical, kind='stable')
        phase = np.exp(1j * (0.27 * j + 0.42 * np.arange(nm) + 0.12 * np.sin(j)))
        a[j] = (expm(j * ha) @ qa)[:, order] * phase[order]
        b[j] = (expm(j * hb) @ qb)[:, order] * phase[order]
        values[j] = physical[order] ** 2
    return dict(nf=nf, na=na, nb=nb, nm=nm, minimum=minimum,
                frequencies=frequencies, valid=valid, stable=stable,
                a=a, b=b, values=values)


def reference(case: dict) -> dict:
    nf, nm = case['nf'], case['nm']
    a, b = case['a'], case['b']
    track_a = a.copy()
    track_b = b.copy()
    modes = case['values'].copy()
    source = np.zeros((nf, nm), dtype=int)
    overlap = np.full((nf, nm), np.nan)
    matched = np.zeros((nf, nm), dtype=bool)
    active = False
    for f in range(nf):
        if not case['valid'][f]:
            active = False
            continue
        source[f] = np.arange(1, nm + 1)
        ready = case['stable'][f].all()
        if ready and active and np.isfinite(case['frequencies'][f:f + 1]).all():
            if not (np.isfinite(case['frequencies'][f - 1]) and
                    case['frequencies'][f] > case['frequencies'][f - 1]):
                active = False
        if ready and active:
            # Independent formulation: normalize along each channel-vector axis.
            pa = track_a[f - 1] / np.linalg.norm(track_a[f - 1], axis=0)
            pb = track_b[f - 1] / np.linalg.norm(track_b[f - 1], axis=0)
            ca = a[f] / np.linalg.norm(a[f], axis=0)
            cb = b[f] / np.linalg.norm(b[f], axis=0)
            h = (pa.conj().T @ ca + pb.conj().T @ cb) / 2.0
            quality = np.abs(h)
            row, col = linear_sum_assignment(-quality)
            permutation = col[np.argsort(row)]
            selected = quality[np.arange(nm), permutation]
            if selected.min() >= case['minimum'] and selected.min() > 64 * np.finfo(float).eps:
                source[f] = permutation + 1
                modes[f] = case['values'][f, permutation]
                overlap[f] = selected
                matched[f] = True
                gauge = np.conj(h[np.arange(nm), permutation]) / selected
                track_a[f] = a[f, :, permutation].T * gauge
                track_b[f] = b[f, :, permutation].T * gauge
        active = bool(ready)
    return dict(source=source, overlap=overlap, matched=matched,
                modes=modes, a=track_a, b=track_b)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261030)
    cases = [generate(rng, i) for i in range(80)]
    lines = [str(len(cases))]
    for case in cases:
        lines.append(f"{case['nf']} {case['na']} {case['nb']} {case['nm']} {case['minimum']:.17g}")
        for f in range(case['nf']):
            lines.append(f"{case['frequencies'][f]:.17g} {'T' if case['valid'][f] else 'F'}")
            lines.append(' '.join(f'{x:.17g}' for x in case['values'][f]))
            lines.append(' '.join('T' if x else 'F' for x in case['stable'][f]))
            for k in range(case['nm']):
                for vec in (case['a'][f, :, k], case['b'][f, :, k]):
                    lines.extend(f'{x.real:.17g} {x.imag:.17g}' for x in vec)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
    result = subprocess.run([str(args.executable)], input='\n'.join(lines) + '\n',
                            text=True, capture_output=True, check=True, env=env)
    output = iter(result.stdout.splitlines())
    max_error = 0.0
    total = 0
    for index, case in enumerate(cases, 1):
        header = next(output).split()
        if header != ['CASE', str(index), '0']:
            raise AssertionError(f'case {index}: unexpected output {header!r}')
        truth = reference(case)
        for f in range(case['nf']):
            for k in range(case['nm']):
                parts = next(output).split()
                row, col, source, matched = parts[:4]
                if (int(row), int(col)) != (f + 1, k + 1):
                    raise AssertionError('out-of-order output')
                if int(source) != truth['source'][f, k]:
                    raise AssertionError(f'case {index} bin {f} mode {k}: source label differs')
                if (matched == 'T') != truth['matched'][f, k]:
                    raise AssertionError(f'case {index} bin {f} mode {k}: matching mask differs')
                sim, value = map(float, parts[4:6])
                exp_sim = truth['overlap'][f, k]
                if np.isnan(exp_sim):
                    if not np.isnan(sim):
                        raise AssertionError('restart must have NaN overlap')
                else:
                    max_error = max(max_error, abs(sim - exp_sim))
                max_error = max(max_error, abs(value - truth['modes'][f, k]))
                numbers = np.array([float(x) for x in parts[6:]])
                expected_count = 2 * (case['na'] + case['nb'])
                if len(numbers) != expected_count:
                    raise AssertionError('unexpected coefficient output length')
                z = numbers[::2] + 1j * numbers[1::2]
                expected = np.r_[truth['a'][f, :, k], truth['b'][f, :, k]]
                max_error = max(max_error, np.max(np.abs(z - expected)))
                total += 1
    if max_error > 2.0e-12:
        raise AssertionError(f'maximum absolute error {max_error:.3g} exceeded tolerance')
    print(f'PASS SciPy/NumPy canonical tracking comparison: {len(cases)} cases, '
          f'{total} modes, max absolute error {max_error:.3g}')


if __name__ == '__main__':
    main()
