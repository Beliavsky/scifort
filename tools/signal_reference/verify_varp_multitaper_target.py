#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent SciPy VAR(p) multitaper spectral and canonical references.

The reference solves the companion Lyapunov equation with SciPy, sums all
sample pairs (O(K*L**2)) for each independent taper, then mixes spectral
matrices BEFORE whitening and singular-value decomposition. Both exact
frequency-band integrals and the complex phase are checked.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.signal.windows import dpss, hann
from verify_varp_taper_spectrum import reference
from verify_varp_canonical_target import canonical


def run(executable: str, cases: int):
    rng = np.random.default_rng(20261031)
    scenarios = []
    for case in range(cases):
        n, p, length = 5, 1 + case % 3, 6 + case % 15
        ntaper, nf, nb = 1 + case % 4, 2 + case % 4, 2 + case % 3
        a = rng.normal(size=(n, n, p))
        a *= 0.35 / sum(np.linalg.norm(a[:, :, j], 2) for j in range(p))
        z = rng.normal(size=(n, n))
        sigma = z @ z.T + 0.6 * np.eye(n)
        if case % 3 == 0:
            # SciPy's discrete prolate spheroidal sequences, independent of
            # SciFort's native dpss generator.
            tapers = dpss(length, NW=min(2.5, length / 3),
                          Kmax=ntaper, sym=False).T
        elif case % 3 == 1:
            tapers = rng.normal(size=(length, ntaper))
        else:
            tapers = hann(length, sym=False)[:, None] * (
                1 + rng.random((length, ntaper)))
        tapers = tapers * rng.uniform(0.7, 2.7, size=(1, ntaper))
        freqs = rng.uniform(-1.5, 1.8, size=nf)
        freqs[0] = 0.0
        bands = np.column_stack([
            rng.uniform(-0.6, -0.05, size=nb),
            rng.uniform(0.03, 0.75, size=nb),
        ])
        if nb > 1:
            bands[1] = [0.17, 0.17001]
        weights = rng.uniform(0.02, 4.0, size=ntaper)
        if case % 6 == 0 and ntaper > 1:
            weights[-1] = 0.0
        normalized = case % 2 == 0
        average = case % 3 == 0
        scenarios.append((a, sigma, tapers, freqs, bands, weights,
                          normalized, average))

    lines = [str(cases)]
    for a, sigma, tapers, freqs, bands, weights, norm, avg in scenarios:
        n, _, p = a.shape
        length, ntaper = tapers.shape
        lines += [f'{n} {p} {length} {ntaper} {len(freqs)} {len(bands)}',
                  ' '.join(f'{v:.17g}' for v in a.flatten(order='F')),
                  ' '.join(f'{v:.17g}' for v in sigma.flatten(order='F')),
                  ' '.join(f'{v:.17g}' for v in tapers.flatten(order='F')),
                  ' '.join(f'{v:.17g}' for v in freqs),
                  ' '.join(f'{v:.17g}' for v in bands.flatten(order='C')),
                  ' '.join(f'{v:.17g}' for v in weights),
                  f'{int(norm)} {int(avg)}']
    process = subprocess.run([executable], input='\n'.join(lines) + '\n',
                             text=True, capture_output=True, check=True)
    rows = process.stdout.splitlines()
    if len(rows) != cases:
        raise AssertionError((len(rows), cases, process.stderr[:2000]))
    max_spectral_error, max_canonical_error = 0.0, 0.0
    ncoeffs = 0
    for case, (row, scenario) in enumerate(zip(rows, scenarios)):
        a, sigma, tapers, freqs, bands, weights, norm, avg = scenario
        nf, nb, n = len(freqs), len(bands), a.shape[0]
        raw = np.fromstring(row, sep=' ')
        if raw.size < 2 or not np.array_equal(raw[:2], [0, 0]):
            raise AssertionError(('invalid result', case, raw[:3]))
        expected_size = 2 + (nf + nb) * (2 * n * n + 3)
        if raw.size != expected_size:
            raise AssertionError(('length', case, raw.size, expected_size))
        cut = 2 + (nf + nb) * 2 * n * n
        # The Fortran driver writes (frequency, column, row) iteration order.
        sp = raw[2:cut].reshape(nf + nb, n, n, 2).transpose(0, 2, 1, 3)
        measured = sp[..., 0] + 1j * sp[..., 1]
        measured_coherence = raw[cut:].reshape(nf + nb, 3)
        weighted = np.zeros_like(measured)
        for k in range(tapers.shape[1]):
            if weights[k] == 0:
                continue
            _, point, integrated = reference(
                a, sigma, tapers[:, k], freqs, bands,
                average=avg, normalized=norm)
            weighted[:nf] += (weights[k] / weights.sum()) * point
            weighted[nf:] += (weights[k] / weights.sum()) * integrated
        target_coherence = canonical(weighted)
        error_s = float(np.max(abs(measured - weighted)))
        error_c = float(np.max(abs(measured_coherence - target_coherence)))
        if error_s > 4e-10 or error_c > 4e-10:
            raise AssertionError((case, 'spectral', error_s, 'canonical', error_c,
                                  measured_coherence, target_coherence))
        max_spectral_error = max(max_spectral_error, error_s)
        max_canonical_error = max(max_canonical_error, error_c)
        ncoeffs += 2 * (nf + nb) * n * n

    print(f'PASS {cases} independent DPSS/random VAR(p) multitaper cases; '
          f'{ncoeffs} real/imag spectral components; '
          f'max spectral error {max_spectral_error:.3e}, '
          f'canonical error {max_canonical_error:.3e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases', type=int, default=40)
    args = parser.parse_args()
    run(args.executable, args.cases)
