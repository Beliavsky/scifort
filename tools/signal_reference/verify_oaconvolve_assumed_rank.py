#!/usr/bin/env python3
"""Seeded SciPy comparisons for rank-1/2/3/5/15 Fortran oaconvolve_nd.

Build the Fortran driver before running; pass its executable as argument.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import tempfile

import numpy as np
from scipy.signal import oaconvolve


def make_cases(seed: int, per_rank: int):
    rng = np.random.default_rng(seed)
    cases = []
    for r in [1, 2, 3, 5, 15]:
        for idx in range(per_rank):
            n, m = np.ones(r, dtype=int), np.ones(r, dtype=int)
            candidate_axes = list(range(r))
            rng.shuffle(candidate_axes)
            if r == 15:
                # Include the final axis frequently without growing the input
                # exponentially in singleton-broadcast dimensions.
                candidate_axes[0] = r - 1
                candidate_axes = list(dict.fromkeys(candidate_axes + list(range(r))))
            num_axes = int(rng.integers(1, min(r, 3) + 1))
            axes = candidate_axes[:num_axes]
            mode = ['full', 'same', 'valid'][idx % 3]
            backwards = (idx % 2) == 1
            blocking = (idx % 4) != 0
            for j, ax in enumerate(axes):
                if j == 0 and blocking:
                    small = int(rng.integers(2, 9))
                    big = int(rng.integers(9 * small, 16 * small))
                else:
                    small = int(rng.integers(1, 4))
                    big = int(rng.integers(small, small + 4))
                n[ax], m[ax] = (small, big) if backwards else (big, small)
            if r > 1:
                others = [j for j in range(r) if j not in axes]
                for j in rng.choice(others, size=min(len(others), 2), replace=False):
                    n[j], m[j] = [(2, 1), (1, 3), (2, 2)][int(rng.integers(3))]
            dtype = idx % 4  # RR, CC, RC, CR
            x = np.asfortranarray(rng.standard_normal(n))
            z = np.asfortranarray(rng.standard_normal(m))
            if dtype in (1, 3):
                x = x + 1j*rng.standard_normal(n)
            if dtype in (1, 2):
                z = z + 1j*rng.standard_normal(m)
            try:
                expected = oaconvolve(x, z, mode=mode, axes=tuple(axes))
            except ValueError as exc:
                raise AssertionError((r, n, m, axes, mode)) from exc
            cases.append((r, dtype, mode, axes, x, z, expected))
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver', type=Path, help='compiled Fortran test driver')
    parser.add_argument('--cases-per-rank', type=int, default=80)
    args = parser.parse_args()
    cases = make_cases(48017, args.cases_per_rank)
    with tempfile.TemporaryDirectory(prefix='scifort-oa-') as folder:
        inp = Path(folder)/'cases.txt'
        out = Path(folder)/'out.txt'
        with inp.open('w', encoding='ascii') as f:
            for r, dtype, mode, axes, x, z, expected in cases:
                f.write(f'{r} {dtype} {["full", "same", "valid"].index(mode)+1} {len(axes)}\n')
                for dims in (x.shape, z.shape, expected.shape):
                    f.write(' '.join(map(str, dims))+'\n')
                f.write(' '.join(str(k+1) for k in axes)+'\n')
                for values in (x, z):
                    cvalues = values.ravel(order='F').astype(complex)
                    for value in cvalues:
                        f.write(f'{value.real:.17g} {value.imag:.17g}\n')
        with inp.open('rb') as source, out.open('wb') as dest:
            subprocess.run([str(args.driver.resolve())], stdin=source, stdout=dest, check=True)
        lines = out.read_text().splitlines()
    line = 0
    worst = 0.0
    for ci, (r, dtype, mode, axes, x, z, ref) in enumerate(cases):
        status, count = map(int, lines[line].split())
        line += 1
        if status != 0 or count != ref.size:
            raise AssertionError(f'case {ci}: status={status}, count={count}, expected={ref.size}; '
                                 f'rank={r}, shapes={x.shape}/{z.shape}, axes={axes}, mode={mode}')
        received = np.array([np.fromstring(s, sep=' ') for s in lines[line:line+count]])
        line += count
        if received.shape != (count, 2):
            raise AssertionError(f'case {ci}: invalid output: {received.shape}')
        result = (received[:, 0]+1j*received[:, 1]).reshape(ref.shape, order='F')
        error = float(np.max(np.abs(result-ref)))
        worst = max(worst, error)
        if error > 2e-10:
            raise AssertionError(f'case {ci}: discrepancy {error:.3g}, rank={r}, '
                                 f'shapes={x.shape}/{z.shape}, axes={axes}, mode={mode}, dtype={dtype}')
    assert line == len(lines), f'{len(lines)-line} extra output lines'
    print(f'PASS: {len(cases)} SciPy oaconvolve cases, ranks 1,2,3,5,15, '
          f'maximum absolute error {worst:.3e}')


if __name__ == '__main__':
    main()
