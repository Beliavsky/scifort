#!/usr/bin/env python3
"""Independent SciPy comparisons for SciFort linalg.cdf2rdf.

Usage:
    python reference/check_linalg_cdf2rdf.py --driver /path/to/driver

The driver is compiled separately using the original SciFort sources. Results
are compared with scipy.linalg.cdf2rdf, not with reproduced SciPy code.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass

import numpy as np
from scipy import linalg


@dataclass
class Case:
    kind: int
    w: np.ndarray
    v: np.ndarray
    invalid: bool = False
    nonfinite: bool = False


def cases(seed: int = 46723) -> list[Case]:
    rng = np.random.default_rng(seed)
    result: list[Case] = []
    for n in range(0, 19):
        for trial in range(29):
            # Eigenpairs of actual real matrices, including many complex pairs.
            a = rng.standard_normal((n, n))
            if trial % 7 == 0:
                a = np.triu(a)
            if trial % 9 == 0:
                a = np.eye(n) * (trial + 1)
            if n == 0:
                w = np.empty(0, dtype=np.complex128)
                v = np.empty((0, 0), dtype=np.complex128)
            else:
                w, v = linalg.eig(a)
            result.append(Case(1, w, v))
    for n in range(0, 20):
        for trial in range(27):
            # Construct conjugate eigenpairs in occurrence order, optionally
            # separated by real eigenvalues, with exactly matching vectors.
            scale = (1e-70, 1, 1e70)[trial % 3]
            w = np.empty(n, dtype=np.complex128)
            v = np.empty((n, n), dtype=np.complex128)
            i = 0
            while i < n:
                if i + 1 == n or (trial + i) % 5 == 0:
                    w[i] = scale * rng.standard_normal()
                    v[:, i] = scale * rng.standard_normal(n)
                    i += 1
                    continue
                imag = (1 + abs(rng.standard_normal())) * scale
                if trial % 2:
                    imag = -imag
                w[i] = scale * rng.standard_normal() + 1j * imag
                # Optionally leave a real eigenvalue between the pair members.
                k = i + 2 if i + 2 < n and (trial + i) % 3 == 0 else i + 1
                if k == i + 2:
                    w[i+1] = scale * rng.standard_normal()
                    v[:, i+1] = scale * rng.standard_normal(n)
                w[k] = np.conj(w[i])
                v[:, i] = scale * (rng.standard_normal(n) + 1j*rng.standard_normal(n))
                v[:, k] = np.conj(v[:, i])
                i = k + 1
            result.append(Case(1, w, v))
            if trial % 6 == 0:
                # When w is real, SciPy accepts real and complex eigenvectors.
                real_w = np.arange(n, dtype=np.float64) + trial
                result.append(Case(2, real_w, rng.standard_normal((n, n))))
                result.append(Case(3, real_w, rng.standard_normal((n,n)) +
                                   1j*rng.standard_normal((n,n))))
    for n in range(2, 18):
        for trial in range(10):
            w = np.arange(n, dtype=np.complex128)
            v = np.eye(n, dtype=np.complex128)
            w[0] = 1 + 3j
            w[1] = 1 - 3j
            if trial % 3 == 0:
                w[1] = 1 - 4j
            elif trial % 3 == 1:
                w[1] = 2 - 3j
            else:
                w[1] = 1 + 3j
            result.append(Case(1, w, v, invalid=True))
    result.append(Case(1, np.array([1+2j, 2+3j, 1-2j, 2-3j]),
                       np.eye(4, dtype=np.complex128), invalid=True))
    result.append(Case(1, np.array([1+2j, 0j, 5+0j]),
                       np.eye(3, dtype=np.complex128), invalid=True))
    result.append(Case(1, np.array([np.nan+0j]),
                       np.eye(1, dtype=np.complex128), nonfinite=True))
    result.append(Case(1, np.array([np.inf+0j]),
                       np.eye(1, dtype=np.complex128), nonfinite=True))
    result.append(Case(1, np.array([1+0j]),
                       np.array([[np.nan+0j]]), nonfinite=True))
    result.append(Case(2, np.array([1.]), np.array([[np.inf]]), nonfinite=True))
    return result


def serialized(case: Case) -> list[str]:
    n = len(case.w)
    lines = [f"{case.kind} {n}"]
    if case.kind == 1:
        lines += [f"{z.real:.17e} {z.imag:.17e}" for z in case.w]
    else:
        lines += [f"{x:.17e}" for x in case.w]
    if case.kind == 2:
        lines += [f"{x:.17e}" for x in case.v.flatten(order='F')]
    else:
        lines += [f"{z.real:.17e} {z.imag:.17e}" for z in case.v.flatten(order='F')]
    return lines


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--driver', required=True)
    parser.add_argument('--seed', type=int, default=46723)
    args = parser.parse_args()
    test_cases = cases(args.seed)
    message = '\n'.join([str(len(test_cases))] +
                        [line for case in test_cases for line in serialized(case)]) + '\n'
    done = subprocess.run([args.driver], input=message, capture_output=True,
                          text=True, check=True)
    output = iter(done.stdout.splitlines())
    maximum_wr = 0.0
    maximum_vr = 0.0
    maximum_residual = 0.0
    n_rejected = 0
    n_valid = 0
    for idx, case in enumerate(test_cases, 1):
        header = next(output).split()
        assert header[:2] == ['CASE', str(idx)], (idx, header)
        status = int(header[2])
        n = len(case.w)
        fw = np.array([float(next(output)) for _ in range(n*n)]).reshape((n,n), order='F')
        fv = np.array([float(next(output)) for _ in range(n*n)]).reshape((n,n), order='F')
        if case.invalid or case.nonfinite:
            expected = 4 if case.nonfinite else 1
            assert status == expected, (idx, status, expected)
            assert np.all(fw == 0) and np.all(fv == 0), idx
            n_rejected += 1
            continue
        assert status == 0, (idx, status, case.w)
        sw, sv = linalg.cdf2rdf(case.w, case.v)
        def scaled_error(x: np.ndarray, y: np.ndarray) -> float:
            if x.size == 0:
                return 0.0
            scale = max(np.max(np.abs(x)), np.max(np.abs(y)), np.finfo(float).tiny)
            return float(np.max(np.abs(x/scale-y/scale)))
        ew = scaled_error(fw, sw)
        ev = scaled_error(fv, sv)
        maximum_wr = max(maximum_wr, ew)
        maximum_vr = max(maximum_vr, ev)
        assert ew < 3e-13, (idx, 'wr', ew)
        assert ev < 3e-13, (idx, 'vr', ev)
        # Extra mathematical test for real eigenstructure recovered by eig.
        # The constructed random w/v cases are not eigenpairs of a matrix.
        if idx <= 19*29 and n:
            vv = np.linalg.norm(fv)
            if vv:
                a = case.v @ np.diag(case.w) @ np.linalg.pinv(case.v)
                if np.linalg.cond(case.v) < 1e8 and np.max(np.abs(a.imag)) < 1e-6:
                    rel = np.linalg.norm(a.real@fv-fv@fw) / max(
                        np.linalg.norm(a.real)*vv, 1e-300)
                    maximum_residual = max(maximum_residual, rel)
                    assert rel < 3e-11, (idx, 'residual', rel)
        n_valid += 1
    assert next(output, None) is None
    print(f'PASS cdf2rdf SciPy comparisons: {len(test_cases)} '
          f'({n_valid} valid, {n_rejected} expected rejections)')
    print(f'maximum scaled matrix differences: wr={maximum_wr:.3e}, '
          f'vr={maximum_vr:.3e}; maximum eig residual={maximum_residual:.3e}')


if __name__ == '__main__':
    main()
