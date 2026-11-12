#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Compile/check SciFort QAP against SciPy and independent permutation tests.

Run from any directory::

    python validation/compare_optimize_quadratic_assignment_scipy.py
    python validation/compare_optimize_quadratic_assignment_scipy.py --optimize

Only requires gfortran, NumPy, and SciPy. It does not modify the source tree.
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.optimize import quadratic_assignment as scipy_qap


def qap_score(a: np.ndarray, b: np.ndarray, perm: np.ndarray) -> float:
    return float(np.sum(a * b[np.ix_(perm, perm)]))


def generate_cases() -> list[dict]:
    rng = np.random.default_rng(20261110)
    cases = []
    for n in (2, 3, 4, 5, 6, 7, 9, 12, 16):
        for symmetric in (False, True):
            for rep in range(4):
                a = rng.integers(-10, 11, size=(n, n)).astype(float)
                b = rng.integers(-10, 11, size=(n, n)).astype(float)
                if symmetric:
                    a = a + a.T
                    b = b + b.T
                target = rng.permutation(n)
                nseed = (rep % 4) if n > 4 else rep % max(1, n)
                seeds = np.array([(i, target[i]) for i in rng.permutation(n)[:nseed]], dtype=int)
                seeds = seeds.reshape(-1, 2)
                fixed_rows = set(seeds[:, 0])
                fixed_cols = set(seeds[:, 1])
                remaining_rows = [i for i in range(n) if i not in fixed_rows]
                remaining_cols = [j for j in range(n) if j not in fixed_cols]
                for maximize in (False, True):
                    for meth in (1, 2):
                        guesses = np.empty((0, 2), dtype=int)
                        if meth == 2 and rep % 2 == 1:
                            shuffled = rng.permutation(remaining_cols)
                            guesses = np.array(list(zip(remaining_rows, shuffled)), dtype=int)
                            guesses = guesses.reshape(-1, 2)
                        p0 = None
                        if meth == 1 and rep == 3 and len(remaining_rows) > 0:
                            # A genuine doubly stochastic starting matrix.
                            p0 = rng.uniform(0.1, 1.0, size=(len(remaining_rows),) * 2)
                            for _ in range(120):
                                p0 /= p0.sum(axis=0, keepdims=True)
                                p0 /= p0.sum(axis=1, keepdims=True)
                        cases.append(dict(n=n, a=a, b=b, meth=meth, maximize=maximize,
                                          seeds=seeds, guesses=guesses, p0=p0, maxiter=30 if meth == 1 else 10000))
    return cases


def case_input(cases: list[dict]) -> str:
    lines = []
    for c in cases:
        n, seeds, guesses, p0 = c['n'], c['seeds'], c['guesses'], c['p0']
        lines.append(f"{n} {c['meth']} {int(c['maximize'])} {len(seeds)} {len(guesses)} "
                     f"{c['maxiter']} {int(p0 is not None)} 0.03")
        for matrix in (c['a'], c['b']):
            lines.extend(' '.join(f'{v:.17g}' for v in row) for row in matrix)
        lines.extend(' '.join(str(v+1) for v in row) for row in seeds)
        lines.extend(' '.join(str(v+1) for v in row) for row in guesses)
        if p0 is not None:
            lines.extend(' '.join(f'{v:.17g}' for v in row) for row in p0)
    lines.append('-1 0 0 0 0 0 0 0.03')
    return '\n'.join(lines) + '\n'


def validate(cases: list[dict], stdout: str) -> dict:
    lines = [line.strip() for line in stdout.splitlines()]
    assert len(lines) == 2 * len(cases), f'expected {2*len(cases)} output lines, got {len(lines)}'
    counts = dict(faq=0, two_opt=0, scipy_faq_equal=0, scipy_faq_diff=0, scipy_2opt_quality_better=0,
                  scipy_2opt_quality_worse=0, scipy_2opt_quality_equal=0)
    checks = 0
    for i, case in enumerate(cases):
        head = lines[2*i].split()
        tail = lines[2*i+1].split()
        assert head[0] == 'RESULT' and tail[0] == 'PERM', (head, tail)
        status, nit, value = int(head[1]), int(head[2]), float(head[3])
        perm = np.array(list(map(int, tail[1:])), dtype=int) - 1
        n, a, b = case['n'], case['a'], case['b']
        assert status == 0, (i, case, head)
        assert len(perm) == n and sorted(perm) == list(range(n)), (i, perm)
        assert all(perm[r] == c for r, c in case['seeds'])
        independent = qap_score(a, b, perm)
        assert np.isclose(value, independent, atol=1e-9, rtol=1e-12), (i, value, independent)
        checks += 4
        opt_sign = -1 if case['maximize'] else 1
        if case['meth'] == 1:
            counts['faq'] += 1
            assert nit <= case['maxiter']
            options = dict(maximize=case['maximize'], partial_match=case['seeds'][np.argsort(case['seeds'][:, 0])],
                           maxiter=case['maxiter'], tol=0.03, rng=np.random.default_rng(33))
            if case['p0'] is not None:
                options['P0'] = case['p0']
            reference = scipy_qap(a, b, method='faq', options=options)
            ref_perm = reference.col_ind
            assert all(ref_perm[r] == c for r, c in case['seeds'])
            assert np.isclose(qap_score(a, b, ref_perm), reference.fun)
            if np.isclose(value, reference.fun):
                counts['scipy_faq_equal'] += 1
            else:
                counts['scipy_faq_diff'] += 1
            checks += 3
        else:
            counts['two_opt'] += 1
            rows = [r for r in range(n) if r not in set(case['seeds'][:, 0])]
            for r, s in itertools.combinations(rows, 2):
                swapped = perm.copy()
                swapped[r], swapped[s] = swapped[s], swapped[r]
                assert opt_sign * (qap_score(a, b, swapped) - value) >= -1e-8, (i, r, s)
                checks += 1
            # Reference is stochastic; compare quality, not exact permutation/iteration count.
            reference = scipy_qap(a, b, method='2opt', options={
                'maximize':case['maximize'], 'partial_match':case['seeds'][np.argsort(case['seeds'][:, 0])],
                'rng':np.random.default_rng(33)})
            assert all(reference.col_ind[r] == c for r, c in case['seeds'])
            assert np.isclose(qap_score(a, b, reference.col_ind), reference.fun)
            rel = opt_sign * (value - reference.fun)
            if np.isclose(rel, 0.0):
                counts['scipy_2opt_quality_equal'] += 1
            elif rel < 0:
                counts['scipy_2opt_quality_better'] += 1
            else:
                counts['scipy_2opt_quality_worse'] += 1
            checks += 2
    counts['assertions'] = checks
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optimize', action='store_true', help='compile with -O3 rather than runtime checks')
    args = parser.parse_args()
    compiler = shutil.which('gfortran')
    if not compiler:
        raise SystemExit('gfortran not available')
    root = Path(__file__).resolve().parents[1]
    files = [root/'src/kinds.f90', root/'src/optimize/assignment.f90',
             root/'src/optimize/quadratic_assignment.f90',
             root/'validation/optimize_quadratic_assignment_driver.f90']
    flags = ['-std=f2008', '-O3'] if args.optimize else [
        '-std=f2008', '-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero']
    cases = generate_cases()
    with tempfile.TemporaryDirectory(prefix='scifort_qap_') as tmp:
        exe = Path(tmp) / 'qap_driver'
        subprocess.run([compiler, *flags, '-J', tmp, '-I', tmp, *map(str, files), '-o', str(exe)],
                       check=True, capture_output=True, text=True)
        process = subprocess.run([str(exe)], input=case_input(cases), check=True,
                                 text=True, capture_output=True)
    counts = validate(cases, process.stdout)
    print('GNU Fortran/SciPy quadratic assignment checks:', 'optimized' if args.optimize else 'checked')
    for label, value in counts.items():
        print(f'  {label}: {value}')
    print(f'  SciPy version: {__import__("scipy").__version__}')


if __name__ == '__main__':
    main()
