#!/usr/bin/env python3
"""Randomized reference tests of CSR expm_multiply versus SciPy and dense expm.

Exercises real CSR operators, vector/matrix right-hand sides, duplicate COO
entries, trace shifting, negative times, forward/reverse time intervals and
endpoint handling. The Fortran executable processes all cases in one run.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply as scipy_expm_multiply

ROOT = Path(__file__).resolve().parents[2]


def make_case(rng: np.random.Generator, i: int):
    n = int(rng.integers(1, 24))
    p = int(rng.integers(1, 5))
    d = rng.standard_normal((n, n)) * (rng.random((n, n)) < .28)
    d *= float(rng.choice([.05, .25, .9, 1.5]))
    kind = i % 11
    if kind == 0:  # A diagonal, all shifts detected accurately.
        d = np.eye(n) * rng.uniform(-9, 9)
    elif kind == 1:  # Large scalar shift should not dictate Taylor degree.
        d += np.eye(n) * 20
    elif kind == 2:  # Nilpotent operator, possibly strongly nonnormal.
        d = np.triu(d, k=1) * 4
    elif kind == 3:  # Skew-symmetric operator (oscillatory exponential).
        d = d - d.T
    elif kind == 4:  # Negative diagonal, contraction.
        d -= np.eye(n) * 4
    elif kind == 5:  # Zero matrix, including completely empty sparsity.
        d = np.zeros((n, n))
    elif kind == 6:  # Non-normal Jordan-like chain.
        d = np.eye(n) * rng.uniform(-1, 1)
        d += np.diag(rng.uniform(-4, 4, max(n-1, 0)), k=1)
    elif kind == 7:  # Sparse symmetric indefinite.
        d = 0.5 * (d + d.T)
    elif kind == 8:
        d *= 0.1
    elif kind == 9:
        d += np.eye(n) * float(rng.choice([-18, 18]))

    if i % 17 == 0:
        p = 1
    b = rng.standard_normal((n, p))
    if i % 15 == 0:
        b.fill(0.0)
    if i % 19 == 0:
        b *= 1e-120
    if i % 29 == 0:
        b *= 1e-250

    time = float(rng.uniform(-2, 2))
    start = float(rng.uniform(-1.5, 1.5))
    stop = float(rng.uniform(-1.5, 1.5))
    if i % 9 == 0:
        time = 0.0
    if i % 7 == 0:
        start = stop
    num = int(rng.integers(1, 7))
    endpoint = i % 2 == 0
    trace_value = float(np.trace(d)) if i % 3 == 0 else 0.0
    use_trace = i % 3 != 1
    return d, b, time, start, stop, num, endpoint, use_trace, trace_value


def fmt(x: float) -> str:
    return f'{float(x):.17g}'


def assert_accurate(actual, expected, message: str):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    if not np.all(np.isfinite(actual)):
        raise AssertionError(f'{message}: nonfinite result')
    scale = max(np.max(np.abs(expected), initial=0.0), 1e-300)
    max_difference = np.max(np.abs(actual - expected), initial=0.0)
    if max_difference > scale * 2e-10:
        raise AssertionError(
            f'{message}: max_err={max_difference:.5g}, scale={scale:.5g}, '
            f'relative={max_difference/scale:.5g}')
    return float(max_difference / scale)


def parse_record(line: str, n: int, case_num: int):
    try:
        vals = np.fromstring(line, sep=' ')
    except ValueError as exc:
        raise AssertionError(f'case {case_num}: invalid Fortran output: {line!r}') from exc
    if vals.size != n:
        raise AssertionError(f'case {case_num}: expected {n} values, got {vals.size}: {line!r}')
    return vals


def check_mode(mode: str, count: int, seed: int) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-expm-') as temp:
        build = Path(temp)
        flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
        cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
        objs = []
        for name, rel in [('kinds', 'src/kinds.f90'),
                          ('sparse', 'src/sparse.f90'),
                          ('lsmr', 'src/optimize/lsmr.f90'),
                          ('sparse_linalg', 'src/sparse_linalg.f90')]:
            obj = build / (name + '.o')
            subprocess.run([*cc, '-c', str(ROOT / rel), '-o', str(obj)],
                           cwd=build, check=True)
            objs.append(obj)
        executable = build / 'sparse_expm_reference'
        subprocess.run([*cc, str(ROOT / 'tools/sparse_linalg_reference'
                                      '/sparse_expm_reference_driver.f90'),
                        *map(str, objs), '-o', str(executable)],
                       cwd=build, check=True)
        rng = np.random.default_rng(seed)
        cases = [make_case(rng, i) for i in range(count)]
        lines = [str(len(cases))]
        for a, b, t, start, stop, num, endpoint, use_trace, trace in cases:
            rows, cols = np.nonzero(a)
            triples = [(int(r+1), int(c+1), float(v)) for r, c, v in
                       zip(rows, cols, a[rows, cols], strict=True)]
            if triples:
                r, c, v = triples[0]
                triples[0] = (r, c, v * .25)
                triples.extend([(r, c, v * .75), (r, c, 0.0)])
            n, p = b.shape
            lines.append(f'{n} {p} {len(triples)} {fmt(t)} {fmt(start)} '
                         f'{fmt(stop)} {num} {int(endpoint)} {int(use_trace)} '
                         f'{fmt(trace)}')
            lines.extend(f'{r} {c} {fmt(v)}' for r, c, v in triples)
            lines.extend(fmt(x) for x in b.T.ravel())

        response = subprocess.run([str(executable)], cwd=build,
                                  input='\n'.join(lines) + '\n', text=True,
                                  capture_output=True, check=True)
        out = iter(response.stdout.splitlines())
        worst_single = worst_grid = 0.0
        interval_disagreements = 0
        interval_exceptions = 0
        reference_disagreements = 0
        for case_num, (a, b, t, start, stop, num, endpoint,
                       use_trace, trace) in enumerate(cases, 1):
            assert next(out).strip() == f'CASE {case_num}'
            state, grid_state, products, grid_products = map(int, next(out).split())
            if state != 0 or grid_state != 0:
                raise AssertionError(f'case {case_num}: status {state} {grid_state} '
                                     f'{next(out)}')
            n, p = b.shape
            single = np.column_stack([
                parse_record(next(out), n, case_num) for _ in range(p)])
            grid = np.empty((num, n, p))
            for grid_index in range(num):
                grid[grid_index] = np.column_stack([
                    parse_record(next(out), n, case_num) for _ in range(p)])
            kwargs = {'traceA': trace} if use_trace else {}
            sp = csr_matrix(a)
            ref = scipy_expm_multiply(t * sp, b,
                    **({'traceA': t*trace} if use_trace else {}))
            expected_single = expm(t*a) @ b
            worst_single = max(worst_single,
                assert_accurate(single, expected_single, f'case {case_num} single'))
            if not np.allclose(ref,expected_single,rtol=2e-10,atol=0):
                reference_disagreements += 1
            # Dense small-matrix actions serve as an independent reference;
            # SciPy is also checked, but inaccurate traceA hints can degrade
            # its accuracy under severe exponential cancellation.
            sample_times = np.linspace(start, stop, num, endpoint=endpoint)
            expected = np.stack([expm(float(tt)*a) @ b for tt in sample_times],axis=0)
            scipy_direct = np.stack([
                scipy_expm_multiply(float(tt)*sp, b,
                    **({'traceA': float(tt)*trace} if use_trace else {}))
                for tt in sample_times], axis=0)
            if not np.allclose(scipy_direct,expected,rtol=2e-10,atol=0):
                reference_disagreements += 1
            if num >= 2:
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', RuntimeWarning)
                        scipy_interval = scipy_expm_multiply(
                            sp,b,start=start,stop=stop,num=num,endpoint=endpoint,**kwargs)
                    if not np.allclose(scipy_interval,expected,rtol=2e-10,atol=0):
                        interval_disagreements += 1
                except (ZeroDivisionError, ValueError, FloatingPointError):
                    interval_exceptions += 1
            worst_grid = max(worst_grid,
                             assert_accurate(grid, expected, f'case {case_num} grid'))
            if t == 0.0 and products != 0:
                raise AssertionError('time-zero calculation must not multiply CSR')
        print(f'{mode}: {count} expm_multiply cases passed; '
              f'max rel error {worst_single:.3g}, grid {worst_grid:.3g}; '
              f'SciPy interval/dense disagreements {interval_disagreements}; '
              f'SciPy interval exceptions {interval_exceptions}; '
              f'SciPy direct/dense disagreements {reference_disagreements}', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['checked', 'optimized', 'both'], default='both')
    parser.add_argument('--cases', type=int, default=180)
    parser.add_argument('--seed', type=int, default=20261104)
    args = parser.parse_args()
    modes = ['checked', 'optimized'] if args.mode == 'both' else [args.mode]
    for mode in modes:
        check_mode(mode, args.cases, args.seed)


if __name__ == '__main__':
    main()
