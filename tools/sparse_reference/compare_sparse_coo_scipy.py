#!/usr/bin/env python3
"""Compare SciFort's one-based COO format and sparse_find with SciPy 1.17.

Compiles the independent Fortran source and runs batched randomized cases.
Requires NumPy, SciPy, and gfortran; no network access or FPM is necessary.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]


def generate_cases(n: int, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    cases = []
    for i in range(n):
        m, nc = (int(rng.integers(0, 11)), int(rng.integers(0, 11)))
        if i % 17 == 0:
            m = 0
        if i % 23 == 0:
            nc = 0
        count = 0 if min(m, nc) == 0 else int(rng.integers(0, min(100, 4 * m * nc) + 1))
        rows = rng.integers(0, m, size=count)
        cols = rng.integers(0, nc, size=count)
        data = rng.integers(-8, 9, size=count).astype(float)
        x = rng.integers(-5, 6, size=nc).astype(float)
        cases.append(dict(m=m, n=nc, rows=rows, cols=cols, data=data, x=x))
    return cases


def compare(mode: str, cases: list[dict]) -> None:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-coo-') as td:
        work = Path(td)
        cmd = ['gfortran', '-std=f2018', *flags, '-J', str(work), '-I', str(work)]
        for name in ('kinds', 'sparse', 'sparse_coo'):
            subprocess.run([*cmd, '-c', str(ROOT / 'src' / f'{name}.f90'),
                            '-o', str(work / f'{name}.o')], cwd=work, check=True)
        exe = work / 'coo_driver'
        subprocess.run([*cmd, str(ROOT / 'tools/sparse_reference/sparse_coo_reference_driver.f90'),
                        *(str(work / f'{n}.o') for n in ('kinds', 'sparse', 'sparse_coo')),
                        '-o', str(exe)], cwd=work, check=True)
        lines = [str(len(cases))]
        for c in cases:
            data = c['data']
            lines.append(f"{c['m']} {c['n']} {len(data)}")
            for r, col, val in zip(c['rows'], c['cols'], data):
                lines.append(f'{r+1} {col+1} {val:.17g}')
            if c['n']:
                lines.append(' '.join(str(v) for v in c['x']))
        result = subprocess.run([str(exe)], cwd=work, input='\n'.join(lines) + '\n',
                                text=True, capture_output=True, check=True, timeout=120)
        output = iter(result.stdout.splitlines())
        for idx, case in enumerate(cases, 1):
            m, n = case['m'], case['n']
            coo = sparse.coo_array((case['data'], (case['rows'], case['cols'])), shape=(m, n))
            canonical = coo.copy()
            canonical.sum_duplicates()
            f_rows, f_cols, f_vals = sparse.find(coo)
            expected_dense = coo.toarray()
            expected_y = np.asarray(coo.toarray() @ case['x'])
            header = next(output).split()
            assert header[0] == 'CASE' and int(header[1]) == idx, (mode,idx,header)
            got_m, got_n, got_orig, got_uniq, got_find = map(int,header[2:])
            assert (got_m, got_n, got_orig, got_uniq, got_find) == (
                m, n, coo.nnz, canonical.nnz, len(f_vals)), (mode,idx,header)
            for k in range(got_uniq):
                tag, row, col, val = next(output).split()
                assert tag == 'SUM'
                assert (int(row),int(col)) == (canonical.row[k]+1, canonical.col[k]+1), (
                    mode,idx,k,'SUM')
                np.testing.assert_allclose(float(val), canonical.data[k], rtol=0, atol=2e-12)
            for k in range(got_find):
                tag, row, col, val = next(output).split()
                assert tag == 'FIND'
                assert (int(row),int(col)) == (f_rows[k]+1, f_cols[k]+1), (
                    mode,idx,k,'FIND')
                np.testing.assert_allclose(float(val), f_vals[k], rtol=0, atol=2e-12)
            for i in range(m):
                for j in range(n):
                    tag, row, col, val = next(output).split()
                    assert (tag,int(row),int(col)) == ('D',i+1,j+1)
                    np.testing.assert_allclose(float(val), expected_dense[i,j], rtol=0, atol=2e-12)
                tag, row, val = next(output).split()
                assert (tag,int(row)) == ('Y',i+1)
                np.testing.assert_allclose(float(val), expected_y[i], rtol=0, atol=2e-12)
            tag, count = next(output).split()
            assert tag == 'NONZERO' and int(count) == np.count_nonzero(canonical.data), (
                mode,idx,'NONZERO')
        assert next(output,None) is None, f'unexpected trailing lines in {mode}'
        print(f'{mode}: {len(cases)} randomized SciPy COO/find cases: PASS', flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    p.add_argument('--trials', type=int, default=200)
    p.add_argument('--seed', type=int, default=20261105)
    args = p.parse_args()
    cases = generate_cases(args.trials,args.seed)
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        compare(mode,cases)


if __name__ == '__main__':
    main()
