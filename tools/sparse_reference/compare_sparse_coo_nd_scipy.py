#!/usr/bin/env python3
"""Independent randomized SciPy comparison for SciFort N-dimensional COO.

Comparison covers dense values, canonical sparse entries, C/F reshaping,
axis permutations, scalar lookup, and 2-D COO conversion. The Fortran driver
receives unsorted one-based coordinates, while SciPy uses zero-based axes.
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
        ndim = int(rng.integers(1, 8))
        shape = tuple(int(v) for v in rng.integers(1, 6, size=ndim))
        if np.prod(shape) > 200:
            shape = tuple(min(v, 2) for v in shape)
        if i % 17 == 0:
            shape = shape[:i % ndim] + (0,) + shape[i % ndim + 1:]
        count = 0 if 0 in shape else int(rng.integers(0, 101))
        coords = np.zeros((ndim, count), dtype=int)
        for ax, extent in enumerate(shape):
            if extent > 0:
                coords[ax] = rng.integers(0, extent, size=count)
        if count > 0:
            # Explicit duplicates, zero entries and cancelling pairs.
            for j in range(1, min(count, 10), 3):
                coords[:, j] = coords[:, 0]
        vals = rng.integers(-8, 9, size=count).astype(float)
        if count > 3:
            vals[1] = 0
            vals[2] = -vals[0]
        perm = tuple(int(v) for v in rng.permutation(ndim))
        total = int(np.prod(shape))
        if i % 3 == 0:
            target = (total,)
        elif i % 3 == 1:
            target = (1, total, 1)
        else:
            target = (2, total // 2) if total % 2 == 0 else (total, 1)
        layout = 'F' if i % 2 else 'C'
        cases.append(dict(shape=shape, coords=coords, data=vals,
                          perm=perm, target=target, layout=layout))
    return cases


def compare(mode: str, cases: list[dict]) -> None:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-coond-') as td:
        work = Path(td)
        cmd = ['gfortran', '-std=f2018', *flags, '-J', str(work), '-I', str(work)]
        for name in ('kinds', 'sparse', 'sparse_coo', 'sparse_coo_nd'):
            subprocess.run([*cmd, '-c', str(ROOT / 'src' / f'{name}.f90'),
                            '-o', str(work / f'{name}.o')], cwd=work, check=True)
        exe = work / 'coond_driver'
        subprocess.run([*cmd, str(ROOT / 'tools/sparse_reference/sparse_coo_nd_reference_driver.f90'),
                        *(str(work / f'{n}.o') for n in ('kinds','sparse','sparse_coo','sparse_coo_nd')),
                        '-o', str(exe)], cwd=work, check=True)
        lines = [str(len(cases))]
        for c in cases:
            shape, coords, target, layout = c['shape'], c['coords'], c['target'], c['layout']
            lines.append(f"{len(shape)} {len(c['data'])} {len(target)} {int(layout == 'F')}")
            lines.append(' '.join(map(str, shape)))
            lines.append(' '.join(map(str, target)))
            lines.append(' '.join(str(p+1) for p in c['perm']))
            for entry in range(len(c['data'])):
                lines.append(' '.join(str(v+1) for v in coords[:,entry]) +
                             f" {c['data'][entry]:.17g}")
        result = subprocess.run([str(exe)], cwd=work, input='\n'.join(lines)+'\n',
                                text=True, capture_output=True, check=True, timeout=120)
        iterator = iter(result.stdout.splitlines())
        for i, c in enumerate(cases, 1):
            header = next(iterator).split()
            assert header[0] == 'CASE' and int(header[1]) == i, (mode, i, header)
            shape, coords, vals = c['shape'], c['coords'], c['data']
            sp = sparse.coo_array((vals, tuple(coords[ax] for ax in range(len(shape)))), shape=shape)
            assert int(header[2]) == len(shape) and int(header[3]) == len(vals)

            def dense(label: str, expected: np.ndarray) -> None:
                row = next(iterator).split()
                assert row[0] == label, (mode, i, label, row)
                count = int(row[1])
                expected = np.asarray(expected).ravel(order='C')
                assert count == len(expected), (mode, i, label, count, len(expected))
                obtained = np.fromiter((float(next(iterator)) for _ in range(count)), float, count=count)
                np.testing.assert_array_equal(obtained, expected, err_msg=f'{mode} case {i} {label}')

            dense('BASE', sp.toarray())
            dense('TRANSPOSE', sp.transpose(c['perm']).toarray())
            dense('RESHAPE', sp.reshape(c['target'], order=c['layout']).toarray())
            canon = sp.copy()
            canon.sum_duplicates()
            row = next(iterator).split()
            assert row[0] == 'CANON', (mode, i, row)
            count = int(row[1])
            assert count == canon.nnz, (mode, i, count, canon.nnz)
            got_coords = np.empty((len(shape), count), dtype=int)
            got_data = np.empty(count, dtype=float)
            for j in range(count):
                parts = next(iterator).split()
                got_coords[:, j] = [int(v)-1 for v in parts[:len(shape)]]
                got_data[j] = float(parts[len(shape)])
            np.testing.assert_array_equal(got_coords,
                                          np.asarray(canon.coords),
                                          err_msg=f'{mode} case {i} CANON coords')
            np.testing.assert_array_equal(got_data, canon.data,
                                          err_msg=f'{mode} case {i} CANON values')
            dense('FROM_FLAT', sp.toarray())
            row = next(iterator).split()
            assert row[0] == 'ROUNDTRIP' and row[1] == '1', (mode, i, row)
            row = next(iterator).split()
            assert row[0] == 'GET', (mode, i, row)
            expected = float(sp.toarray()[tuple(coords[:,0])]) if len(vals) else 0.0
            assert float(row[1]) == expected, (mode, i, 'GET', row[1], expected)
        assert next(iterator, None) is None, (mode, 'extra output')
    print(f'{mode}: {len(cases)} SciPy N-D COO cases PASS', flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    p.add_argument('--trials', type=int, default=250)
    p.add_argument('--seed', type=int, default=20261105)
    args = p.parse_args()
    cases = generate_cases(args.trials,args.seed)
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        compare(mode,cases)


if __name__ == '__main__':
    main()
