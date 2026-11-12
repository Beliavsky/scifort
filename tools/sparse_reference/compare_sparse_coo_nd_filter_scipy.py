#!/usr/bin/env python3
"""Randomized COO N-D sparse mask comparisons with independent NumPy/SciPy references.

The Fortran program is built separately in checked and optimized modes and
accepts all cases in one stream, including dimensions with zero extent.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_array

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ['kinds', 'sparse', 'sparse_coo', 'sparse_coo_nd', 'sparse_coo_nd_filter']
OPS = ['>', '>=', '<', '<=', '==', '/=', 'abs>']


def compiled_driver(build: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    for name in SOURCES:
        subprocess.run([*cc, '-c', str(ROOT / 'src' / f'{name}.f90'),
                        '-o', str(build / f'{name}.o')], cwd=build, check=True)
    exe = build / 'coo_nd_filter_reference_driver'
    subprocess.run([*cc, str(ROOT / 'tools/sparse_reference/sparse_coo_nd_filter_reference_driver.f90'),
                    *(str(build / f'{name}.o') for name in SOURCES), '-o', str(exe)],
                   cwd=build, check=True)
    return exe


def random_array(rng: np.random.Generator, ndim: int, shape: tuple[int, ...], floating: bool):
    nnz = int(rng.integers(0, 48)) if all(shape) else 0
    coords = np.zeros((ndim, nnz), dtype=np.int64)
    for ax, extent in enumerate(shape):
        if nnz:
            coords[ax] = rng.integers(0, extent, size=nnz)
    if nnz >= 3 and rng.random() < .8:
        coords[:, -2:] = coords[:, :1]  # deliberately duplicate positions
    values = rng.integers(-8, 9, size=nnz).astype(float)
    if floating:
        values += rng.normal(0, .1, size=nnz)
        values[rng.random(nnz) < .2] = 0
    return coords, values


def aggregate(coords: np.ndarray, values: np.ndarray) -> dict[tuple[int, ...], float]:
    out: dict[tuple[int, ...], float] = {}
    for k in range(values.size):
        idx = tuple(int(x) for x in coords[:, k])
        out[idx] = out.get(idx, 0) + float(values[k])
    return out


def test_value(x: float, op: str, y: float) -> bool:
    return {
        '>': lambda: x > y, '>=': lambda: x >= y,
        '<': lambda: x < y, '<=': lambda: x <= y,
        '==': lambda: x == y, '/=': lambda: x != y,
        'abs>': lambda: abs(x) > y,
    }[op]()


def put_array(lines, shape, coords, values):
    lines.append(f'{len(shape)} {values.size}')
    lines.append(' '.join(str(x) for x in shape))
    for k in range(values.size):
        lines.append(' '.join(str(int(c)+1) for c in coords[:, k]) + f' {values[k]:.17g}')


def expected_output(shape, coords, values):
    return [(tuple(int(c)+1 for c in coords[:, k]), float(values[k]))
            for k in range(values.size)]


def make_trials(seed: int, trials: int, floating: bool):
    rng = np.random.default_rng(seed)
    lines = [str(trials)]
    expected = []
    for testid in range(trials):
        ndim = int(rng.integers(1, 7))
        shape = tuple(int(x) for x in rng.integers(0, 5, size=ndim))
        if rng.random() < .65:
            shape = tuple(max(1, x) for x in shape)
        coords, values = random_array(rng, ndim, shape, floating)
        # SciPy sparse COO is the independent source for dense aggregation.
        sparse = coo_array((values, tuple(coords)), shape=shape)
        dense = np.asarray(sparse.todense())
        agg = aggregate(coords, values)
        for key, value in agg.items():
            assert np.isclose(value, dense[key]), (key, value, dense[key])
        mode = int(rng.integers(1, 6))
        lines.append(str(mode))
        put_array(lines, shape, coords, values)
        if mode in (1, 4):
            op = OPS[int(rng.integers(len(OPS)))]
            bound = float(rng.choice([-4, -2, -1, 0, 1, 3, 5]))
            if floating and rng.random() < .5:
                bound += .17
            lines.append(f'"{op}" {bound:.17g}')
            if test_value(0., op, bound):
                expected.append(None)  # no sparse false-background representation
                continue
            keep = {key for key, value in agg.items() if test_value(value, op, bound)}
            if mode == 1:
                expected.append((shape, [(tuple(k+1 for k in key), 1.0)
                                         for key in sorted(keep)]))
            else:
                result = [(idx, v) for idx, v in expected_output(shape, coords, values)
                          if tuple(x-1 for x in idx) in keep]
                expected.append((shape, result))
        elif mode == 2:
            mc, mv = random_array(rng, ndim, shape, floating)
            put_array(lines, shape, mc, mv)
            inverse = bool(rng.integers(0, 2))
            lines.append(str(int(inverse)))
            magg = aggregate(mc, mv)
            take = {key for key, val in magg.items() if val != 0}
            result = [(idx, val) for idx, val in expected_output(shape, coords, values)
                      if ((tuple(x-1 for x in idx) in take) != inverse)]
            expected.append((shape, result))
        elif mode == 3:
            npoints = int(rng.integers(0, 50)) if all(shape) else 0
            points = np.zeros((ndim, npoints), dtype=int)
            for ax, extent in enumerate(shape):
                if npoints:
                    points[ax] = rng.integers(0, extent, size=npoints)
            inverse = bool(rng.integers(0, 2))
            lines.append(f'{npoints} {int(inverse)}')
            for k in range(npoints):
                lines.append(' '.join(str(p+1) for p in points[:, k]))
            set_points = {tuple(int(p) for p in points[:, k]) for k in range(npoints)}
            result = [(idx, val) for idx, val in expected_output(shape, coords, values)
                      if ((tuple(x-1 for x in idx) in set_points) != inverse)]
            expected.append((shape, result))
        else:
            expected.append(('count', int(np.count_nonzero(dense))))
    return '\n'.join(lines) + '\n', expected


def run_mode(mode: str, trials: int, seed: int, floating: bool):
    with tempfile.TemporaryDirectory(prefix='scifort-coo-filter-') as temp:
        build = Path(temp)
        exe = compiled_driver(build, mode)
        stream, expected = make_trials(seed, trials, floating)
        finished = subprocess.run([str(exe)], input=stream, text=True,
                                  capture_output=True, check=True)
        rows = finished.stdout.splitlines()
        assert len(rows) == trials, (len(rows), trials, finished.stderr[:2000])
        for i, (row, ref) in enumerate(zip(rows, expected)):
            tokens = row.split()
            status, nd = int(tokens[0]), int(tokens[1])
            if ref is None:
                assert status != 0 and nd == -1, (i, mode, row, ref)
                continue
            assert status == 0, (i, mode, row, ref)
            if ref[0] == 'count':
                assert nd == -2 and int(tokens[2]) == ref[1], (i, row, ref)
                continue
            shape, entries = ref
            n = int(tokens[2])
            assert nd == len(shape) and n == len(entries), (i, row, ref)
            assert tuple(int(x) for x in tokens[3:3+nd]) == shape, (i, row, ref)
            cursor = 3 + nd
            for idx, value in entries:
                assert tuple(int(x) for x in tokens[cursor:cursor+nd]) == idx, (i, row, ref)
                cursor += nd
                actual = float(tokens[cursor]); cursor += 1
                assert np.isclose(actual, value, rtol=1e-12, atol=1e-12), (i, actual, value)
            assert cursor == len(tokens), (i, row, ref)
        print(f'{mode}: PASS {trials} randomized COO predicates/masks '
              f'(seed={seed}, floating={floating})')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['checked','optimized','both'], default='both')
    parser.add_argument('--trials', type=int, default=500)
    parser.add_argument('--seed', type=int, default=20261105)
    parser.add_argument('--floating', action='store_true')
    args = parser.parse_args()
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        run_mode(mode, args.trials, args.seed, args.floating)


if __name__ == '__main__':
    main()
