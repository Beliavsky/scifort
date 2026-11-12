#!/usr/bin/env python3
"""Randomized complement-aware N-D Boolean-mask tests against NumPy/SciPy.

Fortran outputs exact default bits and canonical exception sets; an independent
NumPy truth table is used to check every coordinate (including implicit ones).
Explicit COO zeros, duplicates, uncoalesced order and zero-size dimensions are
covered. These mask operations are a Fortran extension, not an existing SciPy API.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_array

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('kinds', 'sparse', 'sparse_coo', 'sparse_coo_nd', 'sparse_coo_nd_bool')
OPS = ('>', '>=', '<', '<=', '==', '/=', 'abs>')


def build_driver(build: Path, mode: str) -> Path:
    flags = ('-O0', '-g', '-fcheck=all', '-fbacktrace') if mode == 'checked' else ('-O3',)
    cc = ('gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build))
    for name in SOURCES:
        subprocess.run((*cc, '-c', str(ROOT / 'src' / f'{name}.f90'),
                        '-o', str(build / f'{name}.o')), cwd=build, check=True)
    exe = build / 'sparse_coo_nd_bool_reference_driver'
    subprocess.run((*cc, str(ROOT / 'tools/sparse_reference' /
                                  'sparse_coo_nd_bool_reference_driver.f90'),
                    *(str(build / f'{name}.o') for name in SOURCES), '-o', str(exe)),
                   cwd=build, check=True)
    return exe


def point_coords(rng: np.random.Generator, shape: tuple[int, ...], limit: int = 48):
    n = int(rng.integers(0, limit + 1)) if all(shape) else 0
    coords = np.zeros((len(shape), n), dtype=np.int64)
    for axis, dim in enumerate(shape):
        if n:
            coords[axis] = rng.integers(0, dim, size=n)
    if n > 4 and rng.random() < .85:
        coords[:, -3:] = coords[:, :1]
    return coords


def real_coo(rng, shape):
    coords = point_coords(rng, shape)
    values = rng.integers(-6, 7, size=coords.shape[1]).astype(float)
    if values.size >= 3 and rng.random() < .8:
        values[-2:] = (2.0, -2.0)
        coords[:, -2:] = coords[:, :1]
    if values.size:
        values[rng.random(values.size) < .25] = 0.0
    return coords, values


def put_point_arrays(lines, points):
    for pt in points.T:
        lines.append(' '.join(str(int(x) + 1) for x in pt))


def put_coo(lines, shape, coords, values):
    lines.append(f'{len(shape)} {len(values)}')
    lines.append(' '.join(str(i) for i in shape))
    for point, value in zip(coords.T, values):
        lines.append(' '.join(str(int(i)+1) for i in point) + f' {value:.17g}')


def make_mask(shape, coords, default):
    x = np.full(shape, bool(default), dtype=bool)
    for pt in coords.T:
        x[tuple(int(i) for i in pt)] = not bool(default)
    return x


def aggregate_scipy(shape, coords, values):
    x = coo_array((values, tuple(coords)), shape=shape)
    return np.asarray(x.todense())


def compare(dense, op, threshold):
    if op == '>': return dense > threshold
    if op == '>=': return dense >= threshold
    if op == '<': return dense < threshold
    if op == '<=': return dense <= threshold
    if op == '==': return dense == threshold
    if op == '/=': return dense != threshold
    if op == 'abs>': return abs(dense) > threshold
    raise ValueError(op)


def make_cases(seed: int, trials: int):
    rng = np.random.default_rng(seed)
    lines = [str(trials)]
    expected = []
    modes = np.arange(1, 11)
    for t in range(trials):
        mode = int(modes[t % len(modes)] if t < 10 else rng.choice(modes))
        rank = int(rng.integers(1, 6))
        shape = tuple(int(x) for x in rng.integers(0, 5, size=rank))
        if rng.random() < .72:
            shape = tuple(max(1, n) for n in shape)
        a = point_coords(rng, shape)
        b = point_coords(rng, shape)
        da, db = bool(rng.integers(2)), bool(rng.integers(2))
        lines += [f'{mode} {rank}', ' '.join(str(s) for s in shape),
                  f'{a.shape[1]} {b.shape[1]} {int(da)} {int(db)}']
        put_point_arrays(lines, a)
        put_point_arrays(lines, b)
        left = make_mask(shape, a, da)
        right = make_mask(shape, b, db)
        if mode <= 5:
            out = {1: np.logical_and, 2: np.logical_or,
                   3: np.logical_xor, 4: lambda x, y: x & ~y}.get(mode)
            dense = ~left if mode == 5 else out(left, right)
            expected.append(('mask', shape, dense, bool(da ^ True) if mode == 5 else
                             {1: da and db, 2: da or db, 3: da != db,
                              4: da and not db}[mode]))
        elif mode == 6:
            coords, values = real_coo(rng, shape)
            put_coo(lines, shape, coords, values)
            op = str(rng.choice(OPS))
            threshold = float(rng.choice([-3, -1, 0, 1, 4]))
            lines.append(f'"{op}" {threshold:.17g}')
            dense = compare(aggregate_scipy(shape, coords, values), op, threshold)
            background = bool(compare(np.array([0.0]), op, threshold)[0])
            expected.append(('mask', shape, dense, background))
        elif mode == 7:
            coords, values = real_coo(rng, shape)
            put_coo(lines, shape, coords, values)
            keep = [k for k, pt in enumerate(coords.T)
                    if left[tuple(int(i) for i in pt)]]
            expected.append(('coo', shape, [(tuple(int(p)+1 for p in coords[:, k]),
                                             float(values[k])) for k in keep]))
        elif mode == 8:
            expected.append(('count', int(np.count_nonzero(left))))
        elif mode == 9:
            if da and all(shape):
                expected.append(('reject',))
            else:
                selected = np.argwhere(left)
                expected.append(('coo', shape, [(tuple(int(p)+1 for p in pt), 1.0)
                                                 for pt in selected]))
        else:
            coords, values = real_coo(rng, shape)
            put_coo(lines, shape, coords, values)
            dense = aggregate_scipy(shape, coords, values) != 0
            if db:
                dense = ~dense
            expected.append(('mask', shape, dense, db))
    return '\n'.join(lines) + '\n', expected


def verify_output(row: str, ref, mode: str, test_no: int):
    tokens = row.split()
    info, ndim = int(tokens[0]), int(tokens[1])
    if ref[0] == 'reject':
        assert info != 0 and ndim == -1, (mode, test_no, row, ref)
        return
    assert info == 0, (mode, test_no, row, ref)
    if ref[0] == 'count':
        assert ndim == -2 and int(tokens[2]) == ref[1], (mode, test_no, row, ref)
        return
    _, shape = ref[:2]
    assert ndim == len(shape), (mode, test_no, row, ref)
    nexceptions = int(tokens[2])
    if ref[0] == 'mask':
        _, _, dense, default = ref
        assert int(tokens[3]) == int(default), (mode, test_no, row, ref)
        got_shape = tuple(int(x) for x in tokens[4:4+ndim])
        assert got_shape == shape, (mode, test_no, row, ref)
        locations = [tuple(int(x)-1 for x in tokens[k:k+ndim])
                     for k in range(4+ndim,len(tokens),ndim)]
        assert len(locations) == nexceptions, (mode, test_no, row, ref)
        assert locations == sorted(set(locations)), (mode, test_no, row, ref)
        expected_locs = [tuple(int(x) for x in pt)
                         for pt in np.argwhere(dense != default)]
        assert locations == expected_locs, (mode, test_no, row, locations, expected_locs)
    else:
        entries = ref[2]
        assert nexceptions == len(entries), (mode, test_no, row, ref)
        assert tuple(map(int,tokens[3:3+ndim])) == shape, (mode, test_no, row, ref)
        cursor = 3+ndim
        for idx, val in entries:
            actual_idx = tuple(int(x) for x in tokens[cursor:cursor+ndim])
            cursor += ndim
            got = float(tokens[cursor]); cursor += 1
            assert actual_idx == idx and np.isclose(got,val,rtol=0,atol=1e-12), (
                mode,test_no,actual_idx,idx,got,val)
        assert cursor == len(tokens), (mode,test_no,row,ref)


def run_one(mode, trials, seed):
    with tempfile.TemporaryDirectory(prefix='scifort-sparse-bool-') as tmp:
        exe = build_driver(Path(tmp),mode)
        stream, expected = make_cases(seed,trials)
        completed = subprocess.run([str(exe)],input=stream,text=True,
                                   capture_output=True,check=True)
        rows = completed.stdout.splitlines()
        assert len(rows) == trials, (mode,trials,len(rows),completed.stderr[:3000])
        for i, (row, ref) in enumerate(zip(rows,expected)):
            verify_output(row,ref,mode,i)
        print(f'{mode}: PASS {trials} randomized Boolean-mask cases (seed={seed})', flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    p.add_argument('--trials',type=int,default=500)
    p.add_argument('--seed',type=int,default=20261105)
    args = p.parse_args()
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        run_one(mode,args.trials,args.seed)


if __name__ == '__main__':
    main()
