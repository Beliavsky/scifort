#!/usr/bin/env python3
"""Compare Fortran COO N-D masks/scatter against NumPy and SciPy 1.17 dense references.

Cases preserve unsorted duplicate COO entries and explicit zeros, verify stored
entry counts (including mask-preserved duplicates), and apply repeated updates.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_array

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ['kinds', 'sparse', 'sparse_coo', 'sparse_coo_nd',
           'sparse_coo_nd_index', 'sparse_coo_nd_edit']


def compile_driver(build: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    for name in SOURCES:
        subprocess.run([*cc, '-c', str(ROOT / 'src' / f'{name}.f90'),
                        '-o', str(build / f'{name}.o')], cwd=build, check=True)
    exe = build / 'sparse_coo_nd_edit_ref'
    subprocess.run([*cc, str(ROOT / 'tools/sparse_reference/sparse_coo_nd_edit_reference_driver.f90'),
                    *(str(build / f'{name}.o') for name in SOURCES), '-o', str(exe)],
                   cwd=build, check=True)
    return exe


def base_case(rng: np.random.Generator, floating: bool):
    ndim = int(rng.integers(1, 6))
    shape = tuple(int(x) for x in rng.integers(0, 6, size=ndim))
    if rng.random() < .6:
        shape = tuple(max(1, x) for x in shape)
    n = 0 if not all(shape) else int(rng.integers(0, 40))
    coords = np.empty((ndim, n), dtype=int)
    for ax, sz in enumerate(shape):
        coords[ax, :] = rng.integers(0, sz, size=n)
    if n >= 3 and rng.random() < .7:
        coords[:, -2:] = coords[:, :1]
    values = rng.integers(-9, 10, size=n).astype(float)
    if floating:
        values += rng.normal(scale=.05, size=n)
        values[rng.random(n) < .2] = 0
    dense = np.zeros(shape, dtype=float)
    if n:
        np.add.at(dense, tuple(coords), values)
    d = {}
    for j in range(n):
        key = tuple(int(x) for x in coords[:, j])
        d[key] = d.get(key, 0.0) + float(values[j])
    return shape, coords, values, dense, d


def selected_mask(rng: np.random.Generator, n: int) -> np.ndarray:
    return rng.random(n) < rng.uniform(.0, 1.)


def events_for(rng: np.random.Generator, shape: tuple[int, ...], floating: bool):
    n = 0 if not all(shape) else int(rng.integers(0, 30))
    coords = np.empty((len(shape), n), dtype=int)
    for ax, s in enumerate(shape):
        coords[ax, :] = rng.integers(0, s, size=n)
    if n >= 3 and rng.random() < .7:
        coords[:, -3:] = coords[:, :1]
    values = rng.integers(-8, 9, size=n).astype(float)
    if floating:
        values += rng.normal(scale=.05, size=n)
        values[rng.random(n) < .2] = 0.0
    return coords, values


def updated(reference: dict, keys: list[tuple[int, ...]], values, add: bool, keep: bool):
    result = reference.copy()
    touched = set()
    for key, value in zip(keys, values):
        touched.add(key)
        result[key] = result.get(key, 0.0) + float(value) if add else float(value)
    if not keep:
        for key in touched:
            if result[key] == 0.0:
                result.pop(key)
    return result


def dense_of_dict(shape: tuple[int, ...], data: dict) -> np.ndarray:
    out = np.zeros(shape, dtype=float)
    for key, value in data.items():
        out[key] = value
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    ap.add_argument('--trials', type=int, default=650)
    ap.add_argument('--seed', type=int, default=20261105)
    ap.add_argument('--floating', action='store_true')
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    inputs = [str(args.trials)]
    expected = []
    scipy_checks = 0
    scipy_skips = 0
    operation_counts = {i: 0 for i in range(1, 6)}

    for _ in range(args.trials):
        shape, coords, values, dense, store = base_case(rng, args.floating)
        operation = int(rng.integers(1, 6))
        operation_counts[operation] += 1
        ndim = len(shape)
        inputs.extend([str(operation), f'{ndim} {len(values)}',
                       ' '.join(map(str, shape))])
        for j in range(len(values)):
            inputs.append(' '.join(map(str, coords[:, j]+1)) + f' {values[j]:.17g}')
        try:
            sp = coo_array((values, tuple(coords)), shape=shape)
            np.testing.assert_allclose(sp.toarray(), dense, atol=1e-11, rtol=1e-11)
            scipy_checks += 1
        except (TypeError, ValueError, NotImplementedError):
            scipy_skips += 1

        if operation == 1:
            mask = selected_mask(rng, dense.size)
            inputs.append(str(len(mask)))
            if mask.size:
                inputs.append(' '.join('T' if x else 'F' for x in mask))
            wanted = np.asarray(dense.ravel(order='C')[mask], dtype=float)
            nstored = 0
            for j in range(len(values)):
                index = np.ravel_multi_index(tuple(coords[:, j]), shape, order='C')
                nstored += int(mask[index])
        elif operation == 2:
            ax = int(rng.integers(ndim))
            mask = selected_mask(rng, shape[ax])
            inputs.append(f'{ax+1} {len(mask)}')
            if mask.size:
                inputs.append(' '.join('T' if x else 'F' for x in mask))
            wanted = np.compress(mask, dense, axis=ax)
            nstored = sum(bool(mask[coords[ax, j]]) for j in range(len(values)))
        elif operation in (3, 4):
            points, updates = events_for(rng, shape, args.floating)
            keep_zero = bool(rng.integers(0, 2))
            inputs.append(f'{len(updates)} {int(keep_zero)} 1000000')
            for j in range(len(updates)):
                inputs.append(' '.join(map(str, points[:, j]+1)) + f' {updates[j]:.17g}')
            keys = [tuple(int(x) for x in points[:, j]) for j in range(len(updates))]
            if len(updates):
                result = updated(store, keys, updates, add=(operation==4), keep=keep_zero)
            else:
                # No edits is a true no-op; input storage keeps duplicates.
                result = None
            wanted = dense_of_dict(shape, result) if result is not None else dense
            nstored = len(result) if result is not None else len(values)
        else:
            mask = selected_mask(rng, dense.size)
            add = bool(rng.integers(0, 2))
            keep_zero = bool(rng.integers(0, 2))
            scalar = float(rng.integers(-5, 6))
            if args.floating:
                scalar += float(rng.normal(scale=.05))
            inputs.append(f'{len(mask)} {int(add)} {int(keep_zero)} {scalar:.17g} 1000000')
            if mask.size:
                inputs.append(' '.join('T' if x else 'F' for x in mask))
            keys = [tuple(int(x) for x in np.unravel_index(i, shape, order='C'))
                    for i in np.flatnonzero(mask)]
            if keys:
                result = updated(store, keys, [scalar]*len(keys), add=add, keep=keep_zero)
            else:
                result = None
            wanted = dense_of_dict(shape, result) if result is not None else dense
            nstored = len(result) if result is not None else len(values)
        expected.append((wanted, nstored, operation))

    stream = '\n'.join(inputs)+'\n'
    for mode in (('checked', 'optimized') if args.mode == 'both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-coo-nd-edit-') as tmp:
            exe = compile_driver(Path(tmp), mode)
            proc = subprocess.run([str(exe)], input=stream, text=True,
                                  capture_output=True)
            if proc.returncode:
                raise RuntimeError(f"driver exited {proc.returncode}: {proc.stderr}\nstdout tail: {proc.stdout[-500:]}")
            lines = proc.stdout.strip().splitlines()
            if len(lines) != args.trials:
                raise AssertionError(f'{mode} returned {len(lines)} lines, expected {args.trials}: {proc.stderr}')
            for j, (line, (want, nnz, op)) in enumerate(zip(lines, expected)):
                fields = line.split()
                status, ndim, actual_nnz = map(int, fields[:3])
                if status:
                    raise AssertionError(f'{mode} case {j} (operation {op}): status {status}')
                sizes = tuple(map(int, fields[3:3+ndim]))
                got = np.array([float(x) for x in fields[3+ndim:]], dtype=float)
                if sizes != want.shape or actual_nnz != nnz:
                    raise AssertionError(f'{mode} case {j} op {op}: shape {sizes} vs {want.shape}, '
                                         f'nnz {actual_nnz} vs {nnz}')
                got = got.reshape(sizes, order='C')
                try:
                    np.testing.assert_allclose(got,want,rtol=3e-11,atol=3e-11)
                except AssertionError as err:
                    raise AssertionError(f'{mode} case {j} op {op}: {err}') from err
            print(f'{mode}: PASS {args.trials} COO N-D mask/scatter cases; '
                  f'operations {operation_counts}; SciPy source checks {scipy_checks}, '
                  f'unsupported {scipy_skips}',flush=True)


if __name__ == '__main__':
    main()
