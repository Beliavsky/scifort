#!/usr/bin/env python3
"""Randomized COO N-D indexing/slicing comparisons against SciPy and NumPy.

SciPy 1.17 coo_array constructs the N-D source and implements basic indexing;
NumPy independently supplies the reference for orthogonal multi-axis index
arrays, whose Cartesian semantics differ from simultaneous advanced indexing.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_array

ROOT = Path(__file__).resolve().parents[2]


def example(rng: np.random.Generator, floating: bool) -> tuple[tuple[int, ...], np.ndarray, np.ndarray, np.ndarray]:
    ndim = int(rng.integers(1, 6))
    shape = tuple(int(s) for s in rng.choice([0, 1, 1, 2, 3, 4], size=ndim))
    n = int(rng.integers(0, 25)) if all(shape) else 0
    coords = np.empty((ndim, n), dtype=np.int64)
    for i, s in enumerate(shape):
        coords[i] = rng.integers(0, s, size=n)
    if n > 2 and rng.random() < .5:
        coords[:, -1] = coords[:, 0]
    if floating:
        values = rng.normal(0, 2, size=n)
    else:
        values = rng.integers(-5, 6, size=n).astype(float)
    if n and rng.random() < .25:
        values[0] = 0
    # COO constructor intentionally includes stored zeros and duplicate positions.
    sp = coo_array((values, tuple(coords[i] for i in range(ndim))), shape=shape)
    return shape, coords, values, sp.toarray()


def selectors_case(rng: np.random.Generator, shape: tuple[int, ...], dense: np.ndarray):
    lines: list[str] = []
    result = dense
    removed = 0
    scipy_keys = []
    n_arrays = 0
    keepdims = bool(rng.integers(0, 2))
    for ax, extent in enumerate(shape):
        mode = int(rng.choice([0, 1, 2, 3])) if extent else int(rng.choice([0, 2, 3]))
        if mode == 0:
            key = np.arange(extent)
            lines.append('0 1 1 1 0 0 0')
            scipy_key = slice(None)
        elif mode == 1:
            v = int(rng.integers(1, extent+1))
            key = [v-1] if keepdims else v-1
            scipy_key = slice(v-1, v) if keepdims else v-1
            lines.append(f'1 {v} 1 1 1 1 0')
            if not keepdims:
                removed += 1
        elif mode == 2:
            step = int(rng.choice([-3, -2, -1, 1, 2, 3]))
            hf = bool(rng.integers(0, 2)) if extent else False
            hl = bool(rng.integers(0, 2)) if extent else False
            first = int(rng.integers(1, extent+1)) if hf else (extent if step < 0 else 1)
            last = int(rng.integers(1, extent+1)) if hl else (1 if step < 0 else extent)
            key = np.asarray(list(range(first-1, last if step > 0 else last-2, step)), dtype=int)
            if extent == 0:
                key = np.empty(0, dtype=int)
            lines.append(f'2 {first} {last} {step} {int(hf)} {int(hl)} 0')
            # slice(None,None,-1) correctly handles ending at the first element.
            start = first-1 if hf else None
            stop = (last if step > 0 else (last-2 if last > 1 else None)) if hl else None
            scipy_key = slice(start, stop, step)
        else:
            n_arrays += 1
            n = int(rng.integers(0, 6)) if extent else 0
            values = rng.integers(1, extent+1, size=n).tolist() if extent else []
            key = np.asarray(values, dtype=int) - 1
            lines.append(f'3 1 1 1 0 0 {n}')
            if n:
                lines.append(' '.join(map(str, values)))
            scipy_key = key
        if mode != 0:
            result = np.take(result, key, axis=ax-(removed if mode != 1 or keepdims else removed-1))
        scipy_keys.append(scipy_key)
    # Rank-zero representation is a rank-one length-one array in SciFort.
    if np.ndim(result) == 0:
        result = np.asarray(result).reshape(1)
    return [str(int(keepdims)), *lines], np.asarray(result), scipy_keys, n_arrays


def gather_case(rng: np.random.Generator, shape: tuple[int, ...], dense: np.ndarray):
    n = int(rng.integers(0, 20)) if all(shape) else 0
    coords = np.empty((len(shape), n), dtype=int)
    for i, s in enumerate(shape):
        coords[i] = rng.integers(0, s, size=n)
    if n and rng.random() < .4:
        coords[:, -1] = coords[:, 0]
    output = np.array([dense[tuple(coords[:, j])] for j in range(n)], dtype=float)
    lines = [str(n)]
    lines += [' '.join(map(str, coords[:, j]+1)) for j in range(n)]
    return lines, output


def build_driver(build: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    names = ['kinds', 'sparse', 'sparse_coo', 'sparse_coo_nd', 'sparse_coo_nd_index']
    for name in names:
        subprocess.run([*cc, '-c', str(ROOT / 'src' / (name+'.f90')),
                        '-o', str(build / (name+'.o'))], cwd=build, check=True)
    exe = build / 'sparse_index_ref'
    subprocess.run([*cc, str(ROOT / 'tools/sparse_reference/sparse_coo_nd_index_reference_driver.f90'),
                    *(str(build / (name+'.o')) for name in names), '-o', str(exe)],
                   cwd=build, check=True)
    return exe


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    ap.add_argument('--trials', type=int, default=500)
    ap.add_argument('--seed', type=int, default=20261105)
    ap.add_argument('--floating', action='store_true')
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    lines = [str(args.trials)]
    expected = []
    scipy_checks = 0
    scipy_skips = 0
    for trial in range(args.trials):
        shape, coords, values, dense = example(rng, args.floating)
        mode = 1 if rng.random() < .75 else 2
        lines.extend([str(mode), f'{len(shape)} {len(values)}', ' '.join(map(str, shape))])
        for j, v in enumerate(values):
            lines.append(' '.join(str(x+1) for x in coords[:, j]) + f' {v:.17g}')
        if mode == 1:
            addition, result, scipy_keys, n_arrays = selectors_case(rng, shape, dense)
            lines.extend(addition)
            # Check SciPy basic indexing. Advanced indices use orthogonal semantics here.
            if all(shape) and n_arrays == 0:
                try:
                    sp = coo_array((values, tuple(coords[i] for i in range(len(shape)))), shape=shape)
                    retrieved = sp[tuple(scipy_keys)]
                    retrieved = retrieved.toarray() if hasattr(retrieved, 'toarray') else np.asarray(retrieved)
                    if np.asarray(retrieved).ndim == 0:
                        retrieved = np.asarray(retrieved).reshape(1)
                    np.testing.assert_allclose(result, retrieved, atol=5e-13, rtol=5e-13)
                    scipy_checks += 1
                except (TypeError, NotImplementedError, IndexError, ValueError):
                    # Some valid mixed-rank N-D index patterns are not handled by SciPy 1.17.
                    scipy_skips += 1
            expected.append(result)
        else:
            addition, result = gather_case(rng, shape, dense)
            lines.extend(addition)
            expected.append(result)
    data = '\n'.join(lines) + '\n'
    for mode in (('checked', 'optimized') if args.mode == 'both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix='scifort-coo-nd-index-') as temp:
            exe = build_driver(Path(temp), mode)
            out = subprocess.run([str(exe)], input=data, text=True, capture_output=True, check=True)
            results = out.stdout.strip().splitlines()
            if len(results) != len(expected):
                raise AssertionError(f'{mode}: expected {len(expected)} output lines, got {len(results)}; stderr={out.stderr}')
            for case_id, (line, wanted) in enumerate(zip(results, expected)):
                fields = line.split()
                status, rank = map(int, fields[:2])
                if status != 0:
                    raise AssertionError(f'{mode} case {case_id}: invalid status {status}')
                shape = tuple(map(int, fields[2:2+rank]))
                observed = np.asarray([float(x) for x in fields[2+rank:]], dtype=float)
                if shape != wanted.shape:
                    raise AssertionError(f'{mode} case {case_id}: shape {shape} vs {wanted.shape}')
                observed = observed.reshape(shape, order='C')
                try:
                    np.testing.assert_allclose(observed, wanted, atol=5e-12, rtol=5e-12)
                except AssertionError as exc:
                    raise AssertionError(f'{mode} case {case_id}: {exc}') from exc
            print(f'{mode}: PASS {args.trials} randomized N-D indexing/gather comparisons ' 
                  f'(SciPy cross-checks {scipy_checks}, unsupported SciPy patterns {scipy_skips})', flush=True)


if __name__ == '__main__':
    main()
