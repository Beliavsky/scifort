#!/usr/bin/env python3
"""Randomized NumPy checks for complement-aware COO Boolean broadcasting/reductions.

The mask operations extend SciFort beyond SciPy's public sparse-Boolean API.
Reference truth values are evaluated using independent NumPy dense operations
on deliberately small shapes. Stress/overflow cases are native Fortran tests.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('kinds', 'sparse', 'sparse_coo', 'sparse_coo_nd',
           'sparse_coo_nd_bool', 'sparse_coo_nd_bool_reduce')


def compile_driver(build: Path, mode: str) -> Path:
    flags = ('-O0', '-g', '-fcheck=all', '-fbacktrace') if mode == 'checked' else ('-O3',)
    cc = ('gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build))
    for source in SOURCES:
        subprocess.run((*cc, '-c', str(ROOT / 'src' / (source + '.f90')),
                        '-o', str(build / (source + '.o'))), cwd=build, check=True)
    exe = build / 'bool_reduce_reference_driver'
    subprocess.run((*cc, str(ROOT / 'tools/sparse_reference' /
                                  'sparse_coo_nd_bool_reduce_reference_driver.f90'),
                    *(str(build / (source + '.o')) for source in SOURCES), '-o', str(exe)),
                   cwd=build, check=True)
    return exe


def mask_points(rng: np.random.Generator, shape: tuple[int, ...]):
    n = int(rng.integers(0, 14)) if all(shape) else 0
    if n:
        points = np.vstack([rng.integers(1, d + 1, size=n) for d in shape]).T
    else:
        points = np.empty((0, len(shape)), dtype=np.int64)
    return points


def mask_dense(shape: tuple[int, ...], points: np.ndarray, default: bool) -> np.ndarray:
    a = np.full(shape, default, dtype=bool)
    for point in points:
        a[tuple(int(v) - 1 for v in point)] = not default
    return a


def append_mask(lines: list[str], shape: tuple[int, ...], points: np.ndarray, default: bool):
    lines.append(' '.join(str(d) for d in shape))
    lines.append(f'{len(points)} {int(default)}')
    lines.extend(' '.join(map(str, point)) for point in points)


def ensure_one_or_equal(rng, target):
    return tuple(int(v if rng.random() < 0.55 else 1) for v in target)


def make_trials(rng: np.random.Generator, trials: int):
    lines = [str(trials)]
    expected = []
    for k in range(trials):
        mode = 1 + k % 6 if k < 6 else int(rng.integers(1, 7))
        rank = int(rng.integers(1, 6))
        if mode == 5:
            target = tuple(map(int, rng.integers(0, 5, size=rank)))
            if rng.random() < .85:
                target = tuple(max(1, v) for v in target)
            ashape = ensure_one_or_equal(rng, target)
            bshape = ensure_one_or_equal(rng, target)
            # Force at least one operand to carry each zero (otherwise 1,1 -> 1).
            ashape = tuple(0 if v == 0 and rng.random() < .5 else u
                           for v, u in zip(target, ashape))
            bshape = tuple(0 if v == 0 and u != 0 else u
                           for v, u in zip(target, bshape))
            if rng.random() < .5:
                ashape = ashape[int(rng.integers(rank)):]  # right-aligned lower rank
            if rng.random() < .5:
                bshape = bshape[int(rng.integers(rank)):]
            # Broadcasting after dropping leading axes may turn a desired 0
            # extent into 1; the reference uses actual broadcast_shapes.
            points_a, points_b = mask_points(rng, ashape), mask_points(rng, bshape)
            da, db = bool(rng.integers(2)), bool(rng.integers(2))
            a, b = mask_dense(ashape,points_a,da), mask_dense(bshape,points_b,db)
            shape = np.broadcast_shapes(ashape,bshape)
            mode_op = int(rng.integers(1, 5))
            if mode_op == 1: ref = np.logical_and(a,b)
            elif mode_op == 2: ref = np.logical_or(a,b)
            elif mode_op == 3: ref = np.logical_xor(a,b)
            else: ref = np.logical_and(a,np.logical_not(b))
            lines.append(f'5 {len(ashape)}')
            append_mask(lines,ashape,points_a,da)
            lines.append(str(len(bshape)))
            append_mask(lines,bshape,points_b,db)
            lines.append(str(mode_op))
        else:
            shape = tuple(map(int, rng.integers(0, 5, size=rank)))
            if rng.random() < .83:
                shape = tuple(max(1, v) for v in shape)
            points = mask_points(rng,shape)
            default = bool(rng.integers(2))
            dense = mask_dense(shape,points,default)
            lines.append(f'{mode} {rank}')
            append_mask(lines,shape,points,default)
            if mode == 1:
                expand = list(shape)
                for i, d in enumerate(expand):
                    if d == 1 and rng.random() < .7:
                        expand[i] = int(rng.integers(0,5))
                if rng.random() < .65:
                    expand = list(map(int,rng.integers(0,4,size=int(rng.integers(1,3))))) + expand
                target = tuple(expand)
                lines.append(str(len(target)))
                lines.append(' '.join(map(str,target)))
                ref = np.broadcast_to(dense,target)
            elif mode in (2,3,4):
                axis = int(rng.integers(rank))
                keepdims = bool(rng.integers(2))
                lines.append(f'{axis+1} {int(keepdims)}')
                if mode == 2:
                    ref = np.any(dense,axis=axis,keepdims=keepdims)
                elif mode == 3:
                    ref = np.all(dense,axis=axis,keepdims=keepdims)
                else:
                    ref = np.sum(dense,axis=axis,keepdims=keepdims,dtype=np.int64)
                if np.ndim(ref) == 0:
                    ref = np.reshape(ref,(1,))
            elif mode == 6:
                ref = (int(np.any(dense)),int(np.all(dense)))
        expected.append((mode,ref))
    return '\n'.join(lines)+'\n', expected


def check_results(output: str, expected: list[tuple[int,object]]):
    records = output.splitlines()
    if len(records) != len(expected):
        raise AssertionError(f'Expected {len(expected)} output records, received {len(records)}')
    for k, (record,(mode,ref)) in enumerate(zip(records,expected)):
        fields = list(map(int,record.split()))
        if mode == 6:
            if tuple(fields) != ref:
                raise AssertionError(f'Case {k}: full any/all {fields}, expected {ref}')
            continue
        rank = fields[0]
        shape = tuple(fields[1:rank+1])
        n = fields[rank+1]
        data = fields[rank+2:]
        ref = np.asarray(ref,dtype=np.int64)
        if shape != ref.shape or n != ref.size or data != ref.ravel(order='C').tolist():
            raise AssertionError(f'Case {k}, mode {mode}: wrong result, shape={shape}, '
                                 f'expected={ref.shape}, n={n}, expected={ref.size}, '
                                 f'got={data[:32]}, expected={ref.ravel()[:32]}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--trials',type=int,default=300)
    parser.add_argument('--seed',type=int,default=716)
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    payload,expected = make_trials(rng,args.trials)
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-bool-reduce-') as tmp:
            exe = compile_driver(Path(tmp),mode)
            completed = subprocess.run((str(exe),),input=payload,encoding='utf-8',
                                       text=True,capture_output=True,check=True)
            check_results(completed.stdout,expected)
            print(f'{mode}: PASS {args.trials} N-D Boolean broadcasting/reduction cases')


if __name__ == '__main__':
    main()
