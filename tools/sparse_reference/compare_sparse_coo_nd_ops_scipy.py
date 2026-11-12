#!/usr/bin/env python3
"""Compare SciFort N-D COO reductions, broadcasts, arithmetic, contractions with SciPy/NumPy.

Builds one Fortran streaming driver per compiler mode and exercises randomized
explicit-zero, duplicate, empty-shape, broadcast, and contracted-axis cases.
This does not require FPM or any third-party Fortran dependencies.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_array

ROOT = Path(__file__).resolve().parents[2]


@dataclass
class Example:
    shape: tuple[int, ...]
    coords: np.ndarray
    data: np.ndarray

    def dense(self) -> np.ndarray:
        # SciPy's coo_array is the independently implemented sparse reference.
        coords = tuple(self.coords[k] for k in range(len(self.shape)))
        return coo_array((self.data, coords), shape=self.shape).toarray()

    def input_lines(self) -> list[str]:
        lines = [f'{len(self.shape)} {len(self.data)}', ' '.join(map(str, self.shape))]
        for k in range(len(self.data)):
            lines.append(' '.join(map(str, self.coords[:, k]+1)) + f' {self.data[k]:.17g}')
        return lines


def example(rng: np.random.Generator, shape: tuple[int, ...], floating: bool = False) -> Example:
    nnz = int(rng.integers(0, 22)) if all(shape) else 0
    coords = np.zeros((len(shape), nnz), dtype=np.int64)
    for d, extent in enumerate(shape):
        if nnz:
            coords[d, :] = rng.integers(0, extent, size=nnz)
    vals = (rng.normal(0, 2, size=nnz) if floating else
            rng.integers(-4, 5, size=nnz).astype(np.float64))
    if nnz and rng.random() < .2:
        vals[0] = 0.0
    # Force repeated coordinates and explicit zeros in a significant fraction.
    if nnz > 2 and rng.random() < .45:
        coords[:, -1] = coords[:, 0]
    return Example(shape, coords, vals)


def random_shape(rng: np.random.Generator, rank: int | None = None) -> tuple[int, ...]:
    rank = int(rng.integers(1, 6)) if rank is None else rank
    return tuple(int(x) for x in rng.choice([0, 1, 1, 2, 2, 3, 4], size=rank))


def case(rng: np.random.Generator, mode: int, floating: bool = False) -> tuple[list[str], np.ndarray]:
    if mode in (1, 6):
        a = example(rng, random_shape(rng), floating)
        lines = [str(mode), *a.input_lines()]
        dense = a.dense()
        if mode == 6:
            return lines, np.asarray(np.sum(dense))
        k = int(rng.integers(0, len(a.shape)+1))
        axes = rng.choice(len(a.shape), size=k, replace=False).tolist()
        keep = bool(rng.integers(0, 2)) or k == len(a.shape)
        lines.append(f'{k} {int(keep)}')
        if k:
            lines.append(' '.join(str(i+1) for i in axes))
        return lines, np.sum(dense, axis=tuple(axes), keepdims=keep)
    if mode == 2:
        a = example(rng, random_shape(rng), floating)
        prefix = tuple(int(x) for x in rng.integers(1, 4, size=rng.integers(0, 3)))
        target = prefix + tuple(int(rng.integers(0, 4)) if d == 1 else d for d in a.shape)
        lines = [str(mode), *a.input_lines(), str(len(target)), ' '.join(map(str, target))]
        return lines, np.broadcast_to(a.dense(), target)
    if mode in (3, 4):
        common = random_shape(rng)
        # The input ranks may differ, and singleton axes may be broadcast.
        shape_a = tuple(1 if rng.random() < .35 else d for d in common)
        shape_b = tuple(1 if rng.random() < .35 else d for d in common)
        shape_a = shape_a[int(rng.integers(0, len(shape_a))):]
        shape_b = shape_b[int(rng.integers(0, len(shape_b))):]
        a, b = example(rng, shape_a, floating), example(rng, shape_b, floating)
        a_ref, b_ref = a.dense(), b.dense()
        lines = [str(mode), *a.input_lines(), *b.input_lines()]
        return lines, a_ref+b_ref if mode == 3 else a_ref*b_ref
    # Mode 5: general contraction, including outer and full contractions.
    shape_a = random_shape(rng)
    shape_b = list(random_shape(rng))
    k = int(rng.integers(0, min(len(shape_a), len(shape_b))+1))
    axes_a = rng.choice(len(shape_a), size=k, replace=False).tolist()
    axes_b = rng.choice(len(shape_b), size=k, replace=False).tolist()
    for ia, ib in zip(axes_a, axes_b):
        shape_b[ib] = shape_a[ia]
    b_shape = tuple(shape_b)
    a, b = example(rng, shape_a, floating), example(rng, b_shape, floating)
    lines = [str(mode), *a.input_lines(), *b.input_lines(), str(k)]
    if k:
        lines.extend([' '.join(str(i+1) for i in axes_a),
                      ' '.join(str(i+1) for i in axes_b)])
    result = np.tensordot(a.dense(), b.dense(), axes=(axes_a, axes_b))
    # The Fortran COO representation has positive rank; a scalar contraction
    # is represented by a length-one vector with its scalar value.
    if np.ndim(result) == 0:
        result = np.asarray([result.item()])
    return lines, result


def compile_driver(build: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    objects = []
    for name in ('kinds', 'sparse', 'sparse_coo', 'sparse_coo_nd', 'sparse_coo_nd_ops'):
        source = ROOT / 'src' / f'{name}.f90'
        obj = build / f'{name}.o'
        subprocess.run([*cc, '-c', str(source), '-o', str(obj)], check=True, cwd=build)
        objects.append(obj)
    exe = build / 'coo_nd_ops_reference'
    subprocess.run([*cc, str(ROOT / 'tools/sparse_reference/sparse_coo_nd_ops_reference_driver.f90'),
                    *map(str, objects), '-o', str(exe)], check=True, cwd=build)
    return exe


def run(mode: str, trials: int, seed: int, floating: bool) -> None:
    rng = np.random.default_rng(seed)
    inputs, expected = [], []
    for i in range(trials):
        m = i % 6 + 1
        lines, result = case(rng, m, floating)
        inputs.extend(lines)
        expected.append((m, np.asarray(result)))
    text_in = str(trials) + '\n' + '\n'.join(inputs) + '\n'
    with tempfile.TemporaryDirectory(prefix='scifort-coo-nd-ops-') as folder:
        exe = compile_driver(Path(folder),mode)
        process = subprocess.run([str(exe)],input=text_in,capture_output=True,text=True)
        if process.returncode != 0:
            print('Fortran stderr:',process.stderr[-1200:])
            print('Fortran stdout tail:',process.stdout[-1200:])
            raise RuntimeError('reference driver failure')
        lines_out = process.stdout.splitlines()
        if len(lines_out) != trials:
            raise AssertionError(f'{mode} output lines {len(lines_out)} != {trials}: {process.stdout[:300]}')
        for i, ((m, target), line) in enumerate(zip(expected, lines_out)):
            fields = line.split()
            status = int(fields[0])
            if status != 0:
                raise AssertionError(f'{mode} case {i}: operation {m} failed status {status}')
            rank = int(fields[1])
            if rank == 0:
                result = np.asarray(float(fields[2]))
            else:
                dims = tuple(int(x) for x in fields[2:2+rank])
                result = np.array([float(x) for x in fields[2+rank:]],dtype=np.float64).reshape(dims)
            if result.shape != target.shape or not np.allclose(result,target,rtol=2e-12,atol=2e-12):
                raise AssertionError(f'{mode} case {i} op {m}: {result.shape} vs {target.shape}, '
                                     f'{result} vs {target}')
    print(f'{mode}: PASS {trials} randomized SciPy/NumPy reference cases ('
          f'{"floating" if floating else "integer"} values)')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--trials',type=int,default=600)
    parser.add_argument('--seed',type=int,default=20261105)
    parser.add_argument('--floating',action='store_true',help='Use noninteger Gaussian values')
    args = parser.parse_args()
    modes = ('checked','optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        run(mode,args.trials,args.seed,args.floating)


if __name__ == '__main__':
    main()
