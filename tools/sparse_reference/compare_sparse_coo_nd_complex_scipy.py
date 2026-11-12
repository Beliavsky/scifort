#!/usr/bin/env python3
"""Randomized complex N-D COO comparisons against SciPy and NumPy.

Tests complex128 with duplicate coordinates, stored zeros, axis operations,
reshape order, singleton broadcasting, reductions, contractions and vdot.
Requires only NumPy, SciPy and GNU Fortran; FPM is not needed.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_array

ROOT = Path(__file__).resolve().parents[2]


class Tensor:
    def __init__(self, shape: tuple[int, ...], coords: np.ndarray, values: np.ndarray):
        self.shape, self.coords, self.values = shape, coords, values

    def dense(self) -> np.ndarray:
        return coo_array((self.values, tuple(self.coords)), shape=self.shape).toarray()

    def emit(self) -> list[str]:
        lines = [f"{len(self.shape)} {len(self.values)}", " ".join(map(str, self.shape))]
        for i, z in enumerate(self.values):
            lines.append(" ".join(map(str, self.coords[:, i]+1)) +
                         f" {z.real:.17g} {z.imag:.17g}")
        return lines


def random_shape(rng: np.random.Generator, rank: int | None = None) -> tuple[int, ...]:
    n = int(rng.integers(1, 6)) if rank is None else rank
    return tuple(int(x) for x in rng.choice([0, 1, 1, 2, 2, 3, 4], size=n))


def random_tensor(rng: np.random.Generator, shape: tuple[int, ...], *, floating: bool) -> Tensor:
    nnz = int(rng.integers(0, 23)) if all(shape) else 0
    coords = np.empty((len(shape), nnz), dtype=np.int64)
    for ax, extent in enumerate(shape):
        coords[ax] = rng.integers(0, extent, size=nnz) if nnz else []
    if floating:
        values = rng.normal(0, 1.5, nnz) + 1j*rng.normal(0, 1.5, nnz)
    else:
        values = rng.integers(-3, 4, size=nnz) + 1j*rng.integers(-3, 4, size=nnz)
    values = np.asarray(values, dtype=np.complex128)
    if nnz and rng.random() < .33:
        values[0] = 0j
    if nnz > 2 and rng.random() < .55:
        coords[:, -1] = coords[:, 0]  # force duplicate coordinates
    return Tensor(shape, coords, values)


def new_case(rng: np.random.Generator, mode: int, floating: bool) -> tuple[list[str], np.ndarray]:
    if mode in (6, 7):
        shared = random_shape(rng)
        sa = tuple(1 if rng.random() < .35 else x for x in shared)
        sb = tuple(1 if rng.random() < .35 else x for x in shared)
        sa = sa[int(rng.integers(0, len(sa))):]
        sb = sb[int(rng.integers(0, len(sb))):]
        a, b = random_tensor(rng, sa, floating=floating), random_tensor(rng, sb, floating=floating)
        target = a.dense()+b.dense() if mode == 6 else a.dense()*b.dense()
        return [str(mode), *a.emit(), *b.emit()], np.asarray(target)
    if mode == 9:
        shape_a, shape_b = random_shape(rng), list(random_shape(rng))
        k = int(rng.integers(0, min(len(shape_a),len(shape_b))+1))
        axes_a = rng.choice(len(shape_a), size=k, replace=False).tolist()
        axes_b = rng.choice(len(shape_b), size=k, replace=False).tolist()
        for ia, ib in zip(axes_a, axes_b):
            shape_b[ib] = shape_a[ia]
        a = random_tensor(rng, shape_a, floating=floating)
        b = random_tensor(rng, tuple(shape_b), floating=floating)
        lines = [str(mode), *a.emit(), *b.emit(), str(k)]
        if k:
            lines.extend([" ".join(str(x+1) for x in axes_a),
                          " ".join(str(x+1) for x in axes_b)])
        target = np.tensordot(a.dense(), b.dense(), axes=(axes_a, axes_b))
        if np.ndim(target) == 0:
            target = np.asarray([complex(target)])
        return lines, np.asarray(target)
    if mode == 10:
        shape = random_shape(rng)
        a = random_tensor(rng, shape, floating=floating)
        b = random_tensor(rng, shape, floating=floating)
        return [str(mode), *a.emit(), *b.emit()], np.asarray(np.vdot(a.dense(), b.dense()))
    shape = random_shape(rng, rank=2 if mode == 11 else None)
    a = random_tensor(rng, shape, floating=floating)
    lines = [str(mode), *a.emit()]
    dense = a.dense()
    if mode in (1, 13):
        if mode == 13:
            lines.append(rng.choice(['C','F']))
        return lines, np.asarray(dense)
    if mode == 2:
        return lines, np.conjugate(dense)
    if mode == 3:
        perm = rng.permutation(len(shape)).tolist()
        lines.extend([str(len(shape)), " ".join(str(x+1) for x in perm)])
        return lines, np.transpose(dense, perm)
    if mode == 4:
        order = rng.choice(['C','F'])
        total = int(np.prod(shape))
        # Include empty arrays and change rank; avoid zero-dimensional outputs.
        target_shape = (total,) if rng.random() < .5 else (1, total, 1)
        lines.extend([f'{len(target_shape)} {order}'," ".join(map(str,target_shape))])
        return lines, np.reshape(dense,target_shape,order=order)
    if mode == 5:
        prefix = tuple(int(x) for x in rng.integers(1,4,size=int(rng.integers(0,3))))
        target_shape = prefix + tuple(int(rng.integers(0,4)) if x == 1 else x for x in shape)
        lines.extend([str(len(target_shape))," ".join(map(str,target_shape))])
        return lines, np.broadcast_to(dense,target_shape)
    if mode == 8:
        k = int(rng.integers(0,len(shape)+1))
        axes = rng.choice(len(shape),size=k,replace=False).tolist()
        keepdims = bool(rng.integers(0,2)) or k == len(shape)
        lines.append(f'{k} {int(keepdims)}')
        if k:
            lines.append(" ".join(str(x+1) for x in axes))
        return lines,np.sum(dense,axis=tuple(axes),keepdims=keepdims)
    if mode == 11:
        return lines, dense.conjugate().T
    if mode == 12:
        return lines, np.asarray(np.sum(dense))
    raise AssertionError(f"unexpected mode {mode}")


def compile_driver(build: Path, mode: str) -> Path:
    flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
    cc = ['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
    objects = []
    for name in ('kinds','sparse','sparse_coo','sparse_coo_nd',
                 'sparse_coo_nd_complex','sparse_coo_nd_complex_ops'):
        obj = build / f'{name}.o'
        subprocess.run([*cc,'-c',str(ROOT/'src'/f'{name}.f90'),'-o',str(obj)],
                       cwd=build,check=True)
        objects.append(str(obj))
    exe = build/'complex_reference_driver'
    subprocess.run([*cc,str(ROOT/'tools/sparse_reference/sparse_coo_nd_complex_reference_driver.f90'),
                    *objects,'-o',str(exe)],check=True,cwd=build)
    return exe


def run(mode: str, trials: int, seed: int, floating: bool) -> None:
    rng = np.random.default_rng(seed)
    inputs, targets = [], []
    for i in range(trials):
        op = (i % 13) + 1
        lines, result = new_case(rng,op,floating)
        inputs.extend(lines)
        targets.append((op,np.asarray(result,dtype=np.complex128)))
    payload = f'{trials}\n' + '\n'.join(inputs) + '\n'
    with tempfile.TemporaryDirectory(prefix='scifort-complex-coo-') as folder:
        exe = compile_driver(Path(folder),mode)
        output = subprocess.run([str(exe)],input=payload,text=True,capture_output=True)
        if output.returncode != 0:
            print('driver stdout tail:',output.stdout[-1200:])
            print('driver stderr tail:',output.stderr[-1200:])
            raise AssertionError(f'{mode}: Fortran driver failed')
        rows = output.stdout.splitlines()
        if len(rows) != trials:
            raise AssertionError(f'{mode}: expected {trials} output lines, got {len(rows)}')
        for i, ((op,target),row) in enumerate(zip(targets,rows)):
            fields = row.split()
            status, rank = map(int,fields[:2])
            if status != 0:
                raise AssertionError(f'{mode} trial={i} op={op} status={status}')
            if rank:
                dims = tuple(int(x) for x in fields[2:rank+2])
                values = np.asarray([float(x) for x in fields[rank+2:]],dtype=np.float64)
                result = (values[::2] + 1j*values[1::2]).reshape(dims)
            else:
                result = np.complex128(float(fields[2])+1j*float(fields[3]))
            if result.shape != target.shape or not np.allclose(result,target,rtol=1e-11,atol=1e-11):
                raise AssertionError(f'{mode} trial={i} op={op}: result shape {result.shape}, expected '
                                     f'{target.shape}, result={result!r}, expected={target!r}')
    print(f'{mode}: PASS {trials} randomized SciPy/NumPy complex COO cases '
          f'({"floating" if floating else "integer"} values)',flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--trials',type=int,default=390)
    parser.add_argument('--seed',type=int,default=20261105)
    parser.add_argument('--floating',action='store_true')
    args = parser.parse_args()
    for mode in (('checked','optimized') if args.mode == 'both' else (args.mode,)):
        run(mode,args.trials,args.seed,args.floating)


if __name__ == '__main__':
    main()
