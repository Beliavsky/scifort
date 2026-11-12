#!/usr/bin/env python3
"""Independent seeded SciPy comparisons for SciFort block_diag and khatri_rao.

Build the reference Fortran driver, then run:
  python tools/compare_linalg_block_diag_khatri_rao.py --exe /path/to/driver

The driver never imports or copies SciPy code. This script treats installed
SciPy as a numerical oracle for randomized input ranks, dtypes, and batches.
"""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass

import numpy as np
from scipy.linalg import block_diag, khatri_rao


@dataclass
class Case:
    name: str
    payload: list[str]
    expected: np.ndarray | None
    output_shape: tuple[int, ...]
    kind: str
    expected_error: bool = False


KINDS = {'J': np.dtype('int32'), 'I': np.dtype('int64'),
         'R': np.dtype('float64'), 'C': np.dtype('complex128')}


def code_of(array: np.ndarray) -> str:
    return {v: k for k, v in KINDS.items()}[np.asarray(array).dtype]


def raw3(array: np.ndarray) -> np.ndarray:
    array = np.asarray(array)
    if array.ndim == 0:
        return array.reshape(1, 1, 1)
    if array.ndim == 1:
        return array.reshape(1, 1, array.size)
    if array.ndim == 2:
        return array[np.newaxis, :, :]
    if array.ndim == 3:
        return array
    raise ValueError('unsupported case')


def entries(array: np.ndarray) -> list[str]:
    c = np.asarray(array, dtype=np.complex128)
    return [f'({np.real(v):.17e},{np.imag(v):.17e})'
            for v in c.reshape(-1, order='C')]


def values(rng: np.random.Generator, shape: tuple[int, ...], dtype: str) -> np.ndarray:
    if dtype in ('I', 'J'):
        return rng.integers(-13, 14, shape, dtype=KINDS[dtype])
    scale = float(rng.choice((1.0, 1.0e-40, 1.0e40, 1.0e-2)))
    a = rng.normal(size=shape) * scale
    if dtype == 'C':
        a = a + 1j*rng.normal(size=shape)*scale
    return np.asarray(a, dtype=KINDS[dtype])


def make_blocks(rng: np.random.Generator, count: int) -> list[Case]:
    result = []
    for case_id in range(count):
        inputs = []
        nblocks = int(rng.integers(0, 5))
        use_batches = bool(rng.random() < 0.55)
        main_batch = int(rng.choice((0, 1, 2, 3)))
        for _ in range(nblocks):
            dtype = str(rng.choice(tuple(KINDS)))
            rank = int(rng.choice((0, 1, 2, 3) if use_batches else (0, 1, 2)))
            if rank == 0:
                shape = ()
            elif rank == 1:
                shape = (int(rng.integers(0, 7)),)
            elif rank == 2:
                shape = (int(rng.integers(0, 5)), int(rng.integers(0, 6)))
            else:
                b = int(rng.choice((main_batch, main_batch, 1, 2)))
                shape = (b, int(rng.integers(0, 5)), int(rng.integers(0, 6)))
            inputs.append(values(rng, shape, dtype))
        bshape = 1
        nrows = 1 if not inputs else 0
        ncols = 0
        batched = any(a.ndim == 3 for a in inputs)
        for a in inputs:
            z = raw3(a)
            nrows += z.shape[1]
            ncols += z.shape[2]
            if a.ndim == 3 and z.shape[0] != 1 and bshape == 1:
                bshape = z.shape[0]
        out_shape = (bshape, nrows, ncols) if batched else (nrows, ncols)
        expected = None
        expected_error = False
        try:
            expected = block_diag(*inputs)
            out_shape = expected.shape
            dtype_out = code_of(expected)
        except ValueError:
            dtype_out = 'C'
            expected_error = True
        # Occasionally supply an intentionally incorrect caller-allocated shape.
        if rng.random() < 0.10 and not expected_error:
            shape_list = list(out_shape)
            shape_list[-1] += 1
            out_shape = tuple(shape_list)
            expected_error = True
        payload = ['B', f'{len(inputs)} {dtype_out} {len(out_shape)} '
                    + ' '.join(str(x) for x in ((1,)+out_shape if len(out_shape)==2 else out_shape))]
        for a in inputs:
            z = raw3(a)
            payload.append(f'{code_of(a)} {a.ndim} {z.shape[0]} {z.shape[1]} {z.shape[2]}')
            payload.extend(entries(z))
        result.append(Case(f'block#{case_id}', payload, expected, out_shape,
                           dtype_out, expected_error))
    return result


def make_khatri(rng: np.random.Generator, count: int) -> list[Case]:
    result = []
    combos = [(a,b) for a in ('J','I','R','C') for b in ('J','I','R','C')]
    for case_id in range(count):
        ta,tb = combos[int(rng.integers(0,len(combos)))]
        ra = int(rng.choice((2,3)))
        rb = int(rng.choice((2,3)))
        main_batch = int(rng.choice((0, 1, 2, 3)))
        ba = 1 if ra == 2 else int(rng.choice((main_batch,main_batch,1,2)))
        bb = 1 if rb == 2 else int(rng.choice((main_batch,main_batch,1,3)))
        n = int(rng.integers(0,8))
        m = int(rng.integers(0,8))
        k = int(rng.integers(0,6))
        shape_a = (n,k) if ra==2 else (ba,n,k)
        shape_b = (m,k) if rb==2 else (bb,m,k)
        a = values(rng,shape_a,ta)
        b = values(rng,shape_b,tb)
        expected = None
        expected_error = False
        bout = max(ba,bb)
        oshape = (n*m,k) if ra==2 and rb==2 else (bout,n*m,k)
        try:
            expected = khatri_rao(a,b)
            oshape = expected.shape
        except ValueError:
            expected_error = True
        if rng.random() < 0.10 and not expected_error:
            sh = list(oshape)
            sh[-2] += 1
            oshape = tuple(sh)
            expected_error = True
        bout,rows,cols = (1,)+oshape if len(oshape)==2 else oshape
        payload = ['K', f'{ta} {tb} {ra} {rb} {ba} {bb} {n} {m} {k} '
                   f'{bout} {rows} {cols}']
        payload.extend(entries(raw3(a)))
        payload.extend(entries(raw3(b)))
        result.append(Case(f'khatri#{case_id}',payload,expected,oshape,
                           ta+tb,expected_error))
    return result


def compare(executable: str, cases: list[Case]) -> tuple[float, int, int]:
    stream = [str(len(cases))]
    for case in cases:
        stream.extend(case.payload)
    run = subprocess.run([executable], input='\n'.join(stream)+'\n',
                         text=True, capture_output=True, timeout=120, check=False)
    if run.returncode != 0:
        raise RuntimeError(f'Fortran driver failed with {run.returncode}: '
                           f'{run.stderr[-4000:]}\nstdout={run.stdout[-500:]}')
    output = iter(run.stdout.splitlines())
    worst = 0.0
    failures = 0
    expected_errors = 0
    for case in cases:
        try:
            status = int(next(output).strip())
        except (StopIteration,ValueError) as ex:
            raise AssertionError(f'No valid Fortran status for {case.name}') from ex
        if case.expected_error:
            if status == 0:
                raise AssertionError(f'{case.name}: should reject invalid input or output shape')
            expected_errors += 1
            continue
        if status != 0:
            raise AssertionError(f'{case.name}: SciPy succeeded but SciFort returned status {status}')
        need = int(np.prod(case.output_shape, dtype=np.int64))
        got = np.empty(need, dtype=np.complex128)
        for i in range(need):
            try:
                real,imag = next(output).split()
            except (StopIteration,ValueError) as ex:
                raise AssertionError(f'Incomplete result: {case.name}, item {i}') from ex
            got[i] = complex(float(real), float(imag))
        assert case.expected is not None
        expected = np.asarray(case.expected,dtype=np.complex128)
        got = got.reshape(case.output_shape,order='C')
        if expected.shape != got.shape:
            raise AssertionError(f'{case.name}: expected {expected.shape}, got {got.shape}')
        if expected.size:
            scale = max(1.0, float(np.max(np.abs(expected))))
            difference = float(np.max(np.abs(expected-got))/scale)
            worst = max(worst,difference)
            if difference > 1.0e-12:
                raise AssertionError(f'{case.name}: scaled error={difference:.3e}; '
                                     f'expected={expected} got={got}')
        failures += 0
    if next(output, None) is not None:
        raise AssertionError('Unexpected surplus Fortran output')
    return worst, expected_errors, len(cases)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', required=True, help='path to built Fortran reference driver')
    parser.add_argument('--seed', type=int, default=921091)
    parser.add_argument('--blocks', type=int, default=650)
    parser.add_argument('--khatri', type=int, default=650)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    cases = make_blocks(rng,args.blocks) + make_khatri(rng,args.khatri)
    max_err, errors, count = compare(args.exe,cases)
    print(f'{count}/{count} SciPy reference comparisons passed; '
          f'{errors} expected invalid inputs; max scaled error {max_err:.3e}; '
          f'NumPy {np.__version__}')


if __name__ == '__main__':
    main()
