#!/usr/bin/env python3
"""Independent SciPy 1.17 reference checks for SciFort matrix structure/diagsvd.

Usage: python tools/compare_linalg_structure_diagsvd.py <compiled_fortran_driver>
Tests are generated deterministically. Fortran's input/output is fully numeric;
no SciPy implementation source is copied or used during Fortran compilation.
"""
from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass

import numpy as np
from scipy.linalg import bandwidth, diagsvd, ishermitian, issymmetric

RNG = np.random.default_rng(20261109)


@dataclass
class Case:
    op: int
    dtype: int
    ndim: int
    batch: int
    m: int
    n: int
    ns: int
    atol: float
    rtol: float
    data: np.ndarray
    expected: object
    invalid: bool


def integer_data(dtype: int, shape: tuple[int, ...]) -> np.ndarray:
    dtype_np = np.int32 if dtype == 1 else np.int64
    x = RNG.integers(-7, 8, size=shape, dtype=dtype_np)
    if dtype == 2 and RNG.random() < 0.3:
        x = np.where(x > 0, x + 9007199254740993, x - 9007199254740993).astype(dtype_np)
    return x


def make_data(dtype: int, shape: tuple[int, ...]) -> np.ndarray:
    if dtype in (1, 2):
        return integer_data(dtype, shape)
    a = RNG.normal(size=shape)
    if dtype == 4:
        a = a + 1j * RNG.normal(size=shape)
        a = np.array(a, dtype=np.complex128)
    else:
        a = np.array(a, dtype=np.float64)
    if RNG.random() < 0.2:
        a *= 10.0 ** RNG.uniform(-70, 70)
    return a


def reference(op: int, data: np.ndarray, ndim: int, m: int, n: int, atol: float, rtol: float):
    if op == 4:
        return diagsvd(data[0, :] if ndim == 2 else data, m, n)
    a = data[0] if ndim == 2 else data
    if op == 1:
        return bandwidth(a)
    func = issymmetric if op == 2 else ishermitian
    kw = {}
    if atol >= 0:
        kw['atol'] = atol
    if rtol >= 0:
        kw['rtol'] = rtol
    return func(a, **kw)


def cases():
    result = []
    for kind in range(1, 5):
        # Structured comparisons (including deliberately invalid shapes/batches).
        for k in range(550):
            ndim = 2 if k % 3 == 0 else 3
            batch = 1 if ndim == 2 else int(RNG.choice([0, 1, 2, 3, 4], p=[0.07, .26, .3, .27, .1]))
            m = int(RNG.integers(0, 9))
            n = m if RNG.random() < 0.73 else int(RNG.integers(0, 9))
            op = 1 + k % 3
            d = make_data(kind, (batch, m, n))
            if op in (2, 3) and m == n and m > 0:
                mode = k % 9
                if mode in (0, 1, 2, 3):
                    for b in range(batch):
                        aa = d[b]
                        if mode == 3:
                            aa[:] = 0
                        for i in range(m):
                            for j in range(i):
                                aa[i, j] = np.conjugate(aa[j, i]) if op == 3 else aa[j, i]
                            if op == 3 and kind == 4:
                                aa[i, i] = aa[i, i].real
                if mode == 4:
                    d[:] = 0
                if mode == 5 and kind in (3, 4):
                    if batch:
                        d[0, 0, 0] = np.nan
                if mode == 6 and m > 1 and kind in (3, 4):
                    if batch:
                        d[0, 0, 1] = np.inf
                        d[0, 1, 0] = np.inf
                if mode == 7 and kind == 4:
                    if batch:
                        d[0, 0, 0] = 1j * 1e-6
            if op == 1 and k % 7 == 0:
                d[:] = 0
                if m > 0 and n > 0 and batch:
                    d[0, m-1, 0] = 2
            at, rt = [(-1.,-1.), (0.,0.), (1e-13,0.), (0.,1e-8),
                      (1e-5,1e-5), (1e-1,0.), (0.,1e-1)][k % 7]
            try:
                expected = reference(op, d, ndim, m, n, at, rt)
                invalid = False
            except (ValueError, TypeError):
                expected, invalid = None, True
            result.append(Case(op,kind,ndim,batch,m,n,0,at,rt,d,expected,invalid))
        # Diagonal singular value construction, matching only the minimal dimension.
        for k in range(350):
            ndim = 2 if k % 3 == 0 else 3
            batch = 1 if ndim == 2 else int(RNG.choice([0,1,2,3],p=[.07,.31,.32,.30]))
            m, n = int(RNG.integers(0, 12)), int(RNG.integers(0, 12))
            ns = min(m,n)
            if k % 8 == 0:
                ns = ns + 1  # reject wrong number of singular values
            if k % 13 == 0:
                ns = max(0, ns - 1)
            d = make_data(kind, (batch, ns))
            try:
                expected = reference(4, d, ndim, m, n, -1, -1)
                invalid = False
            except (ValueError, TypeError):
                expected, invalid = None, True
            result.append(Case(4,kind,ndim,batch,m,n,ns,-1.,-1.,d,expected,invalid))
    # Extreme int64 exactness checks around the 2**53 boundary.
    for k in range(24):
        d = np.zeros((1,2,2), dtype=np.int64)
        d[0,0,1] = 9007199254740993
        d[0,1,0] = 9007199254740994
        op = 2 if k % 2 else 3
        atol = [-1.,0.1,1.0,2.0][k%4]
        expected = reference(op,d,2,2,2,atol,-1.)
        result.append(Case(op,2,2,1,2,2,0,atol,-1.,d,expected,False))
    return result


def to_line(v, typ):
    if typ == 4:
        return f'{float(v.real):.17e} {float(v.imag):.17e}'
    if typ in (1,2):
        return str(int(v))
    return f'{float(v):.17e}'


def verify(exe: str):
    cs = cases()
    inputs = [str(len(cs))]
    for c in cs:
        inputs.append(f'{c.op} {c.dtype} {c.ndim} {c.batch} {c.m} {c.n} {c.ns} {c.atol:.17e} {c.rtol:.17e}')
        if c.op == 4:
            for b in range(c.batch):
                for i in range(c.ns):
                    inputs.append(to_line(c.data[b,i], c.dtype))
        else:
            for b in range(c.batch):
                for j in range(c.n):
                    for i in range(c.m):
                        inputs.append(to_line(c.data[b,i,j], c.dtype))
    proc = subprocess.run([exe], input='\n'.join(inputs)+'\n', text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    output = iter(proc.stdout.splitlines())
    ok = 0
    rejected = 0
    total_entries = 0
    for index,c in enumerate(cs):
        try:
            status = int(next(output).strip())
        except (ValueError, StopIteration) as e:
            raise AssertionError(f'case {index}: missing status, stdout {proc.stdout[-200:]}, stderr {proc.stderr}') from e
        if c.invalid:
            if status != 1:
                raise AssertionError(f'case {index}: expected invalid status 1, got {status}; {c}')
            rejected += 1
            ok += 1
            continue
        if status != 0:
            raise AssertionError(f'case {index}: unexpected status {status}; {c}')
        if c.op == 1:
            result = [tuple(int(x) for x in next(output).split()) for _ in range(c.batch)]
            expected = c.expected
            if c.ndim == 2:
                good = result[0] == tuple(expected)
            else:
                good = (np.array_equal([x[0] for x in result], expected[0]) and
                        np.array_equal([x[1] for x in result], expected[1]))
            total_entries += 2*c.batch
        elif c.op in (2,3):
            result = [bool(int(next(output))) for _ in range(c.batch)]
            expected = c.expected
            good = result[0] == bool(expected) if c.ndim == 2 else np.array_equal(result, expected)
            total_entries += c.batch
        else:
            if c.dtype in (1,2):
                vals = [int(next(output)) for _ in range(c.batch*c.m*c.n)]
                actual = np.array(vals,dtype=np.int32 if c.dtype==1 else np.int64)
            elif c.dtype==3:
                vals = [float(next(output)) for _ in range(c.batch*c.m*c.n)]
                actual = np.array(vals,dtype=np.float64)
            else:
                vals = []
                for _ in range(c.batch*c.m*c.n):
                    r,im = next(output).split()
                    vals.append(complex(float(r),float(im)))
                actual = np.array(vals,dtype=np.complex128)
            actual = actual.reshape((c.batch,c.n,c.m)).swapaxes(1,2)
            expected = c.expected if c.ndim==3 else c.expected[None,:,:]
            good = np.array_equal(actual, expected, equal_nan=True) if c.dtype >=3 else np.array_equal(actual,expected)
            total_entries += c.batch*c.m*c.n
        if not good:
            raise AssertionError(f'case {index}: output mismatch, op={c.op}, kind={c.dtype}, '
                                 f'shape={c.data.shape}, tolerance={c.atol, c.rtol}, '
                                 f'got={result if c.op != 4 else actual}, expected={expected}')
        ok += 1
    assert next(output, None) is None, 'extra lines in reference driver output'
    print(f'{ok} / {len(cs)} independent SciPy comparison cases passed')
    print(f'{total_entries} output entries compared; {rejected} expected input rejections')
    print(f'NumPy {np.__version__}; SciPy {__import__("scipy").__version__}')
    return ok, total_entries, rejected


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('Usage: python tools/compare_linalg_structure_diagsvd.py <driver>')
    verify(sys.argv[1])
