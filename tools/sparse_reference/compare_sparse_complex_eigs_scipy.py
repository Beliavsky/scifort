#!/usr/bin/env python3
"""Randomized complex sparse eigs comparisons against NumPy and SciPy.

    python tools/sparse_reference/compare_sparse_complex_eigs_scipy.py --mode both --cases 80

NumPy's dense eigenspectrum is the independent accuracy standard.  SciPy's
sparse ARPACK eigs is also used on a subset of well-conditioned cases.
The Fortran routine reports actual residual norms, checked again by a
separate CSR matrix-vector multiply inside the reference driver.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sla

ROOT = Path(__file__).resolve().parents[2]
MODULES = (
    'src/kinds.f90', 'src/sparse.f90', 'src/sparse_coo.f90',
    'src/sparse_coo_nd.f90', 'src/sparse_coo_nd_complex.f90',
    'src/sparse_complex.f90', 'src/sparse_linalg_complex_direct.f90',
    'src/sparse_linalg_complex_eigs.f90',
)


def rank(values: np.ndarray, which: str) -> np.ndarray:
    metric = {
        'LM': abs(values),
        'SM': -abs(values),
        'LR': values.real,
        'SR': -values.real,
        'LI': values.imag,
        'SI': -values.imag,
    }[which]
    return np.argsort(-metric, kind='stable')


def compile_driver(build: Path, flags: list[str]) -> Path:
    compiler = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    objs = []
    for source in MODULES:
        obj = build / (Path(source).stem + '.o')
        subprocess.run([*compiler, '-c', str(ROOT / source), '-o', str(obj)],
                       check=True, cwd=build, capture_output=True, text=True)
        objs.append(obj)
    exe = build / 'complex_eigs_driver'
    subprocess.run([*compiler, str(ROOT / 'tools/sparse_reference/complex_eigs_driver.f90'),
                    *map(str, objs), '-o', str(exe)], check=True, cwd=build,
                   capture_output=True, text=True)
    return exe


def generate_cases(count: int, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    generated = []
    modes = ['LM', 'LR', 'SR', 'LI', 'SI', 'SM']
    for case_id in range(count):
        n = int(rng.integers(5, 15))
        which = modes[case_id % len(modes)]
        k = int(rng.integers(1, min(3, n-1)))
        # Random upper triangular, diagonal, full, and moderately sparse cases.
        # A small diagonal separation keeps selection well conditioned.
        kind = case_id % 4
        if kind == 0:
            a = np.triu(rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))) * 0.16
        elif kind == 1:
            a = np.diag(rng.normal(size=n) + 1j*rng.normal(size=n))
        elif kind == 2:
            a = (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))) / np.sqrt(n)
        else:
            a = (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))) / np.sqrt(n)
            a[rng.uniform(size=(n,n)) < 0.65] = 0
        a += np.diag(np.linspace(-1.2,1.7,n) + 0.32j*np.linspace(1.6,-1.1,n))
        if kind == 0:
            a = np.triu(a)
        rows, cols = np.nonzero(a)
        values = a[rows, cols]
        # Test duplicate-coordinate assembly independently of the dense input.
        if len(values) > 2 and case_id % 3 == 0:
            rows = np.repeat(rows, 2)
            cols = np.repeat(cols, 2)
            values = np.repeat(values / 2, 2)
        start = (rng.normal(size=n) + 1j*rng.normal(size=n)) if case_id % 2 else None
        generated.append(dict(a=a, n=n, nz=len(values), k=k, width=n,
                              cap=3, which=which, tol=2e-8, fmt=case_id%2,
                              start=start, rows=rows, cols=cols, values=values))
    return generated


def input_stream(cases: list[dict]) -> str:
    records = [str(len(cases))]
    for c in cases:
        records.append(f"{c['n']} {c['nz']} {c['k']} {c['width']} {c['cap']} "
                       f"{c['which']} {c['tol']:.17e} {c['fmt']} {int(c['start'] is not None)}")
        for row,col,z in zip(c['rows'], c['cols'], c['values']):
            records.append(f'{row+1} {col+1} {z.real:.17e} {z.imag:.17e}')
        if c['start'] is not None:
            records.extend(f'{z.real:.17e} {z.imag:.17e}' for z in c['start'])
    return '\n'.join(records)+'\n'


def validate(exe: Path, cases: list[dict]) -> tuple[int,int,float]:
    proc = subprocess.run([str(exe)], input=input_stream(cases), text=True,
                          capture_output=True, check=True, timeout=180)
    lines = proc.stdout.strip().splitlines()
    pos = 0
    maximum = 0.0
    scipy_count = 0
    for i,c in enumerate(cases):
        assert lines[pos] == f'CASE {i+1}', (i,lines[pos])
        pos += 1
        status, cycles, products, size = map(int,lines[pos].split())
        pos += 1
        assert status == 0, (i,c['which'],c['n'],status)
        assert size == c['k'] and cycles == 1 and products > 0
        output = []
        for j in range(size):
            real,imag,reported,verified = map(float,lines[pos].split())
            pos += 1
            output.append(complex(real,imag))
            assert np.isclose(reported,verified,atol=2e-10,rtol=2e-8), (i,reported,verified)
            assert verified < 1e-6, (i,verified)
        vals = np.linalg.eigvals(c['a'])
        desired = vals[rank(vals,c['which'])[:c['k']]]
        for x,y in zip(output,desired):
            err = abs(x-y)
            maximum = max(maximum,err)
            assert err <= 2e-6 * max(1.0,abs(y)), (i,c['which'],x,y,err)
        if i % 5 == 0 and c['which'] != 'SM' and c['n'] >= c['k']+3:
            # Explicit sparse-solver cross-check on a subset (dense NumPy is
            # always the primary reference for the full spectrum).
            try:
                comparison = sla.eigs(sparse.csr_array(c['a']),k=c['k'],
                                      which=c['which'],tol=1e-10,maxiter=2000,
                                      return_eigenvectors=False)
                assert np.max(np.abs(np.sort_complex(np.asarray(comparison))
                                     - np.sort_complex(np.asarray(output)))) < 2e-5
                scipy_count += 1
            except sla.ArpackNoConvergence:
                pass
    assert pos == len(lines), (pos,len(lines))
    return len(cases),scipy_count,maximum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=80)
    parser.add_argument('--seed',type=int,default=20261106)
    args = parser.parse_args()
    cases = generate_cases(args.cases,args.seed)
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        flags = (['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3'])
        with tempfile.TemporaryDirectory(prefix='scifort-complex-eigs-') as temp:
            exe = compile_driver(Path(temp),flags)
            completed,scipy_count,maximum = validate(exe,cases)
            print(f'{mode}: {completed} NumPy comparisons, '
                  f'{scipy_count} SciPy eigs cross-checks, '
                  f'largest eigenvalue error {maximum:.3e}',flush=True)


if __name__ == '__main__':
    main()
