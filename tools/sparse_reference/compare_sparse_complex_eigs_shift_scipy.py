#!/usr/bin/env python3
"""Independent randomized shift-invert complex eigs comparisons.

Run as:
    python tools/sparse_reference/compare_sparse_complex_eigs_shift_scipy.py --mode both --cases 80

References are NumPy dense eigenvalues and optionally SciPy's shift-invert
ARPACK eigs.  The Fortran solver factors the shifted CSR matrix only once.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

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


def order_by_shift(eigenvalues: np.ndarray, shift: complex, which: str) -> np.ndarray:
    shifted = 1.0 / (eigenvalues - shift)
    priority = {
        'LM': abs(shifted),
        'SM': -abs(shifted),
        'LR': shifted.real,
        'SR': -shifted.real,
        'LI': shifted.imag,
        'SI': -shifted.imag,
    }[which]
    return np.argsort(-priority, kind='stable')


def compile_driver(build: Path, flags: list[str]) -> Path:
    cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    objs = []
    for source in MODULES:
        obj = build / (Path(source).stem + '.o')
        subprocess.run([*cc, '-c', str(ROOT / source), '-o', str(obj)],
                       check=True, cwd=build, capture_output=True, text=True)
        objs.append(obj)
    exe = build / 'complex_eigs_shift_driver'
    subprocess.run([*cc, str(ROOT / 'tools/sparse_reference/complex_eigs_shift_driver.f90'),
                    *map(str, objs), '-o', str(exe)], check=True,
                   cwd=build, capture_output=True, text=True)
    return exe


def generate_cases(count: int, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    cases = []
    modes = ('LM', 'LR', 'SR', 'LI', 'SI', 'SM')
    for case_id in range(count):
        n = int(rng.integers(6, 13))
        k = int(rng.integers(1, min(3, n-2) + 1))
        kind = case_id % 4
        if kind == 0:
            a = np.triu(rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))) * 0.1
        elif kind == 1:
            a = np.diag(rng.normal(size=n) + 1j*rng.normal(size=n))
        elif kind == 2:
            a = (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))) / np.sqrt(n)
        else:
            a = (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))) / np.sqrt(n)
            a[rng.uniform(size=(n,n)) < 0.6] = 0
        a += np.diag(np.linspace(-2.2,3.1,n) + 0.25j*np.linspace(1.6,-1.1,n))
        if kind == 0:
            a = np.triu(a)
        eigenvalues = np.linalg.eigvals(a)
        anchor = eigenvalues[int(rng.integers(0,n))]
        shift = anchor + (0.2+0.17j) * rng.uniform(0.8,1.2)
        rows, cols = np.nonzero(a)
        values = a[rows, cols]
        if len(values) > 2 and case_id % 3 == 0:
            rows, cols, values = np.repeat(rows,2), np.repeat(cols,2), np.repeat(values/2,2)
        start = (rng.normal(size=n) + 1j*rng.normal(size=n)) if case_id % 2 else None
        cases.append(dict(n=n, k=k, a=a, which=modes[case_id % len(modes)],
                          shift=shift, nz=len(values), rows=rows, cols=cols,
                          values=values, width=n, cap=3, fmt=case_id%2,
                          tol=1e-8, start=start))
    return cases


def input_text(cases: list[dict]) -> str:
    text = [str(len(cases))]
    for c in cases:
        s = c['shift']
        text.append(f"{c['n']} {c['nz']} {c['k']} {c['width']} {c['cap']} "
                    f"{c['which']} {c['tol']:.17e} {c['fmt']} "
                    f"{int(c['start'] is not None)} {s.real:.17e} {s.imag:.17e}")
        for r, col, val in zip(c['rows'],c['cols'],c['values']):
            text.append(f'{r+1} {col+1} {val.real:.17e} {val.imag:.17e}')
        if c['start'] is not None:
            text.extend(f'{val.real:.17e} {val.imag:.17e}' for val in c['start'])
    return '\n'.join(text) + '\n'


def validate(exe: Path, cases: list[dict]) -> tuple[int,int,float]:
    proc = subprocess.run([str(exe)],input=input_text(cases),text=True,
                          capture_output=True,check=True,timeout=120)
    lines = proc.stdout.strip().splitlines()
    pos, scipy_checks, maximum = 0,0,0.0
    for i,c in enumerate(cases):
        assert lines[pos] == f'CASE {i+1}', (i,lines[pos]); pos += 1
        status,iterations,products,k = map(int,lines[pos].split()); pos+=1
        assert status == 0 and iterations == 1 and products >= c['width'] and k == c['k'], (i,status)
        result = []
        for j in range(k):
            re,im,reported,verified = map(float,lines[pos].split()); pos+=1
            result.append(complex(re,im))
            assert np.isclose(reported,verified,atol=5e-10,rtol=1e-7), (i,reported,verified)
            assert verified < 1e-6, (i,verified)
        spectrum = np.linalg.eigvals(c['a'])
        desired = spectrum[order_by_shift(spectrum,c['shift'],c['which'])[:k]]
        for a,b in zip(result,desired):
            err=abs(a-b)
            maximum=max(maximum,err)
            assert err < 3e-6*max(1.0,abs(b)), (i,c['which'],a,b,err)
        if i % 6 == 0 and c['which']=='LM' and c['n'] >= k+3:
            try:
                comparison = sla.eigs(sparse.csr_array(c['a']),k=k,sigma=c['shift'],
                                      which='LM',tol=1e-10,maxiter=2500,
                                      return_eigenvectors=False)
                scipy_result = np.sort_complex(np.asarray(comparison))
                expected = np.sort_complex(np.asarray(result))
                assert np.max(abs(scipy_result-expected)) < 5e-6, (i,scipy_result,expected)
                scipy_checks+=1
            except sla.ArpackNoConvergence:
                pass
    assert pos == len(lines), (pos,len(lines))
    return len(cases),scipy_checks,maximum


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=80)
    parser.add_argument('--seed',type=int,default=20261106)
    args=parser.parse_args()
    cases=generate_cases(args.cases,args.seed)
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
        with tempfile.TemporaryDirectory(prefix='scifort-complex-shift-') as temp:
            result=validate(compile_driver(Path(temp),flags),cases)
        print(f'{mode}: {result[0]} NumPy comparisons, '
              f'{result[1]} SciPy shift-invert eigs cross-checks, '
              f'largest eigenvalue error {result[2]:.3e}',flush=True)


if __name__=='__main__':
    main()
