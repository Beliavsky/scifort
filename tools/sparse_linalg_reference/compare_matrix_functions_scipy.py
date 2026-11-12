#!/usr/bin/env python3
"""Compare SciFort CSR/CSC matrix powers and inverses with SciPy and NumPy.

Runs a single Fortran process per build mode, then checks every case against
SciPy sparse.linalg and independent NumPy dense linear algebra calculations.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
import scipy
from scipy.sparse import coo_array, csc_array, csr_array
from scipy.sparse.linalg import inv as scipy_inv
from scipy.sparse.linalg import matrix_power as scipy_matrix_power

ROOT = Path(__file__).resolve().parents[2]


def get_case(rng: np.random.Generator, case_no: int) -> tuple[np.ndarray, list[tuple[int,int,float]], int]:
    n = int(rng.integers(1, 15))
    a = rng.standard_normal((n, n)) * (rng.random((n, n)) < .22) * .15
    shape = case_no % 10
    if shape == 0:
        a = np.diag(rng.uniform(.5, 1.5, n))
    elif shape == 1:
        a = np.triu(a)
        np.fill_diagonal(a, rng.choice([-1., 1.], n) * rng.uniform(.5, 1.5, n))
    elif shape == 2:
        a = np.tril(a)
        np.fill_diagonal(a, rng.uniform(.5, 1.5, n))
    elif shape == 3:
        # Permutation matrix, including zero diagonal entries requiring pivoting.
        a = np.eye(n)[rng.permutation(n)]
    elif shape == 4:
        a = a + np.eye(n) * rng.uniform(.8, 1.2)
    elif shape == 5:
        a = a - a.T + np.eye(n) * .9
    elif shape == 6:
        a = np.zeros((n, n))
        np.fill_diagonal(a, rng.choice([-1.,1.], n))
        if n > 1:
            a[0,0] = 0.
            a[0,1] = 1.
            a[1,0] = 1.
            a[1,1] = 0.
    elif shape == 7:
        a = a @ a.T + np.eye(n) * .5
    elif shape == 8:
        a = np.diag(rng.uniform(1.3, 2.0, n))
        if n > 1:
            a += np.diag(rng.uniform(-.25,.25,n-1),k=1)
    else:
        a = a + np.eye(n) * rng.choice([-.7,.7])
    powers = [0,1,2,3,4,5,7,8,11,16]
    power = int(powers[case_no % len(powers)])
    rr, cc = np.nonzero(a)
    vals = a[rr,cc].copy()
    triples = [(int(r)+1, int(c)+1, float(v)) for r,c,v in zip(rr,cc,vals)]
    if triples and case_no % 7 == 0:
        k = int(rng.integers(len(triples)))
        r,c,v = triples[k]
        triples[k] = (r,c,v*.7)
        triples.append((r,c,v*.3))
    if case_no % 4 == 0:
        # Explicitly stored zero in a legal CSR/CSC coordinate.
        triples.append((int(rng.integers(n))+1, int(rng.integers(n))+1, 0.0))
    if case_no % 9 == 0:
        rng.shuffle(triples)
    sparse = coo_array((np.array([v for _,_,v in triples], dtype=float),
                        (np.array([r-1 for r,_,_ in triples], dtype=int),
                         np.array([c-1 for _,c,_ in triples], dtype=int))), shape=(n,n)).tocsr()
    return sparse.toarray(), triples, power


def compile_driver(mode: str, build: Path) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    compiler = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    objects = []
    for name, path in (
        ('kinds','src/kinds.f90'), ('sparse','src/sparse.f90'),
        ('lsmr','src/optimize/lsmr.f90'), ('sparse_linalg','src/sparse_linalg.f90')
    ):
        obj = build/f'{name}.o'
        subprocess.run([*compiler, '-c', str(ROOT/path), '-o', str(obj)],
                       cwd=build, check=True, stdout=subprocess.PIPE)
        objects.append(str(obj))
    exe = build/'matrix_functions_driver'
    subprocess.run([*compiler, str(ROOT/'tools/sparse_linalg_reference'/
                      'sparse_matrix_functions_reference_driver.f90'),
                    *objects, '-o', str(exe)], cwd=build, check=True,
                   stdout=subprocess.PIPE)
    return exe


def check(mode: str, cases: list[tuple[np.ndarray,list[tuple[int,int,float]],int]]) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-matrix-functions-') as temp:
        build = Path(temp)
        exe = compile_driver(mode, build)
        payload = [str(len(cases))]
        for a, entries, power in cases:
            payload.append(f'{len(a)} {len(entries)} {power}')
            payload.extend(f'{i} {j} {v:.17g}' for i,j,v in entries)
        proc = subprocess.run([str(exe)], input='\n'.join(payload)+'\n', text=True,
                              capture_output=True, check=True, cwd=build)
        lines = proc.stdout.splitlines()
        if len(lines) != len(cases)*6:
            raise AssertionError(f'expected {6*len(cases)} lines, got {len(lines)}: {proc.stdout[-500:]}')
        pivot_cases = 0
        checked = 0
        for idx, (a, entries, power) in enumerate(cases):
            n = len(a)
            head = [int(x) for x in lines[6*idx].split()]
            head_inv = [int(x) for x in lines[6*idx+3].split()]
            if head != [idx+1,0,0,power] or head_inv != [idx+1,0,0]:
                raise AssertionError(f'case={idx} power={head}, inv={head_inv}')
            result = [np.fromstring(lines[6*idx+j], sep=' ').reshape((n,n),order='F')
                      for j in [1,2,4,5]]
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                csr = csr_array(a)
                csc = csc_array(a)
                reference = np.linalg.matrix_power(a,power)
                sparse_csr = scipy_matrix_power(csr,power).toarray()
                sparse_csc = scipy_matrix_power(csc,power).toarray()
                dense_inverse = np.linalg.inv(a)
                sparse_inverse = scipy_inv(csc)
                if hasattr(sparse_inverse, 'toarray'):
                    sparse_inverse = sparse_inverse.toarray()
                else:
                    sparse_inverse = np.asarray(sparse_inverse)
            scale = max(1.,np.max(np.abs(reference)))
            for j, target in enumerate((reference,reference)):
                error = np.max(np.abs(result[j]-target))
                assert error <= 8e-11*scale, (mode,idx,'power',j,error,scale)
            assert np.allclose(sparse_csr,reference,atol=8e-11*scale,rtol=0)
            assert np.allclose(sparse_csc,reference,atol=8e-11*scale,rtol=0)
            inv_scale = max(1.,np.max(np.abs(dense_inverse)))
            for j in (2,3):
                error = np.max(np.abs(result[j]-dense_inverse))
                assert error < 2e-10*inv_scale, (mode,idx,'inverse',j,error)
                residual = np.max(np.abs(a@result[j]-np.eye(n)))
                assert residual < 2e-10, (mode,idx,'inverse residual',j,residual)
            assert np.allclose(sparse_inverse,dense_inverse,atol=2e-10*inv_scale,rtol=0)
            if np.any((np.diag(a) == 0.) & (np.max(np.abs(a),axis=0) > 0.)):
                pivot_cases += 1
            checked += 1
        print(f'{mode}: {checked} matrices, {4*checked} SciFort outputs, '
              f'{pivot_cases} cases with zero diagonal entries, '
              'SciPy CSR/CSC and NumPy dense references passed', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['checked','optimized','both'],default='both')
    parser.add_argument('--cases',type=int,default=200)
    parser.add_argument('--seed',type=int,default=20261105)
    args = parser.parse_args()
    if args.cases <= 0:
        parser.error('--cases must be positive')
    rng = np.random.default_rng(args.seed)
    cases = [get_case(rng,i) for i in range(args.cases)]
    for mode in (['checked','optimized'] if args.mode == 'both' else [args.mode]):
        check(mode,cases)
    print('SciPy version:',scipy.__version__)

if __name__ == '__main__':
    main()
