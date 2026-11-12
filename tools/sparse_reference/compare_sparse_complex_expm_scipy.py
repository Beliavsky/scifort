#!/usr/bin/env python3
"""Compare complex CSR/CSC/operator expm actions against SciPy and dense expm.

The reference driver is compiled with gfortran. An intentionally small dense
matrix is passed to the driver to assemble sparse storage and an independent
matrix-free callback; the implementation itself never uses dense exponentials.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_array
from scipy.sparse.linalg import expm_multiply

ROOT = Path(__file__).resolve().parents[2]
FILES = [
    'src/kinds.f90', 'src/sparse.f90', 'src/sparse_coo.f90',
    'src/sparse_coo_nd.f90', 'src/sparse_coo_nd_complex.f90',
    'src/sparse_complex.f90', 'src/sparse_linalg_complex_expm.f90',
]


def compile_driver(directory: Path, mode: str) -> Path:
    flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace']
             if mode == 'checked' else ['-O3'])
    objects = []
    for file in FILES:
        obj = directory / (Path(file).stem + '.o')
        subprocess.run(['gfortran', '-std=f2018', *flags, '-J', str(directory),
                        '-I', str(directory), '-c', str(ROOT / file), '-o', str(obj)],
                       check=True, cwd=directory)
        objects.append(str(obj))
    executable = directory / 'complex_expm_driver'
    subprocess.run(['gfortran', '-std=f2018', *flags, '-J', str(directory),
                    '-I', str(directory),
                    str(ROOT / 'tools/sparse_reference/sparse_complex_expm_reference_driver.f90'),
                    *objects, '-o', str(executable)], check=True, cwd=directory)
    return executable


def random_case(rng: np.random.Generator, index: int):
    n = int(rng.integers(1,12))
    nrhs = int(rng.integers(1,4))
    scale = float(rng.choice([0.0,0.08,0.3,1.0,2.0,5.0]))
    a = scale*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    family = index % 7
    if family == 0:  # complex Hermitian
        a = (a + a.conj().T)/2
    if family == 1:  # triangular, often defective
        a = np.triu(a)
        np.fill_diagonal(a,a[0,0])
    if family == 2:  # sparse and singular
        a[rng.uniform(size=(n,n)) < .8] = 0
    if family == 3:  # large imaginary diagonal shift
        a = a+1j*float(rng.uniform(-4,4))*np.eye(n)
    if family == 4:  # diagonal
        a = np.diag(np.diag(a))
    if family == 5:  # nonsymmetric sparse
        a[rng.uniform(size=(n,n)) < .65] = 0
    b = rng.normal(size=(n,nrhs))+1j*rng.normal(size=(n,nrhs))
    if index % 10 == 0: b[:,0] = 0
    if index % 11 == 0: b[:,0] *= 1e-6
    if index % 13 == 0: b[:,-1] *= 1e5
    time = float(rng.uniform(-1.3,1.3))
    if index % 31 == 0: time = 0.0
    return np.asarray(a,np.complex128),np.asarray(b,np.complex128),time


def run(mode: str, cases: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    examples = [random_case(rng,i) for i in range(cases)]
    input_lines = [str(cases)]
    for idx,(a,b,t) in enumerate(examples):
        n,nrhs = b.shape
        start,stop,num,endpoint = (-0.4,0.8,4,idx%2 == 0)
        input_lines.append(f'{n} {nrhs} {t:.17e} {start:.17e} {stop:.17e} {num} {int(endpoint)}')
        for value in a.ravel(order='C'):
            input_lines.append(f'{value.real:.17e} {value.imag:.17e}')
        for value in b.ravel(order='C'):
            input_lines.append(f'{value.real:.17e} {value.imag:.17e}')
    with tempfile.TemporaryDirectory(prefix='scifort-complex-expm-') as d:
        executable = compile_driver(Path(d),mode)
        result = subprocess.run([str(executable)], input='\n'.join(input_lines)+'\n',
                                text=True, check=True, capture_output=True)
    lines = iter(result.stdout.splitlines())
    max_absolute_error = 0.0
    max_relative_error = 0.0
    for case_no,(a,b,t) in enumerate(examples):
        n,nrhs = b.shape
        dense_ref = expm(t*a) @ b
        scipy_ref = expm_multiply(csr_array(t*a),b)
        ref_disagreement = np.max(np.abs(dense_ref-scipy_ref))
        if ref_disagreement > 2e-10*max(1.0,np.max(np.abs(dense_ref))):
            raise AssertionError(f'References disagree in case {case_no}: {ref_disagreement}')
        for method in ('CSR','CSC','operator'):
            values = np.empty((n,nrhs),np.complex128)
            for i in range(n):
                for j in range(nrhs):
                    real,imag = next(lines).split()
                    values[i,j] = float(real)+1j*float(imag)
            error = float(np.max(np.abs(values-dense_ref)))
            rel = error/max(1.0,float(np.max(np.abs(dense_ref))))
            max_absolute_error = max(max_absolute_error,error)
            max_relative_error = max(max_relative_error,rel)
            if not np.isfinite(rel) or rel > 5e-11:
                raise AssertionError(f'{mode} case {case_no} {method}: '
                                     f'absolute {error}, relative {rel}; t={t}, a={a}, b={b}')
        start,stop,num,endpoint = (-0.4,0.8,4,case_no%2 == 0)
        grid = np.linspace(start,stop,num,endpoint=endpoint)
        interval_dense = [expm(tt*a) @ b for tt in grid]
        interval_scipy = expm_multiply(csr_array(a),b,start=start,stop=stop,
                                       num=num,endpoint=endpoint)
        # SciPy's interval stepping may lose some digits for large, strongly
        # nonnormal operators; the independent dense exponential is primary.
        if max(float(np.max(np.abs(interval_scipy[k]-interval_dense[k])))
               for k in range(num)) > 1e-6*max(1.0,max(np.max(np.abs(x)) for x in interval_dense)):
            raise AssertionError(f'Large SciPy interval discrepancy in case {case_no}')
        for method in ('CSR interval','CSC interval','operator interval'):
            values = np.empty((num,n,nrhs),np.complex128)
            for k in range(num):
                for i in range(n):
                    for j in range(nrhs):
                        real,imag = next(lines).split()
                        values[k,i,j] = float(real)+1j*float(imag)
            error = max(float(np.max(np.abs(values[k]-interval_dense[k])))
                        for k in range(num))
            rel = error/max(1.0,max(float(np.max(np.abs(x))) for x in interval_dense))
            max_absolute_error = max(max_absolute_error,error)
            max_relative_error = max(max_relative_error,rel)
            if not np.isfinite(rel) or rel > 5e-11:
                raise AssertionError(f'{mode} case {case_no} {method}: '
                                     f'absolute {error}, relative {rel}; a={a}, b={b}')
    try:
        extra = next(lines)
    except StopIteration:
        pass
    else:
        raise AssertionError(f'Unexpected driver output: {extra!r}')
    print(f'{mode}: {cases} randomized matrices x 6 backend/time variants passed; '
          f'max_abs={max_absolute_error:.3g}, max_rel={max_relative_error:.3g}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    parser.add_argument('--cases',type=int,default=180)
    parser.add_argument('--seed',type=int,default=20261106)
    args = parser.parse_args()
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        run(mode,args.cases,args.seed)
