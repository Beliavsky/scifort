#!/usr/bin/env python3
"""Compare complex Krylov matrix-function actions with dense SciPy references.

The standalone Fortran reference driver exercises vector and block, CSR and
CSC, plus matrix-free paths. On a subset, SciPy's public Krylov method is
cross-checked too. NumPy/SciPy dense matrix functions are the accuracy oracle.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg, sparse
from scipy.sparse import linalg as sparse_linalg

ROOT = Path(__file__).resolve().parents[2]
FILES = (
    'src/kinds.f90', 'src/sparse.f90', 'src/sparse_coo.f90',
    'src/sparse_coo_nd.f90', 'src/sparse_coo_nd_complex.f90',
    'src/sparse_complex.f90', 'src/sparse_linalg_complex_funm.f90',
)


def compile_driver(directory: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    objects = []
    for path in FILES:
        obj = directory / (Path(path).stem + '.o')
        subprocess.run(['gfortran', '-std=f2018', *flags, '-J', str(directory),
                        '-I', str(directory), '-c', str(ROOT / path), '-o', str(obj)],
                       check=True, cwd=directory)
        objects.append(str(obj))
    executable = directory / 'complex_funm_driver'
    subprocess.run(['gfortran', '-std=f2018', *flags, '-J', str(directory),
                    '-I', str(directory),
                    str(ROOT / 'tools/sparse_reference/sparse_complex_funm_reference_driver.f90'),
                    *objects, '-o', str(executable)], check=True, cwd=directory)
    return executable


def sample(rng: np.random.Generator, index: int):
    n = int(rng.integers(1, 14))
    nrhs = int(rng.integers(1, 4))
    function = index % 3
    family = (index // 3) % 7
    scale = float(rng.choice([0.0, 0.08, 0.2, 0.4, 0.9]))
    a = scale * (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)))
    hermitian = False
    if family == 0:
        a = (a + a.conj().T)/2
        hermitian = True
    elif family == 1:
        a = np.diag(np.diag(a))
    elif family == 2:
        a = np.triu(a)
        np.fill_diagonal(a, a[0,0])
    elif family == 3:
        a[rng.random(size=(n,n)) < 0.75] = 0
    elif family == 4:
        a = (a - a.conj().T)/2
    elif family == 5:
        a = np.tril(a, -1)
    b = rng.normal(size=(n,nrhs)) + 1j*rng.normal(size=(n,nrhs))
    if index % 17 == 0:
        b[:,0] = 0
    if index % 29 == 0:
        b *= 1e-4
    time = complex(float(rng.uniform(-0.35,0.35)),float(rng.uniform(-0.25,0.25)))
    if index % 20 == 0:
        time = 0j
    if index % 3 == 0 and n >= 6:
        m = min(n, max(3, n // 2))
        cap = 12
    else:
        m = n
        cap = 2
    # A larger Krylov space protects strongly nonnormal cases.
    return a,b,function,hermitian,m,cap,time


def function_eval(which: int, matrix: np.ndarray) -> np.ndarray:
    if which == 0:
        return matrix@matrix + (0.2-0.3j)*matrix + (0.75+0.1j)*np.eye(matrix.shape[0])
    if which == 1:
        return linalg.expm(matrix)
    if which == 2:
        return linalg.sinm(matrix)
    raise ValueError(which)


def run(mode: str, cases: int, seed: int, scipy_cases: int) -> None:
    rng = np.random.default_rng(seed)
    examples = [sample(rng, i) for i in range(cases)]
    lines = [str(cases)]
    for a,b,function,hermitian,m,cap,time in examples:
        n,nrhs = b.shape
        lines.append(f'{n} {nrhs} {function} {int(hermitian)} {m} {cap} '
                     f'{time.real:.17e} {time.imag:.17e}')
        for z in a.flat:
            lines.append(f'{z.real:.17e} {z.imag:.17e}')
        for z in b.flat:
            lines.append(f'{z.real:.17e} {z.imag:.17e}')
    with tempfile.TemporaryDirectory(prefix='scifort-complex-funm-') as name:
        executable = compile_driver(Path(name), mode)
        result = subprocess.run([str(executable)], text=True, input='\n'.join(lines)+'\n',
                                capture_output=True, check=True)
    records = iter(result.stdout.splitlines())
    max_error = 0.0
    comparison_count = 0
    scipy_count = 0
    for case_no,(a,b,function,hermitian,m,cap,time) in enumerate(examples):
        n,nrhs = b.shape
        expected = function_eval(function,time*a) @ b
        outputs = []
        for method in range(6):
            k = 1 if method < 3 else nrhs
            output = np.empty((n,k),dtype=np.complex128)
            for row in range(n):
                for col in range(k):
                    re,im = next(records).split()
                    output[row,col] = float(re)+1j*float(im)
            outputs.append(output)
        statuses = [int(x) for x in next(records).split()]
        if statuses != [0]*6:
            raise AssertionError(f'{mode} case {case_no}: solver statuses {statuses}, '
                                 f'family {(case_no//3)%7}, n={n}, m={m}, cap={cap}, '
                                 f'function {function}')
        for method,output in enumerate(outputs):
            ref = expected[:,:1] if method < 3 else expected
            error = float(np.max(abs(output-ref)))
            tolerance = 3e-8*max(1.0,float(np.max(abs(ref))))
            if not np.isfinite(error) or error > tolerance:
                raise AssertionError(f'{mode} case {case_no}, method {method}: '
                                     f'error {error:.3g}, tolerance {tolerance:.3g}, '
                                     f'n={n}, m={m}, cap={cap}, f={function}, A={a}')
            max_error=max(max_error,error)
            comparison_count+=1
        if (case_no < scipy_cases and m == n and np.linalg.norm(b[:,0]) > 0
                and np.linalg.norm(a) > 0 and time != 0j):
            # Exercise the public SciPy reference interface with its own
            # callable on the same small input. Some combinations of SciPy
            # 1.17's stopping rule can deviate more from dense evaluation;
            # compare its norm error independently rather than using it as a
            # source of truth for our solver.
            try:
                scipy_result = sparse_linalg.funm_multiply_krylov(
                    lambda mat: function_eval(function,np.asarray(mat)),
                    sparse.csr_array(a), b[:,0],
                    t=time, assume_a='hermitian' if hermitian else 'general',
                    restart_every_m=m,max_restarts=cap,rtol=1e-9)
                if scipy_result.shape != (n,) or not np.all(np.isfinite(scipy_result)):
                    raise AssertionError(f'bad SciPy result case {case_no}')
                if np.max(abs(scipy_result-expected[:,0])) <= 3e-8*max(1.0,float(np.max(abs(expected[:,0])))):
                    scipy_count+=1
            except (ValueError,RuntimeError,NotImplementedError,IndexError):
                # SciPy may decline a supplied complex time or restricted
                # convergence budget; dense reference remains primary.
                pass
    try:
        unexpected = next(records)
        raise AssertionError(f'unexpected trailing driver output: {unexpected}')
    except StopIteration:
        pass
    print(f'{mode}: {cases} random complex matrices, {comparison_count} method comparisons PASS; '
          f'{scipy_count} public SciPy callable cross-checks; max absolute error {max_error:.3e}')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    parser.add_argument('--cases', type=int, default=100)
    parser.add_argument('--seed', type=int, default=20261106)
    parser.add_argument('--scipy-cases', type=int, default=35)
    args = parser.parse_args()
    for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
        run(mode,args.cases,args.seed,args.scipy_cases)


if __name__ == '__main__':
    main()
