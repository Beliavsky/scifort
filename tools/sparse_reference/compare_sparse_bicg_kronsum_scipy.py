#!/usr/bin/env python3
"""Independent SciPy 1.x comparisons for SciFort BiCG, CGS and CSR kronsum.

Compiles a native Fortran reference driver; checks value and sparse structure
for randomized Kronecker sums, and compares iterative solvers on nonsymmetric,
strictly diagonally dominant systems using both true residuals and SciPy.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

ROOT = Path(__file__).resolve().parents[2]


def build_driver(build: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    objects = []
    for name, source in [('kinds', 'src/kinds.f90'), ('sparse', 'src/sparse.f90'),
                         ('lsmr', 'src/optimize/lsmr.f90'),
                         ('sparse_linalg', 'src/sparse_linalg.f90')]:
        obj = build / (name + '.o')
        subprocess.run([*cc, '-c', str(ROOT / source), '-o', str(obj)],
                       check=True, cwd=build, capture_output=True)
        objects.append(obj)
    binary = build / 'sparse_bicg_kronsum_reference_driver'
    subprocess.run([*cc, str(ROOT / 'tools/sparse_reference/sparse_bicg_kronsum_reference_driver.f90'),
                    *(str(obj) for obj in objects), '-o', str(binary)],
                   check=True, cwd=build, capture_output=True)
    return binary


def value_lines(values):
    return ' '.join(f'{v:.17e}' for v in values)


def run_kronsum(exe: Path, rng: np.random.Generator, count: int) -> None:
    cases, lines = [], []
    for _ in range(count):
        m, n = int(rng.integers(0, 8)), int(rng.integers(0, 8))
        def rand_square(d: int):
            rr, cc = np.nonzero(rng.random((d, d)) < 0.32)
            vv = rng.integers(-6, 7, size=len(rr)).astype(float)
            return rr, cc, vv
        ra, ca, va = rand_square(m)
        rb, cb, vb = rand_square(n)
        a = sparse.coo_array((va, (ra, ca)), shape=(m, m)).tocsr()
        b = sparse.coo_array((vb, (rb, cb)), shape=(n, n)).tocsr()
        cases.append((a, b))
        lines.extend(('K', f'{m} {n} {len(ra)} {len(rb)}'))
        lines.extend(f'{r+1} {c+1} {v:.17e}' for r,c,v in zip(ra,ca,va))
        lines.extend(f'{r+1} {c+1} {v:.17e}' for r,c,v in zip(rb,cb,vb))
    result = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                            text=True, capture_output=True, check=True)
    output = iter(result.stdout.splitlines())
    for i, (a, b) in enumerate(cases):
        fields = next(output).split()
        assert fields[0] == 'OK', f'case {i}: {fields}'
        nr, nc, nz = map(int, fields[1:])
        triplets = np.array([np.fromstring(next(output), sep=' ') for _ in range(nz)])
        if nz:
            actual = sparse.coo_array((triplets[:,2],
                       (triplets[:,0].astype(int)-1,triplets[:,1].astype(int)-1)),
                                      shape=(nr,nc)).toarray()
        else:
            actual = np.zeros((nr,nc))
        ref = sparse.kronsum(a, b, format='csr').toarray()
        assert actual.shape == ref.shape, f'case {i}: shape'
        np.testing.assert_allclose(actual, ref, atol=0, rtol=0,
                                   err_msg=f'kronsum #{i}')
    assert next(output, None) is None
    print(f'kronsum: {count} randomized SciPy comparisons PASS')


def run_solvers(exe: Path, rng: np.random.Generator, count: int) -> None:
    cases, lines = [], []
    for i in range(count):
        n = int(rng.integers(2, 40))
        a = rng.normal(size=(n,n))*(rng.random((n,n)) < rng.uniform(0.06,0.28))
        np.fill_diagonal(a, np.sum(np.abs(a),axis=1) + rng.uniform(0.8,2.0,size=n))
        b = rng.normal(size=n)
        x0 = rng.normal(size=n) if i%3 == 0 else np.zeros(n)
        method = 'B' if i%2 else 'C'
        jac = int(i%4 >= 2)
        rr, cc = np.nonzero(a)
        va = a[rr,cc]
        cases.append((a,b,x0,method,jac))
        lines.extend(('S', f'{n} {len(va)} {method} {jac}'))
        lines.extend(f'{r+1} {c+1} {v:.17e}' for r,c,v in zip(rr,cc,va))
        lines.append(value_lines(b))
        lines.append(value_lines(x0))
    result = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                            text=True, capture_output=True, check=True)
    output = iter(result.stdout.splitlines())
    for i, (a,b,x0,method,jac) in enumerate(cases):
        fields = next(output).split()
        assert fields[0] == 'SOLVE', f'case {i}: {fields}'
        status, iters = int(fields[1]),int(fields[2])
        reported = float(fields[3])
        sol = np.fromstring(next(output),sep=' ')
        M = sparse.diags(1/np.diag(a)) if jac else None
        func = spla.bicg if method == 'B' else spla.cgs
        ref, status_s = func(sparse.csr_array(a),b,x0=x0,rtol=1e-10,atol=1e-12,maxiter=300,M=M)
        error = np.linalg.norm(b-a@sol)
        assert status_s == 0, f'SciPy did not converge #{i}: {status_s}'
        assert status == 0 and 0 < iters <= 300, f'{method} #{i}: {status} {iters} {error}'
        assert error <= max(1e-10*np.linalg.norm(b),1e-12)*1.05, f'residual #{i}: {error}'
        np.testing.assert_allclose(reported,error,rtol=2e-4,atol=5e-14,
                                   err_msg=f'reported residual #{i}')
        np.testing.assert_allclose(sol,ref,rtol=5e-9,atol=1e-9,
                                   err_msg=f'{method} against SciPy #{i}')
    assert next(output, None) is None
    print(f'BiCG/CGS: {count} randomized SciPy comparisons PASS')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'),default='both')
    parser.add_argument('--kronsum-count',type=int,default=180)
    parser.add_argument('--solver-count',type=int,default=120)
    args = parser.parse_args()
    modes = ['checked','optimized'] if args.mode == 'both' else [args.mode]
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-bicg-kronsum-') as path:
            exe = build_driver(Path(path),mode)
            rng = np.random.default_rng(20261105)
            print(mode, end=': ', flush=True)
            run_kronsum(exe,rng,args.kronsum_count)
            print(mode, end=': ', flush=True)
            run_solvers(exe,rng,args.solver_count)


if __name__ == '__main__':
    main()
