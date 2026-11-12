#!/usr/bin/env python3
"""Independent SciPy/NumPy checks for CSR LGMRES and GCROT(m,k).

Runs two solves per test system, with reusable augmentation vectors for the
second RHS. True residuals, dense-reference solutions, SciPy solutions, and
all solver return codes are checked. No BLAS/LAPACK dependencies in Fortran.
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
    flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace',
              '-ffpe-trap=invalid,zero,overflow'] if mode == 'checked'
             else ['-O3'])
    command = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
    objects = []
    for name, source in [('kinds', 'src/kinds.f90'), ('sparse', 'src/sparse.f90'),
                         ('lsmr', 'src/optimize/lsmr.f90'),
                         ('sparse_linalg', 'src/sparse_linalg.f90')]:
        obj = build / (name + '.o')
        subprocess.run([*command, '-c', str(ROOT / source), '-o', str(obj)],
                       cwd=build, check=True, capture_output=True)
        objects.append(obj)
    exe = build / 'sparse_recycling_reference_driver'
    subprocess.run([*command, str(ROOT / 'tools/sparse_reference/'
                    'sparse_recycling_reference_driver.f90'),
                    *(str(obj) for obj in objects), '-o', str(exe)],
                   cwd=build, check=True, capture_output=True)
    return exe


def random_system(rng: np.random.Generator, n: int, index: int) -> np.ndarray:
    # Mild nonsymmetric systems, with diagonal scaling and sparse banded cases.
    a = rng.normal(size=(n, n)) * (rng.random((n, n)) < rng.uniform(0.08, 0.28))
    if index % 5 == 0:
        a = np.triu(a, -2)
        a += np.diag(rng.uniform(-1, 1, n-1), 1)
    np.fill_diagonal(a, 0)
    diagonal = np.sum(abs(a), axis=1) + rng.uniform(0.5, 1.5, n)
    if index % 4 == 1:
        scale = np.exp(rng.uniform(-0.8, 0.8, n))
        a *= scale[:, None]
        diagonal *= scale
    np.fill_diagonal(a, diagonal)
    return a


def run_cases(exe: Path, count: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    cases, lines = [], []
    for i in range(count):
        n = int(rng.integers(2, 46))
        a = random_system(rng, n, i)
        x1, x2 = rng.normal(size=(2, n))
        # Related right sides; recycled augmentation should be useful here.
        b1, b2 = a @ x1, a @ (x1 + 0.1*x2)
        x01 = x1 if i % 7 == 0 else rng.normal(size=n) if i % 3 == 0 else np.zeros(n)
        x02 = rng.normal(size=n) if i % 11 == 0 else np.zeros(n)
        code = 'L' if i % 2 == 0 else 'G'
        pre = (i // 2) % 3
        m = int(rng.integers(3, min(13, n+1))) if n >= 4 else n
        k = int(rng.integers(0, min(4, n)))
        rr, cc = np.nonzero(a)
        va = a[rr,cc]
        lines.extend((code, f'{n} {len(va)} {pre} {m} {k}'))
        lines.extend(f'{r+1} {c+1} {v:.17e}' for r,c,v in zip(rr,cc,va))
        for v in (b1,b2,x01,x02):
            lines.append(' '.join(f'{t:.17e}' for t in v))
        cases.append((a,code,pre,m,k,(b1,b2),(x01,x02)))
    result = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',
                            text=True, capture_output=True,check=True)
    output = iter(result.stdout.splitlines())
    successes = 0
    for i,(a,code,pre,m,k,bs,xs) in enumerate(cases):
        for j,(b,x0) in enumerate(zip(bs,xs)):
            fields = next(output).split()
            assert fields[0]=='SOLVE', (i,j,fields)
            status, steps = int(fields[1]),int(fields[2])
            reported, nkept = float(fields[3]),int(fields[4])
            actual = np.fromstring(next(output),sep=' ')
            assert actual.shape==b.shape,(i,j,actual)
            err = np.linalg.norm(b-a@actual)
            target = max(1e-12,1e-10*np.linalg.norm(b))
            assert status == 0, (i,j,code,pre,m,k,status,steps,err,target)
            assert 0 <= nkept <= k,(i,j,k,nkept)
            assert steps >= 0 and steps <= 180*(m+k+3),(i,j,steps)
            assert err <= target*1.002+2e-13,(i,j,code,pre,err,target)
            assert abs(err-reported) <= max(1e-11,1e-4*err),(i,j,err,reported)
            truth = np.linalg.solve(a,b)
            np.testing.assert_allclose(actual,truth,atol=2e-8,rtol=2e-8,
                                       err_msg=f'{code} case {i} RHS {j}')
            # Use SciPy with matching (unpreconditioned) A and tolerances.
            a_sp = sparse.csr_array(a)
            if code == 'L':
                scipy_x, scipy_status = spla.lgmres(a_sp,b,x0=x0,inner_m=m,outer_k=k,
                                                    maxiter=180,rtol=1e-10,atol=1e-12)
            else:
                scipy_x, scipy_status = spla.gcrotmk(a_sp,b,x0=x0,m=m,k=k,
                                                     maxiter=180,rtol=1e-10,atol=1e-12)
            assert scipy_status == 0,(i,j,code,scipy_status)
            np.testing.assert_allclose(actual,scipy_x,atol=2e-8,rtol=2e-8,
                                       err_msg=f'{code} SciPy case {i} RHS {j}')
            successes += 1
    assert next(output,None) is None
    print(f'{count} systems x 2 RHS = {successes} SciPy comparisons PASS')


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--systems',type=int,default=150)
    args=parser.parse_args()
    for mode in ('checked','optimized') if args.mode=='both' else (args.mode,):
        with tempfile.TemporaryDirectory(prefix='scifort-recycling-') as temp:
            exe=build_driver(Path(temp),mode)
            print(f'{mode}: ',end='',flush=True)
            run_cases(exe,args.systems,20261105)


if __name__=='__main__':
    main()
