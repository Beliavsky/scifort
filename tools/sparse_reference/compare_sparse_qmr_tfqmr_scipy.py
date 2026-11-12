#!/usr/bin/env python3
"""Check SciFort CSR QMR and TFQMR against SciPy and dense reference solves.

Use SciPy's qmr for QMR, including right Jacobi preconditioning. SciFort TFQMR
solves the consistent right-preconditioned A*P problem; reference SciPy TFQMR
therefore runs on A*P rather than supplying its M argument, which follows a
different recurrence for nonidentity preconditioners. Every case additionally
verifies the true original-system residual, solution and reported diagnostic.
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
    exe = build / 'sparse_qmr_tfqmr_reference_driver'
    subprocess.run([*command, str(ROOT / 'tools/sparse_reference/'
                    'sparse_qmr_tfqmr_reference_driver.f90'),
                    *(str(obj) for obj in objects), '-o', str(exe)],
                   cwd=build, check=True, capture_output=True)
    return exe


def random_system(rng: np.random.Generator, n: int, index: int) -> np.ndarray:
    a = rng.normal(size=(n, n)) * (rng.random((n, n)) < rng.uniform(0.09, 0.35))
    np.fill_diagonal(a, 0)
    diagonal = np.sum(np.abs(a), axis=1) + rng.uniform(0.5, 2.0, n)
    if index % 4 == 3:
        # Some widely varying diagonal scales test right-preconditioning.
        scale = np.exp(rng.uniform(-1.0, 1.0, n))
        a *= scale[:, None]
        diagonal *= scale
    np.fill_diagonal(a, diagonal)
    return a


def run_cases(exe: Path, count: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    cases, lines = [], []
    for i in range(count):
        n = int(rng.integers(2, 41))
        a = random_system(rng, n, i)
        b = rng.normal(size=n)
        x0 = rng.normal(size=n) if i % 3 == 0 else np.zeros(n)
        method = 'Q' if i % 2 else 'T'
        pre = (i // 2) % (2 if method == 'Q' else 3)
        rr, cc = np.nonzero(a)
        va = a[rr, cc]
        lines.extend((method, f'{n} {len(va)} {pre}'))
        lines.extend(f'{row+1} {col+1} {value:.17e}'
                     for row, col, value in zip(rr, cc, va))
        lines.append(' '.join(f'{v:.17e}' for v in b))
        lines.append(' '.join(f'{v:.17e}' for v in x0))
        cases.append((a, b, x0, method, pre))
    proc = subprocess.run([str(exe)], input='\n'.join(lines)+'\n', text=True,
                          capture_output=True, check=True)
    output = iter(proc.stdout.splitlines())
    tol = 1e-10
    for i, (a,b,x0,method,pre) in enumerate(cases):
        fields = next(output).split()
        assert fields[0] == 'SOLVE', f'case {i}: {fields}'
        info, iters = int(fields[1]), int(fields[2])
        reported = float(fields[3])
        actual = np.fromstring(next(output), sep=' ')
        assert len(actual) == len(b), f'case {i}: wrong solution dimension'
        actual_resid = np.linalg.norm(b-a @ actual)
        norm_target = max(1e-12, tol*np.linalg.norm(b))
        assert info == 0, (i, method, pre, info, iters, actual_resid, norm_target)
        assert 0 < iters <= 350, (i, method, pre, iters)
        assert actual_resid <= norm_target*1.001 + 1e-14, (i,method,pre,actual_resid)
        assert abs(reported-actual_resid) <= max(5e-13,1e-5*actual_resid), (i,reported,actual_resid)
        reference = np.linalg.solve(a,b)
        np.testing.assert_allclose(actual,reference,rtol=2e-8,atol=2e-9,
                                   err_msg=f'{method} case {i}')
        # Independent SciPy solver with mathematically equivalent operator.
        if method == 'Q':
            m2 = spla.aslinearoperator(sparse.diags(1.0/np.diag(a))) if pre else None
            reference_scipy, scipy_status = spla.qmr(sparse.csr_array(a),b,x0=x0,
                                                     M1=spla.aslinearoperator(sparse.eye(len(b),format="csr")) if pre else None,
                                                     M2=m2,rtol=tol,atol=1e-12,
                                                     maxiter=350)
        elif pre == 1:
            p = 1.0/np.diag(a)
            y, scipy_status = spla.tfqmr(sparse.csr_array(a*p[np.newaxis, :]),
                                         b-a @ x0,rtol=tol,atol=1e-12,maxiter=350)
            reference_scipy = x0+p*y
        else:
            reference_scipy, scipy_status = spla.tfqmr(sparse.csr_array(a),b,
                                                      x0=x0,rtol=tol,atol=1e-12,
                                                      maxiter=350)
        if pre != 2:
            assert scipy_status == 0, (i,method,pre,'SciPy did not converge',scipy_status)
            np.testing.assert_allclose(actual,reference_scipy,rtol=4e-8,atol=4e-9,
                                       err_msg=f'{method} SciPy case {i}')
    assert next(output,None) is None
    print(f'{count} randomized QMR/TFQMR comparisons PASS')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'),default='both')
    parser.add_argument('--solver-count',type=int,default=180)
    args = parser.parse_args()
    modes = ['checked','optimized'] if args.mode == 'both' else [args.mode]
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-qmr-tfqmr-') as path:
            exe = build_driver(Path(path),mode)
            print(mode, end=': ', flush=True)
            run_cases(exe,args.solver_count,20261106)


if __name__ == '__main__':
    main()
