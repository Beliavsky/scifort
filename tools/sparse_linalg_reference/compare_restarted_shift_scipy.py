#!/usr/bin/env python3
"""Verify restarted shift-invert Lanczos/Arnoldi vs NumPy and SciPy.

Uses the full complex inverse for nonsymmetric problems, i.e. compares to
SciPy eigs on a COMPLEX cast of the real CSR operator; SciPy's real-matrix
OPpart='r' convention is a different transform for complex sigma.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs, eigsh

ROOT = Path(__file__).resolve().parents[2]


def make_case(rng, idx):
    sym = idx % 2 == 0
    n = (16, 22, 28, 34)[(idx // 2) % 4]
    k = 2 if sym else 1
    mode = ('LM', 'LA', 'SA')[idx // 8 % 3] if sym else 'LM'
    diag = rng.uniform(-4, 4, n)
    if sym:
        # Two isolated interior modes around a real shift, while the rest
        # lie away from zero.  Small couplings avoid exactly singular solves.
        diag = rng.choice([-1, 1], n) * rng.uniform(1.5, 5, n)
        diag[:4] = [-0.25, 0.35, -0.9, 1.05]
        sigma = float(rng.uniform(-0.08, 0.08))
        a = np.diag(diag)
        for i in range(n - 1):
            x = rng.uniform(-0.03, 0.03)
            a[i, i + 1] = a[i + 1, i] = x
        p = min(n, 12)
    else:
        # Isolated nonsymmetric complex pair; complex sigma selects +Im.
        a = np.diag(diag)
        a[0, 0] = a[1, 1] = 0.45
        a[0, 1], a[1, 0] = -1.05, 1.05
        for i in range(2, n - 1):
            a[i, i + 1] = rng.uniform(-.04, .04)
        sigma = complex(rng.uniform(.2, .4), rng.uniform(.75, .95))
        p = min(n, 10)
    return ('S' if sym else 'N', a, k, p, 70, mode, 2e-7, sigma)


def priority(mu, mode):
    return dict(LM=np.abs(mu), SM=-np.abs(mu), LA=mu.real,
                SA=-mu.real, LR=mu.real, SR=-mu.real,
                LI=mu.imag, SI=-mu.imag)[mode]


def compare_case(case, values, errors, vectors, status, cycles, locks, ops, case_id):
    kind, a, k, p, cap, mode, tol, sigma = case
    n = len(a)
    if status != 0:
        raise AssertionError(f'case {case_id}: nonconverged {kind} {mode} {n=} '
                             f'{status=} {cycles=} {errors=}')
    if cycles < 1 or locks != k or ops < p:
        raise AssertionError(f'case {case_id}: bad diagnostics {cycles=} {locks=} {ops=}')
    true_resid = np.linalg.norm(a @ vectors - vectors * values[None, :], axis=0)
    np.testing.assert_allclose(true_resid, errors, rtol=1e-5, atol=5e-11)
    if max(errors) > 1.1 * tol * max(1, np.linalg.norm(a, ord=np.inf), max(abs(values))):
        raise AssertionError(f'case {case_id}: excessive true residual {errors}')
    np.testing.assert_allclose(np.linalg.norm(vectors, axis=0), 1, atol=1e-8)
    dense = np.linalg.eigvalsh(a) if kind == 'S' else np.linalg.eigvals(a)
    mu = 1 / (dense - sigma)
    expected = dense[np.argsort(-priority(mu, mode))[:k]]
    rows, cols = linear_sum_assignment(abs(values[:, None] - expected[None, :]))
    if max(abs(values[rows]-expected[cols])) > 3e-5:
        raise AssertionError(f'case {case_id} {kind} {mode}: spectral mismatch '
                             f'{values=} {expected=}')
    # SciPy uses normal-mode inverse ranking for real symmetric problems.
    # Cast nonsymmetric A to complex to use full complex inverse, not OPpart='r'.
    if kind == 'S':
        ref = eigsh(csr_matrix(a), k=k, sigma=sigma, which=mode,
                    return_eigenvectors=False, tol=1e-11, maxiter=5000)
    else:
        ref = eigs(csr_matrix(a.astype(complex)), k=k, sigma=sigma,
                   which=mode, return_eigenvectors=False, tol=1e-11,
                   maxiter=5000)
    rows, cols = linear_sum_assignment(abs(values[:, None]-ref[None, :]))
    if max(abs(values[rows]-ref[cols])) > 3e-5:
        raise AssertionError(f'case {case_id}: SciPy mismatch {values=} {ref=}')
    return cycles


def run_mode(build_mode, count, seed):
    rng = np.random.default_rng(seed)
    cases = [make_case(rng, i) for i in range(count)]
    with tempfile.TemporaryDirectory(prefix='scifort-shift-restarted-') as folder:
        build = Path(folder)
        flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace']
                 if build_mode == 'checked' else ['-O3'])
        cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
        objs = []
        for name, path in [('kinds', 'src/kinds.f90'),
                           ('sparse', 'src/sparse.f90'),
                           ('lsmr', 'src/optimize/lsmr.f90'),
                           ('sparse_linalg', 'src/sparse_linalg.f90')]:
            obj = build / (name+'.o')
            subprocess.run([*cc,'-c',str(ROOT/path),'-o',str(obj)],cwd=build,check=True)
            objs.append(obj)
        exe = build / 'shift_restarted_reference'
        subprocess.run([*cc, str(ROOT/'tools/sparse_linalg_reference' /
                                  'sparse_restarted_shift_reference_driver.f90'),
                        *(str(o) for o in objs), '-o', str(exe)], cwd=build, check=True)
        lines = [str(len(cases))]
        for kind,a,k,p,cap,which,tol,sigma in cases:
            lines.append(f'{kind} {len(a)} {k} {p} {cap} {which} {tol:.17g} '
                         f'{sigma.real:.17g} {sigma.imag:.17g}')
            lines.append(' '.join(f'{v:.17g}' for v in a.flatten(order='F')))
        proc = subprocess.run([str(exe)],cwd=build,input='\n'.join(lines)+'\n',
                              capture_output=True,text=True,check=True,timeout=150)
        output = iter(proc.stdout.splitlines())
        restarted = 0
        for i,case in enumerate(cases,1):
            if next(output).strip() != f'CASE {i}':
                raise AssertionError(f'case output mismatch {i}')
            status, cycles, locks, ops = map(int, next(output).split())
            kind,a,k,*_ = case
            if status < 0:
                raise AssertionError(f'case {i}: breakdown status={status}')
            values, errors, vectors = [], [], []
            for j in range(k):
                data = np.fromstring(next(output),sep=' ')
                if kind == 'S':
                    assert data.size == 3+len(a)
                    vec = data[3:]
                else:
                    assert data.size == 3+2*len(a)
                    vec = data[3:3+len(a)] + 1j*data[3+len(a):]
                values.append(data[0]+1j*data[1])
                errors.append(data[2])
                vectors.append(vec)
            values = np.array(values)
            errors = np.array(errors)
            vectors = np.column_stack(vectors)
            restarted += compare_case(case,values,errors,vectors,status,
                                      cycles,locks,ops,i) > 1
        assert next(output,None) is None
        print(f'{build_mode}: {count} shifted restarted SciPy/NumPy comparisons passed; '
              f'{restarted} required restarts',flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    ap.add_argument('--cases',type=int,default=48)
    ap.add_argument('--seed',type=int,default=20261104)
    args = ap.parse_args()
    if args.cases < 1: ap.error('--cases must be positive')
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        run_mode(mode,args.cases,args.seed)

if __name__=='__main__': main()
