#!/usr/bin/env python3
"""Compare SciFort's CSR norm estimator and expm against SciPy/NumPy.

The Fortran executable handles every randomized case in a single process.
Checks the one-norm lower bound and its explicit witness, along with a full
matrix exponential against both scipy.sparse.linalg.expm and scipy.linalg.expm.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy.linalg import expm as dense_expm
from scipy.sparse import csr_array
from scipy.sparse.linalg import expm as sparse_expm, onenormest as scipy_onenormest

ROOT = Path(__file__).resolve().parents[2]


def make_case(rng: np.random.Generator, i: int):
    n = int(rng.integers(1, 18))
    a = rng.normal(0, 1, size=(n, n)) * (rng.random((n, n)) < .38)
    scale = float(rng.choice([.02, .1, .5, 1., 2.]))
    a *= scale
    if i % 12 == 0:
        a = np.diag(rng.uniform(-15, 15, n))
    elif i % 12 == 1:
        a += np.eye(n) * rng.uniform(-13, 13)
    elif i % 12 == 2:
        a = np.triu(a * 4, k=1)  # nilpotent
    elif i % 12 == 3:
        a = a - a.T  # rotation
    elif i % 12 == 4:
        a = (a + a.T) / 2
    elif i % 12 == 5:
        a.fill(0.)
    elif i % 12 == 6:
        a -= np.eye(n) * 5
    elif i % 12 == 7:
        a = np.diag(rng.uniform(-1, 1, n))
        if n > 1:
            a += np.diag(rng.uniform(-8, 8, n-1), k=1)
    elif i % 12 == 8:
        a[:, 0] = a[:, -1]  # rank deficient
    elif i % 12 == 9:
        a = a @ a.T * .3
    elif i % 12 == 10:
        a *= .0001

    width = int(rng.integers(1, min(n+2, 7)))
    if i % 17 == 0:
        width = n
    itmax = int(rng.integers(2, 9))
    chunk = int(rng.integers(1, min(n+3, 6)))
    # COO will occasionally contain duplicate or explicitly zero elements.
    rr, cc = np.nonzero(a)
    vv = a[rr, cc].copy()
    rr = rr.tolist()
    cc = cc.tolist()
    vv = vv.tolist()
    if len(rr) and i % 7 == 0:
        index = int(rng.integers(len(rr)))
        half = vv[index] * .34
        vv[index] -= half
        rr.append(rr[index])
        cc.append(cc[index])
        vv.append(half)
    if i % 6 == 0:
        rr.append(int(rng.integers(n)))
        cc.append(int(rng.integers(n)))
        vv.append(0.)
    a = csr_array((vv, (rr, cc)), shape=(n,n)).toarray()
    return a, rr, cc, vv, width, itmax, chunk


def run_suite(mode: str, cases: list, verbose: bool = False) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-onenorm-expm-') as td:
        build = Path(td)
        flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
        compiler = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
        objects = []
        for short_name, filename in (
            ('kinds', 'src/kinds.f90'),
            ('sparse', 'src/sparse.f90'),
            ('lsmr', 'src/optimize/lsmr.f90'),
            ('sparse_linalg', 'src/sparse_linalg.f90'),
        ):
            obj = build / (short_name+'.o')
            subprocess.run([*compiler, '-c', str(ROOT/filename), '-o', str(obj)],
                           cwd=build, check=True)
            objects.append(str(obj))
        exe = build/'reference_driver'
        subprocess.run([*compiler, str(ROOT/'tools/sparse_linalg_reference'/
                             'sparse_onenorm_expm_reference_driver.f90'),
                        *objects, '-o', str(exe)], cwd=build, check=True)

        input_lines = [str(len(cases))]
        for a, rows, columns, values, width, limit, chunk in cases:
            # COO can be empty even for nonempty matrices.
            n = a.shape[0]
            input_lines.append(f'{n} {len(rows)} {width} {limit} {chunk}')
            input_lines.extend(f'{r+1} {c+1} {float(v):.17g}'
                               for r, c, v in zip(rows, columns, values))
        proc = subprocess.run([str(exe)], cwd=build,
                              input='\n'.join(input_lines)+'\n', text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              timeout=180, check=True)
        lines = iter(proc.stdout.splitlines())
        worst_expm = 0.
        worst_est_ratio = 1.
        full_width = 0
        scipy_norm_ratio = 1.0
        for i, (a, _, _, _, width, limit, _) in enumerate(cases):
            n = a.shape[0]
            info_n, info_e, est = next(lines).split()
            info_n, info_e, est = int(info_n), int(info_e), float(est)
            assert info_n == 0 and info_e == 0, (mode, i, info_n, info_e)
            v = np.fromstring(next(lines), sep=' ')
            w = np.fromstring(next(lines), sep=' ')
            output = np.fromstring(next(lines), sep=' ').reshape((n, n), order='F')
            assert v.shape == (n,) and w.shape == (n,)
            exact_norm = np.linalg.norm(a, ord=1)
            reference_est = float(scipy_onenormest(csr_array(a),t=min(width,n),itmax=limit))
            assert 0 <= reference_est <= exact_norm*(1+2e-12)+1e-12, (mode,i,'SciPy norm')
            if exact_norm > 0:
                scipy_norm_ratio = min(scipy_norm_ratio, reference_est/exact_norm)
            assert -1e-10 <= est <= exact_norm*(1+2e-12)+1e-12, \
                (mode, i, est, exact_norm)
            assert abs(np.linalg.norm(v, ord=1) - 1) <= 1e-12, (mode, i)
            assert np.allclose(w, a @ v, rtol=2e-12, atol=2e-12), (mode, i)
            assert abs(np.linalg.norm(w, ord=1)-est) <= 2e-12*max(1,est), (mode, i)
            if n <= width:
                full_width += 1
                assert abs(est-exact_norm) <= 4e-12*max(1,exact_norm), (mode,i,est,exact_norm)
            if exact_norm > 0:
                worst_est_ratio = min(worst_est_ratio, est/exact_norm)
            reference = dense_expm(a)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                scipy_result = sparse_expm(csr_array(a)).toarray()
            ref_scale = max(1.,np.max(np.abs(reference)))
            scipy_err = np.max(np.abs(scipy_result-reference)) / ref_scale
            assert scipy_err <= 5e-10, (mode,i,'SciPy reference',scipy_err)
            err = np.max(np.abs(output-reference)) / ref_scale
            worst_expm = max(worst_expm, float(err))
            assert err <= 5e-10, (mode, i, err)
            if verbose and i % 50 == 0:
                print(f'{mode} {i}: est={est:g} exact={exact_norm:g} exp_err={err:.3g}', flush=True)
        assert next(lines, None) is None
        print(f'{mode}: {len(cases)} randomized matrices PASSED; '
              f'{full_width} exact-width norm checks; worst estimate/exact '
              f'{worst_est_ratio:.3f} (SciPy worst ratio {scipy_norm_ratio:.3f}); '
              f'largest scaled expm error {worst_expm:.3e}', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    parser.add_argument('--cases', type=int, default=200)
    parser.add_argument('--seed', type=int, default=20261104)
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    cases = [make_case(rng, i) for i in range(args.cases)]
    for mode in (('checked','optimized') if args.mode == 'both' else (args.mode,)):
        run_suite(mode, cases, verbose=args.verbose)


if __name__ == '__main__':
    main()
