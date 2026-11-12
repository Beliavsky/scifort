#!/usr/bin/env python3
"""Check SciFort real Schur decomposition against independent SciPy results.

Compare similarity, orthogonality, real quasi-triangular structure, eigenvalue
multisets, and whole-conjugate-pair selection, not individual Schur vectors.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_matrix_functions.f90",
           "src/linalg_schur.f90"]
RULES = ("none", "lhp", "rhp", "iuc", "ouc", "positive-imag")


def build(work: Path, mode: str, source: str, stem: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = work / f"real_schur_{stem}_{mode}"
    cmd = ["gfortran", "-std=f2018", *flags, *(str(ROOT / x) for x in SOURCES),
           str(ROOT / source), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=120)
    if p.returncode:
        raise RuntimeError(f"compile failed ({mode}):\n{p.stderr}")
    return exe


def matrices(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = [0, 1, 2, 3, 4, 6, 8, 12, max_order]
    for k in range(count):
        n = sizes[k % len(sizes)]
        a = rng.standard_normal((n, n))
        family = (k // len(sizes)) % 13
        if family == 1:
            a = np.triu(a)
        elif family == 2:
            a = np.tril(a)
        elif family == 3:
            a = np.diag(np.linspace(-2.5, 2.5, n))
        elif family == 4:
            a = a + a.T
        elif family == 5:
            a = a - a.T
        elif family == 6:
            a *= 1e110
        elif family == 7:
            a *= 1e-110
        elif family == 8:
            a = np.eye(n) * 2.3
        elif family == 9:
            a = np.triu(a, 1)
        elif family == 10:
            if n:
                b = rng.standard_normal((n, n))
                a = np.linalg.solve(b, np.diag(np.arange(n)) @ b)
        elif family == 11:
            a = np.diag(np.repeat(np.arange((n+1)//2), 2)[:n]) + np.triu(a, 1)*1e-3
        elif family == 12:
            a = 2.3*np.eye(n) + np.triu(a, 1)*0.1
        rule = RULES[(k // len(sizes)) % len(RULES)]
        if n:
            eig = linalg.eigvals(a)
            scale = max(float(np.max(np.abs(eig))), 1e-300)
            if rule in ('lhp', 'rhp') and np.any(np.abs(eig.real) <= 3e-10 * scale):
                rule = 'none'
            if rule in ('iuc', 'ouc') and np.any(np.abs(np.abs(eig)-1) <= 3e-10):
                rule = 'none'
            if rule == 'positive-imag' and np.any(
                    (np.abs(eig.imag) > 3e-10 * scale) &
                    (np.abs(eig.imag) < 1e-7 * scale)):
                rule = 'none'
        yield a, rule


def compare(exe: Path, cases, tolerance: float):
    input_lines = [str(len(cases))]
    for a, rule in cases:
        input_lines.append(f"{len(a)} {RULES.index(rule)}")
        input_lines.extend(f"{x:.17e}" for x in a.ravel(order="F"))
    env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'}
    p = subprocess.run([str(exe)], input='\n'.join(input_lines)+'\n', text=True,
                       capture_output=True, timeout=180, env=env)
    if p.returncode:
        raise RuntimeError(f"Fortran driver failure:\n{p.stderr[-3000:]}")
    lines = iter(p.stdout.splitlines())
    worst = dict(reconstruction=0., orthogonality=0., structure=0., eigenvalues=0.)
    sorted_cases = 0
    for case_index, (a, rule) in enumerate(cases):
        n = len(a)
        status, sdim = [int(x) for x in next(lines).split()]
        if status:
            raise AssertionError(f"case {case_index} n={n} family={(case_index // 9) % 13} "
                                 f"rule={rule}: returned status={status}")
        if not n:
            if sdim: raise AssertionError('nonzero selected dimension for empty matrix')
            continue
        output = np.asarray([[float(v) for v in next(lines).split()]
                             for _ in range(n*n)])
        t = output[:, 0].reshape((n,n), order='F')
        z = output[:, 1].reshape((n,n), order='F')
        scale = max(float(np.max(np.abs(a))), 1e-300)
        denominator = max(float(linalg.norm(a / scale)), 1.)
        reconstruction = float(linalg.norm(z@(t / scale)@z.T - a / scale))/denominator
        orthogonality = float(linalg.norm(z.T@z - np.eye(n))) / n
        structure = float(linalg.norm(np.tril(t / scale, -2))) / denominator
        for idx in range(n-2):
            if abs(t[idx+1,idx]) > tolerance*scale and abs(t[idx+2,idx+1]) > tolerance*scale:
                raise AssertionError(f"case {case_index} has overlapping 2x2 blocks")
        eig = np.asarray(linalg.eigvals(t))
        reference_eig = linalg.eigvals(a)
        cost = np.abs(eig[:,None] / scale - reference_eig[None,:] / scale)
        rr,cc = linear_sum_assignment(cost)
        eigenvalues = float(cost[rr,cc].max()) / denominator
        if rule != 'none':
            reference_sort = (lambda real, imag: imag > 1e-12*scale) if rule == 'positive-imag' else rule
            _,_,reference_sdim = linalg.schur(a,output='real',sort=reference_sort)
            if sdim != reference_sdim:
                raise AssertionError(f"case {case_index} sort={rule}: sdim={sdim} "
                                     f"but SciPy={reference_sdim}")
            sorted_cases+=1
            if sdim:
                invariant = float(linalg.norm((a/scale)@z[:,:sdim] -
                                              z[:,:sdim]@(t[:sdim,:sdim]/scale))) / denominator
                if invariant > tolerance:
                    raise AssertionError(f"case {case_index} invariant residual={invariant}")
        for key,val in dict(reconstruction=reconstruction, orthogonality=orthogonality,
                            structure=structure, eigenvalues=eigenvalues).items():
            worst[key] = max(worst[key], val)
            if not np.isfinite(val) or val > tolerance:
                raise AssertionError(f"case {case_index} n={n} rule={rule} {key}={val:.3g}")
    return worst, sorted_cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['checked','optimized','both'], default='both')
    parser.add_argument('--cases', type=int, default=250)
    parser.add_argument('--max-order', type=int, default=40)
    parser.add_argument('--seed', type=int, default=20261108)
    parser.add_argument('--tolerance', type=float, default=3e-8)
    args = parser.parse_args()
    if not 2 <= args.max_order <= 40:
        parser.error('max-order must lie between 2 and 40')
    cases = list(matrices(args.cases,args.max_order,args.seed))
    with tempfile.TemporaryDirectory(prefix='scifort_real_schur_') as dirname:
        work = Path(dirname)
        for mode in (['checked','optimized'] if args.mode == 'both' else [args.mode]):
            native = build(work,mode,'test/test_linalg_real_schur.f90','native')
            p = subprocess.run([str(native)],cwd=work,capture_output=True,text=True,timeout=60)
            if p.returncode: raise RuntimeError(f'Native {mode}: {p.stdout}\n{p.stderr}')
            print(f'{mode}: {p.stdout.strip()}')
            exe = build(work,mode,'tools/linalg_reference/driver_real_schur.f90','driver')
            worst, sorted_cases = compare(exe,cases,args.tolerance)
            print(f'{mode}: {len(cases)} SciPy comparisons passed ({sorted_cases} sorted); '
                  + ', '.join(f'{k}={v:.3e}' for k,v in worst.items()))

if __name__ == '__main__':
    main()
