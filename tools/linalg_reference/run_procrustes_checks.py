#!/usr/bin/env python3
"""Compile and compare SciFort's orthogonal Procrustes solution with SciPy.

Checks real/complex/mixed arrays, rank deficiency, empty arrays, orthogonality,
objective optimality and returned nuclear norms. Requires NumPy and SciPy.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    'src/kinds.f90', 'src/linalg.f90', 'src/linalg_svd.f90',
    'src/linalg_procrustes.f90',
)


def make_cases(count: int, seed: int, max_order: int):
    rng = np.random.default_rng(seed)
    sizes = sorted(set([0, 1, 2, 3, 4, 5, 8, 12, 17, max_order]))
    cases = []
    for idx in range(count):
        m = int(rng.choice(sizes))
        n = int(rng.choice(sizes))
        kind_case = 1 + idx % 4
        family = idx % 13
        def random_array(shape, complex_output):
            z = rng.standard_normal(shape)
            if complex_output:
                z = z + 1j*rng.standard_normal(shape)
            return z

        a = random_array((m,n),kind_case in (2,4))
        b = random_array((m,n),kind_case in (2,3))
        if family in (1,4,8) and m and n:
            z = random_array((n,n),kind_case in (2,3))
            q, _ = linalg.qr(z)
            if kind_case in (1,2,3):
                b = a@q
        if family == 2 and n >= 2:
            a[:,-1] = a[:,0]
            b[:,-1] = b[:,0]
        if family == 3:
            a[:] = 0
        if family == 5:
            b[:] = 0
        if family == 6 and m >= 2*n and n:
            a[:] = 0
            b[:] = 0
            for j in range(n):
                a[j,j] = 1
                b[n+j,j] = 1
        if family == 7 and n >= 2:
            a[:,-1] = a[:,0]*(1+1e-10)
        if family == 9 and m and n:
            a[:] = 1
            b[:] = -1
        if family == 10:
            a *= 1e140
            b *= 1e-140
        if family == 11:
            a *= 1e140
            b *= 1e140
        if family == 12:
            a *= 1e-140
            b *= 1e-140
        cases.append((kind_case, np.asarray(a), np.asarray(b)))
    return cases


def compile_driver(tmp: Path, mode: str) -> Path:
    flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
    exe = tmp/('procrustes_'+mode)
    command = ['gfortran','-std=f2018',*flags,*(str(ROOT/p) for p in SOURCES),
               str(ROOT/'tools/linalg_reference/driver_procrustes.f90'),'-o',str(exe)]
    proc = subprocess.run(command,cwd=tmp,capture_output=True,text=True,check=False)
    if proc.returncode:
        raise RuntimeError('driver compile failed:\n'+proc.stderr)
    return exe


def run_cases(exe: Path, cases, tolerance: float):
    lines = [str(len(cases))]
    for kind_case,a,b in cases:
        m,n = a.shape
        lines.append(f'{kind_case} {m} {n}')
        for av,bv in zip(a.ravel(order='F'),b.ravel(order='F')):
            lines.append(f'{av.real:.17e} {av.imag:.17e} {bv.real:.17e} {bv.imag:.17e}')
    proc = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,
                          capture_output=True,timeout=200,check=False)
    if proc.returncode:
        raise RuntimeError('Fortran driver failed:\n'+proc.stderr[-2500:])
    output = iter(proc.stdout.splitlines())
    worst = dict(scale=0.,unitarity=0.,objective=0.,rotation=0.)
    for idx, (_,a,b) in enumerate(cases):
        m,n = a.shape
        stat_str,scale_str = next(output).split()
        stat = int(stat_str)
        scale = float(scale_str)
        if stat:
            raise AssertionError(f'case {idx} returned status {stat}, shape {(m,n)}')
        r = np.array([complex(*map(float,next(output).split()))
                      for _ in range(n*n)]).reshape((n,n),order='F')
        reference, ref_scale = linalg.orthogonal_procrustes(a,b)
        # Non-unique minimizers must be compared by their objectives, not R.
        unit_err = linalg.norm(r.conj().T@r-np.eye(n),'fro')/max(1,n)
        worst['unitarity'] = max(worst['unitarity'],unit_err)
        if unit_err > 2*tolerance:
            raise AssertionError(f'case {idx}: nonunitary solution, {unit_err}')
        common = max(float(np.max(np.abs(a))) if a.size else 0.,
                     float(np.max(np.abs(b))) if b.size else 0.)
        if common:
            aa,bb = a/common,b/common
            objective = linalg.norm(aa@r-bb,'fro')
            objective_ref = linalg.norm(aa@reference-bb,'fro')
            objective_err = max(0.,objective-objective_ref)
            worst['objective'] = max(worst['objective'],objective_err)
            if objective_err > tolerance*max(1.,linalg.norm(aa,'fro'),linalg.norm(bb,'fro')):
                raise AssertionError(f'case {idx}: objective {objective} worse than {objective_ref}')
        denom = max(abs(ref_scale),np.finfo(float).tiny)
        scale_err = abs(scale-ref_scale)/denom
        worst['scale'] = max(worst['scale'],scale_err)
        if scale_err > 4*tolerance:
            raise AssertionError(f'case {idx}: scale {scale} != {ref_scale}, error {scale_err}')
        if m >= n and n and idx%13 in (1,4,8) and np.linalg.matrix_rank(a) == n:
            # The polar factor is unique for invertible cross covariance.
            rotation_err = linalg.norm(r-reference,'fro')/max(1,np.sqrt(n))
            worst['rotation'] = max(worst['rotation'],rotation_err)
            if rotation_err > 20*tolerance:
                raise AssertionError(f'case {idx}: rotation difference {rotation_err}')
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    parser.add_argument('--cases',type=int,default=300)
    parser.add_argument('--seed',type=int,default=1117)
    parser.add_argument('--max-order',type=int,default=24)
    parser.add_argument('--tol',type=float,default=2e-8)
    args = parser.parse_args()
    if args.cases < 1 or args.max_order < 1:
        parser.error('--cases and --max-order must be positive')
    cases = make_cases(args.cases,args.seed,args.max_order)
    with tempfile.TemporaryDirectory(prefix='scifort_procrustes_') as directory:
        tmp = Path(directory)
        modes = ['checked','optimized'] if args.mode == 'both' else [args.mode]
        for mode in modes:
            exe = compile_driver(tmp,mode)
            worst = run_cases(exe,cases,args.tol)
            print(f'{mode}: {len(cases)} SciPy comparisons passed; '
                  +', '.join(f'{key}={value:.3e}' for key,value in worst.items()))


if __name__ == '__main__':
    main()
