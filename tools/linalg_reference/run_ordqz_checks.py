#!/usr/bin/env python3
"""Compare SciFort ordered complex QZ against SciPy ordqz on regular pencils.

The two QZ decompositions are nonunique, so validation uses both scaled
reconstruction identities, unitarity, triangularity, projective eigenvalues,
and selected/unselected subspace ordering rather than identical factors.
Also tests independently supplied QZ factors without the SciFort QZ engine.
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy import linalg
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_dense.f90",
           "src/linalg_svd.f90", "src/linalg_qr.f90", "src/linalg_matrix_functions.f90",
           "src/linalg_schur.f90", "src/linalg_qz.f90", "src/linalg_ordqz.f90"]
SORTS = ("lhp", "rhp", "iuc", "ouc", "none")


def compile_driver(directory: Path, mode: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace",
              "-ffpe-trap=invalid,zero,overflow"] if mode == "checked" else ["-O3"])
    exe = directory / f"ordqz_{mode}"
    cmd = ["gfortran", "-std=f2018", *flags, *(str(ROOT / s) for s in SOURCES),
           str(ROOT / "tools/linalg_reference/driver_ordqz.f90"), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=directory, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"Compilation failed:\n{p.stderr}")
    return exe


def make_cases(cases: int, seed: int, max_order: int, factors_only: bool=False):
    rng = np.random.default_rng(seed)
    sizes = [0, 1, 2, 3, 4, 6, 8, min(12,max_order), max_order]
    for k in range(cases):
        n = sizes[k % len(sizes)]
        family = (k // len(sizes)) % 4 + 1
        sort = SORTS[(k // 2) % len(SORTS)]
        supplied = 1 if factors_only else (k // 3) % 2
        a = rng.standard_normal((n,n))
        b = rng.standard_normal((n,n))
        if family in (3,4):
            a = a + 1j*rng.standard_normal((n,n))
        if family in (2,4):
            b = b + 1j*rng.standard_normal((n,n))
        if n:
            b = b*0.15 + np.eye(n)
        pattern = (k // (len(sizes)*4)) % 11
        if pattern == 1:  # pretriangular, finite
            a = np.triu(a)
            b = np.triu(b)
        elif pattern == 2:  # pretriangular, selected infinite eigenvalues
            a = np.triu(a)
            b = np.triu(b)
            for j in range(0,n,3):
                b[j,j] = 0
                a[j,j] = 1 + 0.2j if family in (3,4) else 1
        elif pattern == 3:  # singular B, regular pencil
            b = np.diag(np.linspace(0.5,1.5,n)).astype(b.dtype)
            if n:
                for j in range(0,n,3):
                    b[j,j] = 0
        elif pattern == 4:  # repeated and opposite signs
            a = np.diag(np.linspace(-3.0,3.0,n)).astype(a.dtype)
            if n > 3:
                a[2,2] = a[1,1]
            b = np.eye(n,dtype=b.dtype)
        elif pattern == 5:  # diagonal with known sorting order
            a = np.diag(np.linspace(-2.0,2.0,n)).astype(a.dtype)
            b = np.diag(np.linspace(0.8,1.3,n)).astype(b.dtype)
        elif pattern == 6:  # scaled
            a *= 1e90
            b *= 1e-90
        elif pattern == 7:
            a *= 1e-90
            b *= 1e90
        elif pattern == 8:
            a = np.triu(a)*0.07
            b = np.eye(n,dtype=b.dtype) + np.triu(b)*0.05
        elif pattern == 9:  # few isolated infinite eigenvalues
            a = np.triu(a) + np.eye(n,dtype=a.dtype)*2
            b = np.diag(np.r_[np.ones(max(0,n-1)), np.zeros(1)][:n]).astype(b.dtype)
        elif pattern == 10:  # zero A, nonsingular B
            a *= 0.0
        if family in (1,2):
            a = np.real(a)
        if family in (1,3):
            b = np.real(b)
        yield family, supplied, sort, a, b


def selected(alpha, beta, sort):
    if sort == 'none':
        return np.zeros_like(alpha,dtype=bool)
    abs_a, abs_b = np.abs(alpha), np.abs(beta)
    b_floor = 64*np.finfo(float).eps * (max(abs_b,default=0.0))
    nonzero = (abs_a != 0) & (abs_b > b_floor)
    if sort in ('lhp','rhp'):
        prod = np.zeros(len(alpha))
        prod[nonzero] = np.real((alpha[nonzero]/abs_a[nonzero]) *
                                np.conj(beta[nonzero]/abs_b[nonzero]))
        return ((prod < 0) if sort == 'lhp' else (prod > 0)) & nonzero
    if sort == 'iuc':
        return (abs_b > b_floor) & (abs_a < abs_b)
    return ((abs_b <= b_floor) & (abs_a > 0)) | (abs_a > abs_b)


def projective_diff(xa, xb, ya, yb):
    # Chordal/projective distance on homogeneous complex coordinates.
    xnorm = np.hypot(np.abs(xa),np.abs(xb))
    ynorm = np.hypot(np.abs(ya),np.abs(yb))
    x1,x2 = xa/xnorm,xb/xnorm
    y1,y2 = ya/ynorm,yb/ynorm
    distances = abs(x1[:,None]*y2[None,:]-x2[:,None]*y1[None,:])
    row,col = linear_sum_assignment(distances)
    return float(max(distances[row,col],default=0.0))


def to_lines(inputs):
    parts = []
    reference = []
    for family,supplied,sort,a,b in inputs:
        n = len(a)
        if supplied:
            if n:
                aa,bb,q,z = linalg.qz(a,b,output='complex')
            else:
                aa = np.zeros((0,0),complex)
                bb = np.zeros((0,0),complex)
                q = np.zeros((0,0),complex)
                z = np.zeros((0,0),complex)
            x,y = aa,bb
        else:
            x,y = a,b
            q=z=None
        parts.append(f"{n} {family} {supplied} {SORTS.index(sort)+1}\n")
        for j in range(n):
            for i in range(n):
                parts.append(f"{x[i,j].real:.17e} {x[i,j].imag:.17e} "
                             f"{y[i,j].real:.17e} {y[i,j].imag:.17e}\n")
        if supplied:
            for j in range(n):
                for i in range(n):
                    parts.append(f"{q[i,j].real:.17e} {q[i,j].imag:.17e} "
                                 f"{z[i,j].real:.17e} {z[i,j].imag:.17e}\n")
        reference.append((sort,a,b))
    return ''.join(parts),reference


def parse_output(output, reference):
    lines = iter(output.splitlines())
    for sort,a,b in reference:
        n = len(a)
        status,sdim = map(int,next(lines).split())
        data = np.array([[float(v) for v in next(lines).split()] for _ in range(n*n)])
        data = data.reshape(n,n,8).transpose(1,0,2)
        aa = data[:,:,0]+1j*data[:,:,1]
        bb = data[:,:,2]+1j*data[:,:,3]
        q  = data[:,:,4]+1j*data[:,:,5]
        z  = data[:,:,6]+1j*data[:,:,7]
        ab = np.array([[float(v) for v in next(lines).split()] for _ in range(n)]).reshape(n,4)
        alpha = ab[:,0]+1j*ab[:,1]
        beta  = ab[:,2]+1j*ab[:,3]
        yield status, sdim, aa,bb,q,z,alpha,beta


def evaluate(case, actual):
    sort,a,b = case
    status,sdim,aa,bb,q,z,alpha,beta = actual
    n=len(a)
    if status:
        raise AssertionError(f"ordqz returned status {status} for n={n}, sort={sort}")
    if n==0:
        assert sdim==0
        return 0.0,0.0,0.0
    qa = q@aa@z.conj().T
    qb = q@bb@z.conj().T
    rec = max(np.linalg.norm(qa-a)/max(1,np.linalg.norm(a)),
              np.linalg.norm(qb-b)/max(1,np.linalg.norm(b)))
    ident=np.eye(n)
    unit=max(np.linalg.norm(q.conj().T@q-ident)/n,
             np.linalg.norm(z.conj().T@z-ident)/n)
    tri=max(np.linalg.norm(np.tril(aa,-1))/max(1,np.linalg.norm(aa)),
            np.linalg.norm(np.tril(bb,-1))/max(1,np.linalg.norm(bb)))
    assert rec < 2e-7, (rec,n,sort)
    assert unit < 2e-9, (unit,n,sort)
    assert tri < 1e-9, (tri,n,sort)
    assert np.allclose(alpha,np.diag(aa),rtol=1e-13,atol=0)
    assert np.allclose(beta,np.diag(bb),rtol=1e-13,atol=0)
    mask=selected(alpha,beta,sort)
    assert sdim==int(np.sum(mask)), (sdim,mask,sort)
    assert np.all(mask[:sdim]) and not np.any(mask[sdim:]), (mask,sort)
    if sort != 'none':
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',category=RuntimeWarning)
            _,_,ref_a,ref_b,_,_=linalg.ordqz(a,b,sort=sort,output='complex')
    else:
        # Standard generalized eigenvalues are enough for no-sort mode.
        ab_ref=linalg.eigvals(a,b,homogeneous_eigvals=True)
        ref_a,ref_b=ab_ref
    if sort != 'none':
        ref_mask=selected(ref_a,ref_b,sort)
        assert sdim == int(np.count_nonzero(ref_mask)), (sdim,ref_mask,sort)
    spec=projective_diff(alpha,beta,ref_a,ref_b)
    assert spec < 1e-6, (spec,n,sort)
    return rec,unit,spec


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    parser.add_argument('--cases',type=int,default=200)
    parser.add_argument('--max-order',type=int,default=18)
    parser.add_argument('--seed',type=int,default=80515)
    parser.add_argument('--factors-only',action='store_true',
                        help='test the factor-only reorderer beyond the order-40 QZ limit')
    args=parser.parse_args()
    if args.max_order < 2 or args.max_order>(256 if args.factors_only else 40):
        parser.error('max-order must be 2..40, or 2..256 with --factors-only')
    env=os.environ.copy()
    env['OPENBLAS_NUM_THREADS']='1'
    env['OMP_NUM_THREADS']='1'
    env['MKL_NUM_THREADS']='1'
    inputs=list(make_cases(args.cases,args.seed,args.max_order,args.factors_only))
    payload,refs=to_lines(inputs)
    modes=['checked','optimized'] if args.mode=='both' else [args.mode]
    with tempfile.TemporaryDirectory() as temp:
        directory=Path(temp)
        for mode in modes:
            exe=compile_driver(directory,mode)
            proc=subprocess.run([str(exe)],input=payload,text=True,capture_output=True,
                                cwd=directory,env=env,timeout=180)
            if proc.returncode:
                raise RuntimeError(f'{mode} Fortran failed:\n{proc.stderr}')
            results=list(parse_output(proc.stdout,refs))
            maxrec=maxunit=maxspec=0.0
            for idx,(case,actual) in enumerate(zip(refs,results)):
                try:
                    rec,unit,spec=evaluate(case,actual)
                except Exception as exc:
                    raise AssertionError(f'{mode} case #{idx} failed: {exc}') from exc
                maxrec=max(maxrec,rec)
                maxunit=max(maxunit,unit)
                maxspec=max(maxspec,spec)
            if len(results)!=len(refs):
                raise AssertionError('not all Fortran cases returned')
            print(f'{mode}: {len(refs)}/{len(refs)} passed; max reconstruction {maxrec:.3e}; '
                  f'unitarity {maxunit:.3e}; projective eigenvalue difference {maxspec:.3e}',flush=True)


if __name__=='__main__':
    main()
