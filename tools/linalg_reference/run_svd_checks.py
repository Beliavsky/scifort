#!/usr/bin/env python3
"""Compile and compare SciFort economy SVD, pseudoinverse, and least squares with SciPy."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('src/kinds.f90', 'src/linalg.f90', 'src/linalg_svd.f90',
           'src/linalg_qr.f90')


@dataclass
class Case:
    a: np.ndarray
    b: np.ndarray
    rcond: float
    kind: str


def cases_random(count: int, seed: int, max_order: int, stress: bool) -> list[Case]:
    rng = np.random.default_rng(seed)
    sizes = [1, 2, 3, 4, 5, 7, 11, 18, max_order]
    if stress:
        sizes += [30, 45, 64]
    result = []
    for case_id in range(count):
        m = int(rng.choice(sizes))
        n = int(rng.choice(sizes))
        k = min(m,n)
        nrhs = int(rng.choice([0, 1, 2, 3]))
        complex_case = case_id % 2 == 1
        def matrix(shape):
            z = rng.normal(size=shape)
            if complex_case:
                z = z + 1j*rng.normal(size=shape)
            return z
        case_type = ['random', 'exact_rank', 'zero', 'diagonal_spectrum',
                     'repeated_column', 'noisy_rank'][case_id % 6]
        if case_type in ('exact_rank', 'zero', 'diagonal_spectrum', 'noisy_rank'):
            q, _ = linalg.qr(matrix((m,k)),mode='economic')
            v, _ = linalg.qr(matrix((n,k)),mode='economic')
            rank = int(rng.integers(0,k+1)) if case_type != 'zero' else 0
            singular = np.ones(k)
            if case_type == 'diagonal_spectrum':
                singular = np.geomspace(1.,1e-10,k)
                rank=k
            if case_type == 'noisy_rank':
                rank = max(0, k-1)
                singular = np.geomspace(1.,1e-4,k)
            singular[rank:] = 0
            a = (q*singular) @ v.conj().T
            if case_type == 'zero':
                a[:] = 0
            if case_type == 'noisy_rank' and k > 1:
                a += matrix((m,n))*1e-16
        else:
            a = matrix((m,n))/np.sqrt(k)
            if case_type == 'repeated_column' and n > 1:
                a[:,-1] = a[:,0]
        b = matrix((m,nrhs))
        if case_id % 23 == 0:
            a *= 1e-100
            b *= 1e-100
        elif case_id % 29 == 0:
            a *= 1e100
            b *= 1e100
        rcond = max(m,n)*np.finfo(float).eps
        if case_id % 17 == 0:
            rcond = 1e-9
        result.append(Case(a,b,rcond,case_type))
    return result


def compile_exe(work: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    exe = work/f'svd_{mode}'
    command = ['gfortran','-std=f2018',*flags,*(str(ROOT/s) for s in SOURCES),
               str(ROOT/'tools/linalg_reference/driver_svd.f90'),'-o',str(exe)]
    output = subprocess.run(command,cwd=work,capture_output=True,text=True)
    if output.returncode:
        raise RuntimeError(' '.join(command)+'\n'+output.stderr)
    return exe


def run(exe: Path, cases: list[Case], tol: float) -> dict[str, float]:
    lines = [str(len(cases))]
    for case in cases:
        a,b=case.a,case.b
        m,n=a.shape
        lines.append(f'{2 if np.iscomplexobj(a) else 1} {m} {n} {b.shape[1]} {case.rcond:.17e}')
        for matrix in (a,b):
            for z in matrix.ravel(order='F'):
                lines.append(f'{z.real:.17e} {z.imag:.17e}')
    run_result = subprocess.run([str(exe)],input='\n'.join(lines)+'\n', text=True,
                                capture_output=True,timeout=240)
    if run_result.returncode:
        raise RuntimeError(f'Fortran failure {run_result.returncode}:\n{run_result.stderr[-3000:]}')
    output=iter(run_result.stdout.splitlines())
    worst=dict(reconstruction=0.,unitarity=0.,singular_values=0.,pinv=0.,
               lstsq=0.,auto=0.)
    count_rankdef=0
    for index,case in enumerate(cases):
        a,b=case.a,case.b
        m,n=a.shape
        k=min(m,n)
        nrhs=b.shape[1]
        status,pstatus,xstatus,astatus,prank,xrank,arank=map(int,next(output).split())
        if (status,pstatus,xstatus,astatus) != (0,0,0,0):
            raise AssertionError(f'Case {index} {case.kind} {m}x{n}: bad statuses '
                                 f'{status,pstatus,xstatus,astatus}')
        s=np.array([float(next(output)) for _ in range(k)])
        def read_matrix(shape):
            count=np.prod(shape)
            values=np.fromiter((complex(*map(float,next(output).split()))
                                for _ in range(count)),dtype=np.complex128,count=count)
            return values.reshape(shape,order='F')
        u=read_matrix((m,k))
        vh=read_matrix((k,n))
        ap=read_matrix((n,m))
        x=read_matrix((n,nrhs))
        xa=read_matrix((n,nrhs))
        singular=linalg.svdvals(a,check_finite=True)
        if nrhs:
            xref,_,rank_ref,_=linalg.lstsq(a,b,cond=case.rcond,lapack_driver='gelsd')
        else:
            rank_ref=int(np.count_nonzero(singular > singular[0]*case.rcond))
            xref=np.empty((n,0),dtype=a.dtype)
        pref=linalg.pinv(a,atol=0,rtol=case.rcond)
        if not np.iscomplexobj(a):
            if any(np.any(np.imag(q)) for q in (u,vh,ap,x,xa)):
                raise AssertionError(f'real case {index} returned imaginary values')
        if (prank,xrank) != (rank_ref,rank_ref):
            raise AssertionError(f'Case {index} ({m}x{n}, {case.kind}): rank '
                                 f'{prank,xrank} != SciPy {rank_ref}: {s!r} vs {singular!r}')
        if rank_ref < k:
            count_rankdef += 1
        eps=np.finfo(float).tiny
        relative = lambda x,y: float(np.linalg.norm(x-y)/max(eps,np.linalg.norm(y)))
        reconstruction=relative((u*s)@vh,a) if np.linalg.norm(a)>0 else np.linalg.norm((u*s)@vh)
        unitarity=max(np.linalg.norm(u.conj().T@u-np.eye(k))/max(1,k),
                      np.linalg.norm(vh@vh.conj().T-np.eye(k))/max(1,k))
        serror=relative(s,singular) if np.linalg.norm(singular)>0 else np.linalg.norm(s)
        perror=relative(ap,pref) if np.linalg.norm(pref)>0 else np.linalg.norm(ap)
        xerror=relative(x,xref) if np.linalg.norm(xref)>0 else np.linalg.norm(x)
        aerror=relative(xa,xref) if np.linalg.norm(xref)>0 else np.linalg.norm(xa)
        errors=dict(reconstruction=reconstruction,unitarity=unitarity,singular_values=serror,
                    pinv=perror,lstsq=xerror,auto=aerror)
        for key,value in errors.items():
            worst[key]=max(worst[key],value)
        threshold=tol if case.kind not in ('diagonal_spectrum', 'noisy_rank') else max(tol,2e-5)
        if (max(reconstruction,unitarity,serror,perror,xerror) > threshold or aerror > max(threshold,3e-5)):
            raise AssertionError(f'Case {index} m={m} n={n} {case.kind} rank {rank_ref}, '
                                 f'error={errors}, s={s}, scipy={singular}')
    if next(output,None) is not None:
        raise AssertionError('trailing driver output')
    print(f'{len(cases)} SciPy cases passed; {count_rankdef} deficient; '
          + ', '.join(f'{key}={value:.3g}' for key,value in worst.items()))
    return worst


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['both','checked','optimized'],default='both')
    parser.add_argument('--cases',type=int,default=240)
    parser.add_argument('--max-order',type=int,default=32)
    parser.add_argument('--seed',type=int,default=710)
    parser.add_argument('--tolerance',type=float,default=2e-7)
    parser.add_argument('--stress',action='store_true')
    opts=parser.parse_args()
    cases=cases_random(opts.cases,opts.seed,opts.max_order,opts.stress)
    modes=('checked','optimized') if opts.mode=='both' else (opts.mode,)
    with tempfile.TemporaryDirectory(prefix='scifort-svd-') as directory:
        work=Path(directory)
        for mode in modes:
            exe=compile_exe(work,mode)
            print(mode+': ',end='',flush=True)
            run(exe,cases,opts.tolerance)


if __name__ == '__main__':
    main()
