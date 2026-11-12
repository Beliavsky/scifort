#!/usr/bin/env python3
"""Compile dense Householder QR and compare real/complex solves against SciPy."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('src/kinds.f90', 'src/linalg.f90', 'src/linalg_svd.f90', 'src/linalg_qr.f90')


@dataclass
class Case:
    a: np.ndarray
    b: np.ndarray
    pivoted: bool
    complete: bool
    ill_conditioned: bool = False


def build_cases(count: int, seed: int, max_order: int, stress: bool) -> list[Case]:
    rng = np.random.default_rng(seed)
    sizes = [1, 2, 3, 4, 6, 9, 15, max_order]
    if stress:
        sizes += [40, 64, 96, 128, 180]
    cases: list[Case] = []
    for index in range(count):
        m = int(rng.choice(sizes if index % 4 else list(range(1, max_order + 1))))
        n = int(rng.choice(sizes if index % 5 else list(range(1, max_order + 1))))
        nrhs = int(rng.choice([1, 2, 3, 5]))
        complex_case = bool(index % 2)
        a = rng.standard_normal((m, n)) / np.sqrt(max(m, n))
        b = rng.standard_normal((m, nrhs))
        if complex_case:
            a = a + 1j*rng.standard_normal((m, n))/np.sqrt(max(m, n))
            b = b + 1j*rng.standard_normal((m, nrhs))
        # Force the pivoting to move at least one column in many cases.
        if n > 2 and index % 7 == 0:
            a[:, 0] *= 0.01
        ill = index % 29 == 0 and min(m, n) > 1
        if ill:
            u, _ = linalg.qr(rng.standard_normal((m, min(m, n))) +
                             (1j*rng.standard_normal((m, min(m, n))) if complex_case else 0),
                             mode='economic')
            v, _ = linalg.qr(rng.standard_normal((n, min(m, n))) +
                             (1j*rng.standard_normal((n, min(m, n))) if complex_case else 0),
                             mode='economic')
            eig = np.geomspace(1.0, 1e-8, min(m, n))
            a = (u * eig) @ v.conj().T
            b = a @ (rng.standard_normal((n, nrhs)) +
                     (1j*rng.standard_normal((n, nrhs)) if complex_case else 0))
        if index % 23 == 0:
            a = a * 1e-80
            b = b * 1e-80
        elif index % 31 == 0:
            a = a * 1e80
            b = b * 1e80
        cases.append(Case(a, b, pivoted=bool((index // 2) % 2),
                          complete=bool((index // 4) % 2), ill_conditioned=ill))
    return cases


def compile_file(work: Path, source: str, name: str, flags: list[str]) -> Path:
    exe = work / name
    cmd = ['gfortran', '-std=f2018', *flags, *(str(ROOT / s) for s in SOURCES),
           str(ROOT / source), '-o', str(exe)]
    result = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'{" ".join(cmd)}\n{result.stdout}\n{result.stderr}')
    return exe


def run(exe: Path, cases: list[Case], tol: float) -> dict[str, float]:
    lines = [str(len(cases))]
    for case in cases:
        a, b = case.a, case.b
        m, n = a.shape
        lines.append(f'{1 if not np.iscomplexobj(a) else 2} {m} {n} {b.shape[1]} '
                     f'{int(case.pivoted)} {int(case.complete)}')
        for array in (a,b):
            for value in array.ravel(order='F'):
                lines.append(f'{value.real:.17e} {value.imag:.17e}')
    result = subprocess.run([str(exe)], input='\n'.join(lines)+'\n', text=True,
                            capture_output=True, timeout=180)
    if result.returncode:
        raise RuntimeError(f'Fortran driver exited {result.returncode}:\n{result.stderr[-4000:]}')
    output = result.stdout.splitlines()
    pos = 0
    worst = dict(qr=0., orthogonal=0., solution=0., optimality=0.)
    for index, case in enumerate(cases):
        a,b=case.a,case.b
        m,n=a.shape
        nrhs=b.shape[1]
        k=min(m,n)
        qcols=m if case.complete else k
        if pos >= len(output):
            raise AssertionError(f'Missing status for case {index}')
        qr_status,ls_status,rank=map(int,output[pos].split())
        pos+=1
        if qr_status or ls_status or rank != k:
            raise AssertionError(f'Case {index} m={m} n={n}: '
                                 f'qr_status={qr_status}, lstsq_status={ls_status}, rank={rank}')
        piv=np.fromstring(output[pos],dtype=int,sep=' ')-1
        pos+=1
        if piv.size!=n or sorted(piv.tolist())!=list(range(n)):
            raise AssertionError(f'Invalid permutation case {index}')
        if not case.pivoted and not np.array_equal(piv,np.arange(n)):
            raise AssertionError(f'Unpivoted QR permuted columns in case {index}')
        n_q=m*qcols
        n_r=qcols*n
        n_x=n*nrhs
        end=pos+n_q+n_r+n_x
        if end > len(output):
            raise AssertionError(f'Truncated output in case {index}')
        vals=np.fromiter((complex(*map(float, line.split())) for line in output[pos:end]),
                         dtype=np.complex128, count=n_q+n_r+n_x)
        pos=end
        q=vals[:n_q].reshape((m,qcols),order='F')
        r=vals[n_q:n_q+n_r].reshape((qcols,n),order='F')
        x=vals[n_q+n_r:].reshape((n,nrhs),order='F')
        if not np.iscomplexobj(a) and max(np.max(np.abs(q.imag)),np.max(np.abs(r.imag)),
                                        np.max(np.abs(x.imag))) > 0.:
            raise AssertionError(f'Real case {index} returned imaginary components')
        scale=np.linalg.norm(a)
        qrerr=np.linalg.norm(q@r-a[:,piv])/max(np.finfo(float).tiny,scale)
        orth=np.linalg.norm(q.conj().T@q-np.eye(qcols))/max(1.,float(qcols))
        if np.any(np.tril(r,-1)):
            raise AssertionError(f'R is not upper trapezoidal in case {index}')
        x_ref,_,rank_ref,_=linalg.lstsq(a,b,cond=None)
        if rank_ref != k:
            raise AssertionError(f'SciPy reported rank {rank_ref}, unexpected in case {index}')
        xerr=np.linalg.norm(x-x_ref)/max(np.finfo(float).tiny,np.linalg.norm(x_ref))
        res=b-a@x
        if m>=n:
            optimality=np.linalg.norm(a.conj().T@res)/max(np.finfo(float).tiny,
                np.linalg.norm(a)*(np.linalg.norm(a)*np.linalg.norm(x)+np.linalg.norm(b)))
        else:
            optimality=np.linalg.norm(res)/max(np.finfo(float).tiny,
                np.linalg.norm(a)*np.linalg.norm(x)+np.linalg.norm(b))
        worst['qr']=max(worst['qr'],float(qrerr))
        worst['orthogonal']=max(worst['orthogonal'],float(orth))
        worst['solution']=max(worst['solution'],float(xerr))
        worst['optimality']=max(worst['optimality'],float(optimality))
        soltol = 3e-6 if case.ill_conditioned else tol
        if qrerr>5e-11 or orth>5e-11 or xerr>soltol or optimality>soltol:
            raise AssertionError(f'Case {index} m={m} n={n}, pivot={case.pivoted} '
                                 f'complex={np.iscomplexobj(a)}: QR={qrerr:.2e}, '
                                 f'orth={orth:.2e}, x={xerr:.2e}, optimality={optimality:.2e}')
    if pos!=len(output):
        raise AssertionError(f'Unexpected {len(output)-pos} lines of Fortran output')
    print(f'passed {len(cases)}/{len(cases)} SciPy cases; '
          + ', '.join(f'{key}={val:.3e}' for key,val in worst.items()))
    return worst


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('both','checked','optimized'), default='both')
    parser.add_argument('--cases', type=int, default=300)
    parser.add_argument('--max-order', type=int, default=40)
    parser.add_argument('--seed', type=int, default=431)
    parser.add_argument('--tolerance', type=float, default=2e-9)
    parser.add_argument('--stress', action='store_true')
    args=parser.parse_args()
    if args.cases<=0 or args.max_order<=0 or args.tolerance<=0:
        parser.error('cases, max-order and tolerance must be positive')
    cases=build_cases(args.cases,args.seed,args.max_order,args.stress)
    modes=('checked','optimized') if args.mode=='both' else (args.mode,)
    with tempfile.TemporaryDirectory(prefix='scifort-qr-') as dirname:
        work=Path(dirname)
        for mode in modes:
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            native=compile_file(work,'test/test_linalg_qr.f90','native_'+mode,flags)
            output=subprocess.run([str(native)],cwd=work,check=True,capture_output=True,text=True)
            print(f'{mode}: {output.stdout.strip()}')
            driver=compile_file(work,'tools/linalg_reference/driver_qr.f90',
                                'driver_'+mode,flags)
            print(f'{mode}: ',end='',flush=True)
            run(driver,cases,args.tolerance)


if __name__=='__main__':
    main()
