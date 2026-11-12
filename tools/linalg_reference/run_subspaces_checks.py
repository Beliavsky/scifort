#!/usr/bin/env python3
"""Compare SciFort orth, null_space and subspace_angles with SciPy on random matrices.

Compiles a standalone Fortran reference driver; no fpm or external BLAS is needed.
The actual subspace bases are not unique, so validation compares orthogonal
projectors (not raw basis entries), rank, principal angles and residuals.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    'src/kinds.f90', 'src/linalg.f90', 'src/linalg_svd.f90',
    'src/linalg_subspaces.f90',
)


@dataclass
class Case:
    a: np.ndarray
    b: np.ndarray
    kind: int
    rcond: float
    description: str


def generate_cases(count: int, seed: int, max_order: int, stress: bool) -> list[Case]:
    rng = np.random.default_rng(seed)
    sizes = [0, 1, 2, 3, 4, 5, 8, 12, max_order]
    if stress:
        sizes += [16, 24, 32, 48]
    samples = []
    for index in range(count):
        m = int(rng.choice(sizes))
        n = int(rng.choice(sizes))
        p = int(rng.choice(sizes))
        kind = (index % 4) + 1
        complex_a = kind in (2, 4)
        complex_b = kind in (2, 3)

        def random_matrix(shape: tuple[int, int], complex_values: bool) -> np.ndarray:
            z = rng.normal(size=shape)
            if complex_values:
                z = z + 1j * rng.normal(size=shape)
            return z.astype(np.complex128 if complex_values else np.float64)

        rcond = (max(m,n,p,1) * np.finfo(float).eps) if index % 3 else 1e-8
        description = ['full', 'zero', 'exact-rank', 'duplicate',
                       'graded', 'contained', 'near-parallel',
                       'near-perpendicular'][(index // 4 + index % 7) % 8]
        a = random_matrix((m,n),complex_a)
        b = random_matrix((m,p),complex_b)
        if description == 'zero':
            a[...] = 0
        elif description == 'exact-rank' and m*n:
            r = int(rng.integers(0,min(m,n)+1))
            a = random_matrix((m,r),complex_a) @ random_matrix((r,n),complex_a)
        elif description == 'duplicate' and m*n and n > 1:
            a[:,-1] = a[:,0]
        elif description == 'graded' and min(m,n) > 0:
            k = min(m,n)
            u, _ = linalg.qr(random_matrix((m,k),complex_a),mode='economic')
            v, _ = linalg.qr(random_matrix((n,k),complex_a),mode='economic')
            a = (u * np.geomspace(1.,1e-6,k)) @ v.conj().T
        elif description == 'contained' and n and p:
            mix = random_matrix((n,p),complex_b)
            b = a @ mix
            b = b if complex_b else np.real(b).astype(np.float64)
        elif description in ('near-parallel','near-perpendicular') and m >= 2 and n and p:
            a = np.zeros((m,n),dtype=np.complex128 if complex_a else np.float64)
            b = np.zeros((m,p),dtype=np.complex128 if complex_b else np.float64)
            a[0,0] = 1
            delta = 1e-9 if index % 2 else 1e-5
            if description == 'near-parallel':
                b[0,0] = np.cos(delta)
                b[1,0] = np.sin(delta)
            else:
                b[0,0] = np.sin(delta)
                b[1,0] = np.cos(delta)
        if index % 31 == 0:
            a *= 1e-100
            b *= 1e-100
        elif index % 37 == 0:
            a *= 1e100
            b *= 1e100
        samples.append(Case(a,b,kind,float(rcond),description))
    return samples


def compile_driver(directory: Path, mode: str) -> Path:
    exe = directory / f'subspaces_{mode}'
    flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
    command = [
        'gfortran','-std=f2018',*flags,
        *(str(ROOT/s) for s in SOURCES),
        str(ROOT/'tools/linalg_reference/driver_subspaces.f90'),
        '-o',str(exe),
    ]
    proc = subprocess.run(command,cwd=directory,capture_output=True,text=True)
    if proc.returncode:
        raise RuntimeError('Fortran compilation failed:\n' + proc.stderr)
    return exe


def compare(exe: Path, cases: list[Case], tolerance: float) -> dict[str,float]:
    lines = [str(len(cases))]
    for case in cases:
        m,n = case.a.shape
        p = case.b.shape[1]
        lines.append(f'{case.kind} {m} {n} {p} {case.rcond:.18e}')
        for matrix in (case.a,case.b):
            for value in matrix.ravel(order='F'):
                lines.append(f'{value.real:.18e} {value.imag:.18e}')
    proc = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',capture_output=True,
                          text=True,timeout=240)
    if proc.returncode:
        raise RuntimeError(f'Fortran failure {proc.returncode}: {proc.stderr[-3000:]}')
    iterator = iter(proc.stdout.splitlines())
    errors = dict(range_projector=0.,null_projector=0.,
                  null_residual=0.,orthogonality=0.,angles=0.)
    counts = dict(wide=0,rank_deficient=0,mixed=0,empty=0)

    for index,case in enumerate(cases):
        a,b = case.a,case.b
        m,n = a.shape
        so,sz,sa,rank,dim,angle_count = map(int,next(iterator).split())
        if (so,sz,sa) != (0,0,0):
            raise AssertionError(f'Case {index} {case.description} ({m},{n}) '
                                 f'status={so,sz,sa}')

        def read_complex_matrix(shape: tuple[int,int]) -> np.ndarray:
            vals = (complex(*map(float,next(iterator).split()))
                    for _ in range(shape[0]*shape[1]))
            return np.fromiter(vals,dtype=np.complex128,count=shape[0]*shape[1]).reshape(
                shape,order='F')

        q=read_complex_matrix((m,rank))
        z=read_complex_matrix((n,dim))
        angles=np.array([float(next(iterator)) for _ in range(angle_count)])
        qr=linalg.orth(a,rcond=case.rcond)
        zr=linalg.null_space(a,rcond=case.rcond)
        ar=linalg.orth(a,rcond=case.rcond)
        br=linalg.orth(b,rcond=case.rcond)
        theta_ref=linalg.subspace_angles(ar,br)
        if (rank,dim,angle_count) != (qr.shape[1],zr.shape[1],len(theta_ref)):
            raise AssertionError(f'Case {index} {case.description} shape rank '
                                 f'{rank,dim,angle_count} != SciPy '
                                 f'{qr.shape[1],zr.shape[1],len(theta_ref)}')
        if n > m:
            counts['wide'] += 1
        if rank < min(m,n):
            counts['rank_deficient'] += 1
        if case.kind in (3,4):
            counts['mixed'] += 1
        if min(m,n,b.shape[1]) == 0:
            counts['empty'] += 1

        def relative(a1: np.ndarray,a2: np.ndarray) -> float:
            return float(np.linalg.norm(a1-a2)/max(1.,np.linalg.norm(a2)))

        p_q=q@q.conj().T
        p_qr=qr@qr.conj().T
        p_z=z@z.conj().T
        p_zr=zr@zr.conj().T
        ep=relative(p_q,p_qr)
        ez=relative(p_z,p_zr)
        eorth=max(relative(q.conj().T@q,np.eye(rank)),
                  relative(z.conj().T@z,np.eye(dim)))
        enorm=float(np.linalg.norm(a@z)/max(np.linalg.norm(a),np.finfo(float).tiny))
        etheta=float(np.max(np.abs(angles-theta_ref),initial=0))
        vals=dict(range_projector=ep,null_projector=ez,null_residual=enorm,
                  orthogonality=eorth,angles=etheta)
        for key,value in vals.items():
            errors[key]=max(errors[key],value)
        if (max(ep,ez,eorth,etheta) > tolerance or
                enorm > max(tolerance,4*case.rcond) or
                np.any(angles < -1e-14) or np.any(angles > np.pi/2+1e-14) or
                np.any(np.diff(angles)>1e-13)):
            raise AssertionError(f'Case {index} {case.description} dims='
                                 f'{a.shape,b.shape} kind={case.kind} '
                                 f'rcond={case.rcond}: {vals}, '
                                 f'computed angles={angles}, SciPy={theta_ref}')
    if next(iterator,None) is not None:
        raise AssertionError('Unexpected trailing Fortran output')
    print(f'{len(cases)} SciPy cases passed '
          f'({counts["wide"]} wide, {counts["rank_deficient"]} rank-deficient, '
          f'{counts["mixed"]} mixed-type, {counts["empty"]} empty); '
          + ', '.join(f'{key}={value:.3g}' for key,value in errors.items()))
    return errors


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('both','checked','optimized'),default='both')
    parser.add_argument('--cases',type=int,default=240)
    parser.add_argument('--max-order',type=int,default=32)
    parser.add_argument('--seed',type=int,default=1019)
    parser.add_argument('--tolerance',type=float,default=2e-7)
    parser.add_argument('--stress',action='store_true')
    opts=parser.parse_args()
    if opts.cases < 1 or opts.max_order < 1 or opts.tolerance <= 0:
        parser.error('cases, max-order, and tolerance must be positive')
    cases=generate_cases(opts.cases,opts.seed,opts.max_order,opts.stress)
    modes=('checked','optimized') if opts.mode == 'both' else (opts.mode,)
    with tempfile.TemporaryDirectory(prefix='scifort-subspaces-') as tmp:
        work=Path(tmp)
        for mode in modes:
            exe=compile_driver(work,mode)
            print(mode+': ',end='',flush=True)
            compare(exe,cases,opts.tolerance)


if __name__ == '__main__':
    main()
