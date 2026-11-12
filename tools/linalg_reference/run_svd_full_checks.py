#!/usr/bin/env python3
"""Verify SciFort's full real/complex SVD and svdvals against SciPy.

Compile stand-alone GNU Fortran programs, stream test matrices, and compare
reconstruction, both square-unitary factors, spectra, and unused null spaces.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('src/kinds.f90', 'src/linalg.f90', 'src/linalg_svd.f90')


def make_cases(count: int, seed: int, max_order: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    sizes = sorted(set([1, 2, 3, 4, 5, 7, 11, 17, max_order]))
    cases = []
    for i in range(count):
        m, n = (int(rng.choice(sizes)), int(rng.choice(sizes)))
        k = min(m,n)
        complex_case = (i % 2 == 1)
        def noise(shape):
            a = rng.standard_normal(shape)
            if complex_case:
                a = a + 1j*rng.standard_normal(shape)
            return a
        family = i % 7
        if family in (2, 3, 4, 5):
            q, _ = linalg.qr(noise((m,k)),mode='economic')
            v, _ = linalg.qr(noise((n,k)),mode='economic')
            if family == 2:
                singular = np.r_[np.ones(max(1,k//2)), np.zeros(k-max(1,k//2))]
            elif family == 3:
                singular = np.zeros(k)
            elif family == 4:
                singular = np.geomspace(1,1e-9,k)
            else:
                singular = np.ones(k)
            a = (q*singular)@v.conj().T
            if family == 3:
                a[:] = 0
        else:
            a = noise((m,n))/np.sqrt(k)
            if family == 1 and n > 1:
                a[:,-1] = a[:,0]
            if family == 6 and m > 1:
                a[-1,:] = a[0,:]
        if i % 31 == 0:
            a = a*1e-110
        elif i % 37 == 0:
            a = a*1e110
        cases.append(np.asarray(a))
    return cases


def compile_driver(tmp: Path, mode: str) -> Path:
    flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
    exe = tmp/('full_svd_'+mode)
    cmd = ['gfortran','-std=f2018',*flags,*(str(ROOT/p) for p in SOURCES),
           str(ROOT/'tools/linalg_reference/driver_svd_full.f90'),'-o',str(exe)]
    proc = subprocess.run(cmd,cwd=tmp,check=False,capture_output=True,text=True)
    if proc.returncode:
        raise RuntimeError('compile failed:\n'+proc.stderr)
    return exe


def run(exe: Path, cases: list[np.ndarray], tol: float) -> dict[str,float]:
    lines = [str(len(cases))]
    for a in cases:
        m,n = a.shape
        lines.append(f'{2 if np.iscomplexobj(a) else 1} {m} {n}')
        for z in a.ravel(order='F'):
            lines.append(f'{z.real:.17e} {z.imag:.17e}')
    proc = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,
                          capture_output=True,timeout=240)
    if proc.returncode:
        raise RuntimeError('Fortran failure:\n'+proc.stderr[-3000:])
    output = iter(proc.stdout.splitlines())
    worst = dict(reconstruction=0.,unitarity=0.,svdvals=0.,spectrum=0.,
                 diagsvd=0.,extra_kernel=0.,full_projector=0.)
    for index,a in enumerate(cases):
        m,n = a.shape
        k = min(m,n)
        st, stvals, stdiag = map(int,next(output).split())
        if (st,stvals,stdiag) != (0,0,0):
            raise AssertionError(f'case {index}: statuses {st,stvals,stdiag}, shape {a.shape}')
        s = np.array([float(next(output)) for _ in range(k)])
        sval = np.array([float(next(output)) for _ in range(k)])
        def read_matrix(shape):
            count = int(np.prod(shape))
            values = np.fromiter((complex(*map(float,next(output).split()))
                                  for _ in range(count)),dtype=np.complex128,count=count)
            return values.reshape(shape,order='F')
        u = read_matrix((m,m))
        vh = read_matrix((n,n))
        d = np.array([float(next(output)) for _ in range(m*n)]).reshape((m,n),order='F')
        if not np.iscomplexobj(a) and (np.any(u.imag) or np.any(vh.imag)):
            raise AssertionError(f'case {index}: imaginary results for real data')
        u_ref, s_ref, vh_ref = linalg.svd(a,full_matrices=True)
        anorm = np.linalg.norm(a)
        safe = max(np.finfo(float).tiny, anorm)
        recons = np.linalg.norm((u[:,:k]*s)@vh[:k,:]-a)/safe
        decomp = np.linalg.norm(u@d@vh-a)/safe
        diag_err = np.linalg.norm(d-linalg.diagsvd(s,m,n))/max(1.,np.linalg.norm(s))
        uerr = np.linalg.norm(u.conj().T@u-np.eye(m))/max(1,m)
        verr = np.linalg.norm(vh@vh.conj().T-np.eye(n))/max(1,n)
        serr = np.linalg.norm(s-s_ref)/max(np.finfo(float).tiny,np.linalg.norm(s_ref))
        sverr = np.linalg.norm(sval-s_ref)/max(np.finfo(float).tiny,np.linalg.norm(s_ref))
        ker = 0.
        projector = 0.
        if m > k:
            ker = max(ker,np.linalg.norm(a.conj().T@u[:,k:])/safe)
            if s_ref[-1] > 1e-8*s_ref[0]:
                projector = np.linalg.norm(u[:,k:]@u[:,k:].conj().T -
                                           u_ref[:,k:]@u_ref[:,k:].conj().T)
        if n > k:
            ker = max(ker,np.linalg.norm(a@vh[k:,:].conj().T)/safe)
            if s_ref[-1] > 1e-8*s_ref[0]:
                projector = np.linalg.norm(vh[k:,:].conj().T@vh[k:,:] -
                                           vh_ref[k:,:].conj().T@vh_ref[k:,:])
        errors = dict(reconstruction=max(recons,decomp),unitarity=max(uerr,verr),
                      svdvals=sverr,spectrum=serr,diagsvd=diag_err,extra_kernel=ker,
                      full_projector=projector)
        for key,val in errors.items():
            worst[key] = max(worst[key],float(val))
        if any(val > tol for val in errors.values()):
            raise AssertionError(f'case {index} {a.shape}: '+str(errors))
    if next(output,None) is not None:
        raise AssertionError('Unexpected trailing driver output')
    print(f'{len(cases)} full-SVD SciPy comparisons passed; '+
          ', '.join(f'{k}={v:.3g}' for k,v in worst.items()))
    return worst


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',type=int,default=250)
    parser.add_argument('--seed',type=int,default=805)
    parser.add_argument('--max-order',type=int,default=30)
    parser.add_argument('--mode',choices=['both','checked','optimized'],default='both')
    parser.add_argument('--tol',type=float,default=2e-8)
    args = parser.parse_args()
    cases = make_cases(args.cases,args.seed,args.max_order)
    with tempfile.TemporaryDirectory(prefix='scifort_full_svd_') as work:
        tmp = Path(work)
        for mode in ('checked','optimized') if args.mode == 'both' else (args.mode,):
            print(mode+':',end=' ',flush=True)
            run(compile_driver(tmp,mode),cases,args.tol)


if __name__ == '__main__':
    main()
