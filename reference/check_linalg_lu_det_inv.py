#!/usr/bin/env python3
"""Reproducible independent scipy.linalg.lu/det/inv comparisons.

Usage: python reference/check_linalg_lu_det_inv.py /path/to/compiled/driver
The Fortran driver is compiled from reference/linalg_lu_det_inv_driver.f90.
"""
from __future__ import annotations
import argparse
import math
import subprocess
import warnings
import numpy as np
from scipy import linalg


def make_cases():
    rng = np.random.default_rng(20261109)
    cases = []
    for m in range(0, 12):
        for n in range(0, 12):
            for variant in range(6):
                for iscomplex in (False, True):
                    a = rng.standard_normal((m, n))
                    if iscomplex:
                        a = a + 1j*rng.standard_normal((m, n))
                    if variant == 1:
                        a *= 0.0
                        for i in range(min(m,n)):
                            a[i,i] = (i+1)/5 * (1+0.75j if iscomplex else 1)
                    elif variant == 2 and n:
                        a[:, 0] = 0.0
                    elif variant == 3 and m > 1:
                        a[-1, :] = a[0, :]
                    elif variant == 4:
                        a *= 1e-12 if (m+n)%2 else 1e12
                    elif variant == 5 and m > 1 and n > 1:
                        a[0, 0] = 0.0
                        a[-1, 0] = 2.0+0j if iscomplex else 2.0
                    cases.append(a)
    for _ in range(160):
        m, n = rng.integers(0, 25, size=2)
        iscomplex = bool(rng.integers(0, 2))
        a = rng.normal(size=(m,n))
        if iscomplex:
            a = a+1j*rng.normal(size=(m,n))
        cases.append(a)
    return cases


def encode(cases):
    stream = [str(len(cases))]
    for a in cases:
        m,n = a.shape
        stream.append(f"{2 if np.iscomplexobj(a) else 1} {m} {n}")
        for x in a.flatten(order='F'):
            if np.iscomplexobj(a):
                stream.append(f"{x.real:.18e} {x.imag:.18e}")
            else:
                stream.append(f"{float(x):.18e}")
    return "\n".join(stream) + "\n"


def unpack(lines, pos, count, dtype):
    if count == 0:
        return np.empty(0,dtype=dtype),pos
    if np.issubdtype(dtype,np.complexfloating):
        out=np.array([complex(*map(float,lines[i].split())) for i in range(pos,pos+count)],dtype=dtype)
    else:
        out=np.array([float(lines[i]) for i in range(pos,pos+count)],dtype=dtype)
    return out, pos+count


def run(driver):
    cases=make_cases()
    cmd=[driver]
    result=subprocess.run(cmd,input=encode(cases),text=True,capture_output=True,check=True)
    lines=result.stdout.splitlines()
    pos=0
    errors=[]
    max_recon = 0.0
    max_det = 0.0
    max_inv = 0.0
    max_factor = 0.0
    compared_factors = 0
    permutations_match = 0
    singular_count = 0
    det_cases = 0
    inv_cases = 0
    near_singular_skipped = 0
    for icase,a in enumerate(cases,1):
        m,n=a.shape
        k=min(m,n)
        dtype=np.complex128 if np.iscomplexobj(a) else np.float64
        line=lines[pos].split();pos+=1
        if line[0]!='CASE' or int(line[1])!=icase:
            raise AssertionError(f'bad record for {icase}: {line}')
        st_lu,st_det,st_inv=map(int,line[2:])
        if st_lu!=0:
            raise AssertionError(f'LU status {st_lu} for case {icase}, shape {a.shape}')
        perm=np.array([int(lines[i]) for i in range(pos,pos+m)],dtype=int);pos+=m
        lv,pos=unpack(lines,pos,m*k,dtype)
        uv,pos=unpack(lines,pos,k*n,dtype)
        l=lv.reshape((m,k),order='F')
        u=uv.reshape((k,n),order='F')
        if m and (sorted(perm.tolist()) != list(range(1,m+1))):
            raise AssertionError(f'invalid permutation in case {icase}')
        if m*n:
            recon=l[perm-1]@u
            err=np.max(np.abs(recon-a)) / max(1e-300, np.max(np.abs(a)))
            max_recon=max(max_recon,float(err))
            if err>2e-11:
                raise AssertionError(f'LU reconstruction error {err} on {icase}, {a.shape}')
        if k:
            if not np.array_equal(np.diag(l[:k,:]),np.ones(k)):
                raise AssertionError('L not unit lower')
            if np.any(np.triu(l,k=1)) or np.any(np.tril(u,k=-1)):
                raise AssertionError('triangular structure violation')
        if m and n:
            p_s,l_s,u_s=linalg.lu(a, p_indices=True)
            if np.array_equal(p_s,perm-1):
                permutations_match+=1
                err_l=np.max(np.abs(l-l_s)) if l.size else 0.0
                err_u=np.max(np.abs(u-u_s))/max(1e-300,np.max(np.abs(u_s))) if u.size else 0.0
                max_factor=max(max_factor,float(max(err_l,err_u)))
                if err_l>2e-11 or err_u>2e-11:
                    raise AssertionError(f'LU factors differ for case {icase} ({err_l} {err_u})')
                compared_factors+=1
        if m==n:
            det_value,pos=unpack(lines,pos,1,dtype)
            inverse_values,pos=unpack(lines,pos,n*n,dtype)
            inverse=inverse_values.reshape((n,n),order='F')
            if st_det!=0:
                raise AssertionError(f'det status {st_det} in case {icase}')
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',linalg.LinAlgWarning)
                reference_det=linalg.det(a)
            # Analytically singular complex matrices may give a nonzero
            # ~epsilon-sized determinant under different pivot algorithms.
            # Scale by the product of the row norms (Hadamard bound), not
            # by a reference determinant that may be rounding noise.
            row_bound = float(np.prod(np.max(np.abs(a),axis=1))) if n else 1.0
            d_err = abs(det_value[0]-reference_det)/max(1e-300,abs(reference_det),row_bound)
            max_det = max(max_det,float(d_err))
            if d_err > 2e-10:
                raise AssertionError(f'det mismatch {d_err} in case {icase}')
            det_cases+=1
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',linalg.LinAlgWarning)
                try:
                    ref_inv=linalg.inv(a)
                    scipy_singular=False
                except linalg.LinAlgError:
                    scipy_singular=True
            # For matrices with effectively zero reciprocal condition number,
            # numerical pivoting can disagree about exact singularity. Such
            # matrices are intentionally excluded from inverse comparisons.
            near_singular = n > 0 and np.linalg.cond(a) > 1.0e10
            if scipy_singular:
                if st_inv == 5:
                    singular_count += 1
                    if inverse.size and np.any(inverse):
                        raise AssertionError(f'inv output not cleared for singular {icase}')
                elif near_singular:
                    near_singular_skipped += 1
                else:
                    raise AssertionError(f'inv failed to report singular {icase}, status={st_inv}')
            elif near_singular:
                near_singular_skipped += 1
            else:
                if st_inv!=0:
                    raise AssertionError(f'inv unexpected failure {icase}, status={st_inv}')
                scale=max(1e-300,float(np.max(np.abs(ref_inv))) if n else 1.0)
                err=float(np.max(np.abs(inverse-ref_inv))/scale) if n else 0.0
                max_inv=max(max_inv,err)
                if err>1e-7:
                    raise AssertionError(f'inv mismatch {err}, case={icase}, shape={a.shape}')
                inv_cases+=1
        else:
            if (st_det,st_inv)!=(-1,-1):
                raise AssertionError(f'non-square det/inv were unexpectedly called for {icase}')
    if pos!=len(lines):
        raise AssertionError(f'unparsed output: pos={pos}, total={len(lines)}')
    print(f'PASS {len(cases)} SciPy LU/det/inv reference cases')
    print(f'LU factor permutations matched in {permutations_match} of '
          f'{sum(x.shape[0]*x.shape[1]>0 for x in cases)} nonempty cases')
    print(f'Individual factors compared in {compared_factors} cases')
    print(f'Determinants: {det_cases}; nonsingular inverses: {inv_cases}; '
          f'expected singular: {singular_count}; near-singular inverses excluded: {near_singular_skipped}')
    print(f'Max relative LU reconstruction: {max_recon:.3e}')
    print(f'Max normalized factor difference: {max_factor:.3e}')
    print(f'Max determinant difference / Hadamard bound: {max_det:.3e}')
    print(f'Max normalized inverse difference: {max_inv:.3e}')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('driver',help='path to compiled Fortran reference driver')
    run(p.parse_args().driver)
