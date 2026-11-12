#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Validate funm without derivatives on nonnormal diagonalizable repeated spectra.

SciPy's specialized matrix functions and similarity formulas provide two
independent references.  These cases should not require analytic derivatives.
"""
from __future__ import annotations

import argparse
import subprocess
import numpy as np
from scipy import linalg


def make_cases(seed: int, count: int):
    rng = np.random.default_rng(seed)
    records = []
    functions = {
        1: (np.exp, linalg.expm),
        2: (np.log, linalg.logm),
        3: (np.sin, linalg.sinm),
        4: (lambda z: z*z, lambda a: a@a),
        6: (np.cos, linalg.cosm),
    }
    for index in range(count):
        n = int(rng.integers(3, 17)) if index % 7 else int(rng.integers(17, 41))
        is_complex = bool(index % 2)
        fid = sorted(functions)[index % len(functions)]
        # Several groups with interleaved exact repeated eigenvalues.
        nclusters = min(n-1, max(2, n//3))
        eig = np.linspace(0.4, 2.1, nclusters)
        if is_complex:
            eig = eig+1j*rng.uniform(-0.35, 0.35, nclusters)
        eig = eig[np.arange(n) % nclusters]
        if index % 11 == 0:
            eig[::nclusters] += 5e-13  # closely clustered but not identical
        p = np.eye(n, dtype=complex if is_complex else float)
        p += np.triu(rng.normal(size=(n,n)),1)*(0.12 if n <= 16 else 0.025)
        if is_complex:
            p += 1j*np.triu(rng.normal(size=(n,n)),1)*(0.07 if n <= 16 else 0.015)
        pinv = np.linalg.inv(p)
        a = (p*eig) @ pinv
        exact = (p*functions[fid][0](eig)) @ pinv
        similarity = index % 4 == 0
        if similarity:
            q = rng.normal(size=(n,n))
            if is_complex:
                q = q+1j*rng.normal(size=(n,n))
            q, _ = np.linalg.qr(q)
            a = q@a@q.conj().T
            exact = q@exact@q.conj().T
        if not is_complex:
            a = a.real
        # SciPy's generic funm can be ill-conditioned with repeated spectra;
        # use the dedicated algorithms and verify the exact similarity too.
        scipy_result = functions[fid][1](a)
        discrepancy = np.linalg.norm(scipy_result-exact, ord='fro')/max(1.0,np.linalg.norm(exact,ord='fro'))
        if not np.isfinite(discrepancy) or discrepancy > 1e-8:
            raise AssertionError(f'bad reference at {index}: {discrepancy}')
        records.append((a, fid, exact, scipy_result, similarity))
    return records


def compare(executable: str, ncases=250, seed=20261107) -> None:
    records=make_cases(seed,ncases)
    lines=[str(len(records))]
    for a,fid,_,_,_ in records:
        n=a.shape[0]
        lines.append(f'{2 if np.iscomplexobj(a) else 1} {n} {fid} 0')
        for j in range(n):
            for i in range(n):
                z=complex(a[i,j])
                lines.append(f'{z.real:.17e} {z.imag:.17e}')
    proc = subprocess.run([executable], input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
    out=proc.stdout.splitlines()
    pos=0
    failures=[]
    largest=0.0
    for k,(a,fid,exact,scipy_result,similarity) in enumerate(records):
        n=a.shape[0]
        status, residual=out[pos].split()
        pos+=1
        result=np.empty((n,n), dtype=complex)
        for j in range(n):
            for i in range(n):
                re,im=map(float,out[pos].split())
                result[i,j]=re+1j*im
                pos+=1
        err=np.linalg.norm(result-scipy_result,ord='fro')/max(1.0,np.linalg.norm(scipy_result,ord='fro'))
        largest=max(largest,err)
        if int(status)!=0 or not np.isfinite(err) or err>2e-8:
            failures.append(f'#{k} n={n} fid={fid} dense={similarity} status={status} error={err:.3e}')
    print(f'Diagonalizable repeated spectra without derivatives: {len(records)-len(failures)}/{len(records)} passed; '
          f'max relative error {largest:.3e}',flush=True)
    if failures:
        print('\n'.join(failures[:30]),flush=True)
        raise AssertionError(f'{len(failures)} failures')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('--cases',type=int,default=250)
    parser.add_argument('--seed',type=int,default=20261107)
    args=parser.parse_args()
    compare(args.executable,args.cases,args.seed)
