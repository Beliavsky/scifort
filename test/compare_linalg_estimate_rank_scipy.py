#!/usr/bin/env python3
"""Independent SciPy/NumPy reference comparison for rank estimation.

Usage: python test/compare_linalg_estimate_rank_scipy.py ./driver

SciPy's ndarray estimator typically returns actual rank + 7; the
LinearOperator estimator typically returns actual rank + 1, and either may
saturate at min(m, n). SciFort returns the accepted randomized basis rank.
Exact integer equality to SciPy's randomized estimates is NOT expected.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import scipy
from scipy.linalg import interpolative
from scipy.sparse.linalg import aslinearoperator


def make_case(m, n, k, kind, eps, variant):
    """Build a known-spectrum real or complex matrix using sine bases."""
    r = min(m, n)
    if m == 0 or n == 0:
        a = np.zeros((m, n))
    else:
        u = np.sqrt(2/(m+1)) * np.sin(
            np.pi*np.outer(np.arange(1,m+1),np.arange(1,r+1))/(m+1))
        v = np.sqrt(2/(n+1)) * np.sin(
            np.pi*np.outer(np.arange(1,n+1),np.arange(1,r+1))/(n+1))
        s = np.zeros(r)
        if variant == 0:
            s[:k] = 1.0
        elif variant == 1:
            s[:k] = np.geomspace(1.0, 0.35, max(k,1))[:k]
        elif variant == 2:
            s[:k] = 1.0
            if k < r:
                s[k] = eps*1.e-3
        a = (u * s) @ v.T
    if kind == 1:
        left = np.exp(0.29j*np.arange(1,m+1))
        right = np.exp(-0.37j*np.arange(1,n+1))
        a = (left[:,None]*a)*right[None,:]
    return a


def cases():
    for kind in (0, 1):
        for m, n in ((0,0),(0,5),(7,0),(3,5),(12,9),
                     (24,35),(35,24),(45,55),(55,45),(60,60)):
            for eps in (1.e-6, 1.e-9):
                r = min(m, n)
                ks = sorted(set((0,1,2,4,r//2, max(0,r-12),r)))
                for k in ks:
                    if k > r:
                        continue
                    for variant in (0,1,2):
                        if variant == 2 and (k == 0 or k == r):
                            continue
                        yield kind, m, n, k, eps, variant, make_case(
                            m,n,k,kind,eps,variant)


def prepare(data):
    lines = [str(len(data))]
    for num,(kind,m,n,k,eps,variant,a) in enumerate(data):
        lines.append(f"{kind} {m} {n} {eps:.18e} {987654321+num}")
        for value in a.flat:
            lines.append(f"{value.real:.18e} {value.imag:.18e}")
    return '\n'.join(lines)+'\n'


def main(driver):
    data=list(cases())
    result = subprocess.run([driver], input=prepare(data), text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            check=True, timeout=240)
    lines=result.stdout.splitlines()
    if len(lines) != len(data):
        raise AssertionError((len(data),len(lines),result.stderr[:500]))
    dense_scipy_delta=[]
    operator_scipy_delta=[]
    counts=[0,0]
    ncheck=0
    for i,(test,line) in enumerate(zip(data,lines),1):
        kind,m,n,k,eps,variant,a=test
        fields=line.split()
        if fields[0] != 'CASE' or int(fields[1]) != i:
            raise AssertionError(('malformed driver line',i,line))
        fd,fo,sd,so=map(int,fields[2:])
        if sd != 0 or so != 0:
            raise AssertionError(('Fortran status',i,fields))
        if fd != k or fo != k:
            raise AssertionError(('rank estimate',i,(m,n),kind,k,fd,fo,eps,variant))
        if m and n:
            # SVD provides an independent numerical-rank oracle for the
            # sharply separated spectra used in this reference suite.
            s=np.linalg.svd(a,compute_uv=False)
            true_rank=int(np.count_nonzero(s>eps*s[0])) if s[0] else 0
            if true_rank != k:
                raise AssertionError(('bad test reference',i,k,true_rank))
            # SciPy provides randomized approximations, not a strict equality
            # target. Check that SciPy returns a bounded integer and record
            # its known systematic overestimation for these matrix classes.
            try:
                scipy_dense=int(interpolative.estimate_rank(
                    a,eps,rng=np.random.default_rng(20000+i)))
                scipy_operator=int(interpolative.estimate_rank(
                    aslinearoperator(a),eps,rng=np.random.default_rng(30000+i)))
            except Exception as e:
                raise AssertionError(('SciPy reference failed',i,str(e))) from e
            if not (0<=scipy_dense<=min(m,n)) or not (0<=scipy_operator<=min(m,n)):
                raise AssertionError(('SciPy rank out of bounds',i,scipy_dense,scipy_operator))
            if k>0 and (scipy_dense<k or scipy_operator<k):
                raise AssertionError(('SciPy unexpectedly underestimates',i,k,
                                      scipy_dense,scipy_operator))
            dense_scipy_delta.append(scipy_dense-k)
            operator_scipy_delta.append(scipy_operator-k)
            ncheck+=1
        counts[kind]+=1
    print(f'SciPy {scipy.__version__}: {len(data)} Fortran/reference test cases passed')
    print('  Real/complex test counts:',counts)
    print('  Exact agreement with NumPy known-spectrum SVD rank:',ncheck)
    print('  Dense SciPy minus exact rank range:',min(dense_scipy_delta),max(dense_scipy_delta))
    print('  Operator SciPy minus exact rank range:',min(operator_scipy_delta),max(operator_scipy_delta))
    print('  SciPy randomized outputs are not expected to match exact rank or SciFort.')


if __name__=='__main__':
    if len(sys.argv)!=2:
        raise SystemExit('Usage: compare_linalg_estimate_rank_scipy.py ./driver')
    main(sys.argv[1])
