#!/usr/bin/env python3
"""Independent SciPy comparisons for SciFort's direct QR multiplication.

Usage:
    python reference/check_linalg_qr_multiply.py /path/to/qr_multiply_driver

The Fortran and SciPy QR routines use different legal phase/sign choices.
For nonsingular test matrices, align the corresponding R diagonal phases
before comparing the Q-products. Test the reconstruction separately.
"""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
import numpy as np
from scipy.linalg import qr, qr_multiply


@dataclass
class Case:
    a: np.ndarray
    c: np.ndarray
    side: int
    conj: int
    pivot: int
    vector: int
    da: int
    dc: int


def make_cases() -> list[Case]:
    rng = np.random.default_rng(20261109)
    cases: list[Case] = []
    # Random full-rank input; each pairing and each multiplication convention.
    for m in range(1, 13):
        for n in range(1, 13):
            k = min(m, n)
            for da in (0, 1):
                a = rng.normal(size=(m, n))
                if da:
                    a = a + 1j*rng.normal(size=(m,n))
                for dc in (0, 1):
                    for side in (0, 1):
                        for pivot in (0, 1):
                            conj = (m + 2*n + da + dc + side + pivot) % 2
                            vector = (m + n + da + dc + pivot) % 2
                            p = 1 if vector else 1+(m+n) % 3
                            shape = (k, p) if side == 0 else (p, m)
                            c = rng.normal(size=shape)
                            if dc:
                                c = c + 1j*rng.normal(size=shape)
                            cases.append(Case(a,c,side,conj,pivot,vector,da,dc))
    # Targeted stress: larger rectangular panels, optional conjugation modes,
    # and scales covering many orders of magnitude.
    for t in range(176):
        m = int(rng.integers(1, 25))
        n = int(rng.integers(1, 25))
        k = min(m, n)
        da, dc, side, conj, pivot = [int(x) for x in rng.integers(0,2,size=5)]
        vector = t % 5 == 0
        p = 1 if vector else int(rng.integers(1, 8))
        scale = 10.0**int(rng.integers(-30, 31))
        a = rng.normal(size=(m,n))*scale
        if da:
            a = a + 1j*rng.normal(size=(m,n))*scale
        c = rng.normal(size=(k,p) if side == 0 else (p,m))
        if dc:
            c = c + 1j*rng.normal(size=c.shape)
        cases.append(Case(a,c,side,conj,pivot,int(vector),da,dc))
    return cases


def encode(cases: list[Case]) -> str:
    lines = [str(len(cases))]
    for case in cases:
        a, c = case.a, case.c
        lines.append(f"{case.da} {case.dc} {case.side} {case.conj} "
                     f"{case.pivot} {case.vector} {a.shape[0]} {a.shape[1]} "
                     f"{1 if case.vector else (c.shape[1] if case.side == 0 else c.shape[0])}")
        for arr in (a, c):
            for x in np.asarray(arr).ravel(order='F'):
                lines.append(f"{np.real(x):.17e} {np.imag(x):.17e}")
    return '\n'.join(lines)+'\n'


def compare(exe: str, cases: list[Case]) -> dict[str, float | int]:
    p = subprocess.run([exe], input=encode(cases), text=True, capture_output=True,
                       check=True, timeout=150)
    lines = iter(p.stdout.splitlines())
    max_product = 0.0
    max_r = 0.0
    max_reconstruction = 0.0
    mismatches = []
    for t, case in enumerate(cases):
        a, c = case.a, case.c
        m, n = a.shape
        k = min(m,n)
        status_and_perm = list(map(int, next(lines).split()))
        status, perm = status_and_perm[0], np.array(status_and_perm[1:])-1
        if status != 0:
            raise AssertionError(f"case {t} unexpectedly failed with status {status}")
        if len(perm) != n or sorted(perm) != list(range(n)):
            raise AssertionError(f"case {t}: invalid permutation {perm}")
        pcols = c.shape[1] if case.side == 0 else c.shape[0]
        out_shape = (m,pcols) if case.side == 0 else (pcols,k)
        def fetch(rows,cols):
            return np.array([complex(*map(float,next(lines).split()))
                             for _ in range(rows*cols)]).reshape((rows,cols),order='F')
        rf = fetch(k,n)
        yf = fetch(*out_shape)
        if case.vector:
            c_arg = c[:,0] if case.side == 0 else c[0,:]
        else:
            c_arg = c
        mode = 'left' if case.side == 0 else 'right'
        sc = qr_multiply(a,c_arg,mode=mode,pivoting=bool(case.pivot),conjugate=bool(case.conj))
        if case.pivot:
            ys,rs,ps = sc
            if not np.array_equal(perm,ps):
                raise AssertionError(f"case {t}: permutation differs: Fortran {perm} vs SciPy {ps}")
        else:
            ys,rs = sc
        ys = np.asarray(ys).reshape(out_shape) if case.vector else ys
        qs,rs_ref,*rest = qr(a,mode='economic',pivoting=bool(case.pivot))
        if case.pivot and not np.array_equal(perm,rest[0]):
            raise AssertionError(f"case {t}: QR permutation differs")
        # Right-unitary phases D align the two valid Householder conventions:
        # Q_fortran = Q_scipy D, R_fortran = D^H R_scipy.
        diag_r = np.diag(rf[:k,:k])
        diag_s = np.diag(rs[:k,:k])
        if np.any(np.abs(diag_r) < 1e-200) or np.any(np.abs(diag_s) < 1e-200):
            raise AssertionError(f"case {t}: near-zero diagonal test fixture")
        d = diag_s / diag_r
        d = d / np.abs(d)
        expected_r = d.conj()[:,None] * rs
        # Cross-validate with both scipy.linalg.qr and scipy.linalg.qr_multiply.
        if case.side == 0:
            expected = np.conj(qs) @ (d.conj()[:,None]*c) if case.conj else qs @ (d[:,None]*c)
            if case.vector:
                expected = expected.reshape(out_shape)
        else:
            expected = (c @ np.conj(qs))*d.conj()[None,:] if case.conj else (c @ qs)*d[None,:]
        prod_diff = np.linalg.norm(yf-expected)/(1.0+np.linalg.norm(expected))
        r_diff = np.linalg.norm(rf-expected_r)/(1.0+np.linalg.norm(expected_r))
        perm_a = a[:,perm]
        qf = qs*d[None,:]
        reconstruct = np.linalg.norm(qf@rf - perm_a)/(1.0+np.linalg.norm(perm_a))
        max_product = max(max_product,prod_diff)
        max_r = max(max_r,r_diff)
        max_reconstruction = max(max_reconstruction,reconstruct)
        if prod_diff > 8e-12 or r_diff > 8e-12 or reconstruct > 8e-12:
            mismatches.append((t,(m,n),case.da,case.dc,mode,case.conj,
                               prod_diff,r_diff,reconstruct))
        # The raw SciPy product must itself match the conventional Q product.
        if case.side == 0:
            scipy_direct = np.conj(qs)@c if case.conj else qs@c
        else:
            scipy_direct = c@np.conj(qs) if case.conj else c@qs
        # SciPy 1.17.0 silently discards Im(C) for real A, complex C.
        # Check this combination against SciPy qr() + NumPy matmul instead.
        if not (case.da == 0 and case.dc == 1) and \
                np.linalg.norm(ys-scipy_direct) > 1e-11*(1+np.linalg.norm(scipy_direct)):
            raise AssertionError(f"case {t}: inconsistent SciPy reference")
    if mismatches:
        raise AssertionError(f"{len(mismatches)} mismatches, first 10: {mismatches[:10]}")
    return {'cases':len(cases),'max_product':max_product,'max_r':max_r,
            'max_reconstruction':max_reconstruction}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver',help='Path to compiled Fortran driver')
    args = parser.parse_args()
    cases = make_cases()
    print(compare(args.driver,cases))


if __name__ == '__main__':
    main()
