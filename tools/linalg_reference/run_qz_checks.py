#!/usr/bin/env python3
"""Compile/test SciFort's nonsingular-B generalized complex Schur decomposition.

This is an independent comparison of reconstruction, unitarity, triangularity,
and unordered generalized eigenvalues against scipy.linalg.qz/eigvals. It does
not require identical Schur factors, which are nonunique.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_dense.f90",
           "src/linalg_svd.f90", "src/linalg_qr.f90", "src/linalg_matrix_functions.f90",
           "src/linalg_schur.f90", "src/linalg_qz.f90"]


def cases(count: int, seed: int, max_order: int):
    rng = np.random.default_rng(seed)
    sizes = [0, 1, 2, 3, 4, 6, 9, 12, max_order]
    for k in range(count):
        n = sizes[k % len(sizes)]
        dtype = (k // len(sizes)) % 4 + 1
        complex_a = dtype in (3, 4)
        complex_b = dtype in (2, 4)
        a = rng.standard_normal((n, n))
        b = rng.standard_normal((n, n))
        if complex_a:
            a = a + 1j * rng.standard_normal((n, n))
        if complex_b:
            b = b + 1j * rng.standard_normal((n, n))
        style = (k // (len(sizes)*4)) % 17
        eye = np.eye(n)
        if n:
            b = eye + b / (3 * np.sqrt(n))
        if style == 1:
            a = np.triu(a)
            b = np.triu(b)
            if n:
                b[np.diag_indices(n)] += 0.5
        elif style == 2:
            a = np.tril(a)
            b = np.tril(b)
        elif style == 3:
            a = np.diag(np.diag(a))
            b = np.diag(np.linspace(0.6, 1.5, n)).astype(b.dtype)
        elif style == 4:
            a = np.zeros_like(a)
        elif style == 5:
            b = b * 1e120
            a = a * 1e-120
        elif style == 6:
            b = b * 1e-120
            a = a * 1e120
        elif style == 7:
            b = b * 1e220
            a = a * 1e220
        elif style == 8:
            b = b * 1e-220
            a = a * 1e-220
        elif style == 9:
            if n:
                b = b @ np.diag(np.geomspace(1, 1e-5, n))
        elif style == 10:
            if n:
                a = b @ np.triu(a)
        elif style == 11:
            if n:
                a = 1e-4 * eye + 0.1 * a
        elif style == 12:
            if n:
                a = a + 4.0 * eye
                b = np.diag(np.r_[np.ones(n-1), 0.0])
        elif style == 13:
            a = np.diag(np.linspace(1.0, 2.0, n)).astype(a.dtype)
            b = np.zeros_like(b)
        elif style == 14:
            if n:
                a = a + 5.0 * eye
                d = np.ones(n)
                d[::3] = 0.0
                b = np.diag(d).astype(b.dtype)
        elif style == 15 and n > 1:
            da = np.linspace(1.0, 2.0, n)
            db = np.linspace(1.5, 2.5, n)
            da[0] = 0.0
            db[-1] = 0.0
            a = np.diag(da).astype(a.dtype)
            b = np.diag(db).astype(b.dtype)
        elif style == 16 and n > 1:
            da = np.linspace(1.0, 2.0, n)
            db = np.linspace(1.5, 2.5, n)
            da[0] = 0.0
            db[-1] = 0.0
            a = np.triu(a * 0.03) + np.diag(da - np.diag(a)*0.03)
            b = np.triu(b * 0.03) + np.diag(db - np.diag(b)*0.03)
        # Respect the explicitly advertised real/complex argument types.
        if not complex_a:
            a = np.real(a)
        if not complex_b:
            b = np.real(b)
        yield dtype, np.asarray(a), np.asarray(b)


def compile_exe(directory: Path, mode: str) -> Path:
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = directory / f"qz_{mode}"
    cmd = ["gfortran", "-std=f2018", *flags, *(str(ROOT / s) for s in SOURCES),
           str(ROOT / "tools/linalg_reference/driver_qz.f90"), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=directory, text=True, capture_output=True)
    if p.returncode:
        raise RuntimeError(f"Compilation failed:\n{p.stderr}")
    return exe


def run_native(directory: Path, mode: str):
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
             if mode == "checked" else ["-O3"])
    exe = directory / f"qz_native_{mode}"
    cmd = ["gfortran", "-std=f2018", *flags, *(str(ROOT / s) for s in SOURCES),
           str(ROOT / "test/test_linalg_qz.f90"), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=directory, text=True, capture_output=True)
    if p.returncode:
        raise RuntimeError(f"Native compilation failed:\n{p.stderr}")
    p = subprocess.run([str(exe)], cwd=directory, text=True, capture_output=True)
    if p.returncode:
        raise RuntimeError(f"Native tests failed:\n{p.stdout}\n{p.stderr}")
    print(f"{mode}: {p.stdout.strip()}")


def run_compare(exe: Path, data, method: str = "auto"):
    entries = [str(len(data))]
    for kind, a, b in data:
        n = a.shape[0]
        entries.append(f"{kind} {n}")
        entries.extend(f"{x.real:.17e} {x.imag:.17e} {y.real:.17e} {y.imag:.17e}"
                       for x, y in zip(a.ravel(order='F'), b.ravel(order='F')))
    payload = '\n'.join(entries) + '\n'
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
    p = subprocess.run([str(exe), method], input=payload, capture_output=True, text=True,
                       timeout=180, env=env)
    if p.returncode:
        raise RuntimeError(f"Driver failed:\n{p.stderr[-2000:]}")
    lines = iter(p.stdout.splitlines())
    worst = {"reconstruct": 0.0, "unitary": 0.0, "eigenvalue": 0.0,
             "ill_conditioned_spectra": 0, "backward_eig": 0.0}
    for k, (kind, a, b) in enumerate(data):
        n = a.shape[0]
        status = int(next(lines))
        if status:
            raise AssertionError(f"case {k}: kind={kind}, n={n}, status={status}")
        vals = np.array([[float(s) for s in next(lines).split()] for _ in range(n*n)])
        if n:
            def from_columns(i):
                return (vals[:, i]+1j*vals[:, i+1]).reshape((n,n),order='F')
            aa, bb, q, z = (from_columns(i) for i in (0,2,4,6))
            ab = np.array([[float(s) for s in next(lines).split()] for _ in range(n)])
            alpha = ab[:,0] + 1j*ab[:,1]
            beta = ab[:,2] + 1j*ab[:,3]
        else:
            continue
        scale_a = max(float(np.max(np.abs(a.real))),float(np.max(np.abs(a.imag))),1e-300)
        scale_b = max(float(np.max(np.abs(b.real))),float(np.max(np.abs(b.imag))),1e-300)
        norm_a = max(float(np.linalg.norm(a / scale_a)),1.)
        norm_b = max(float(np.linalg.norm(b / scale_b)),1.)
        da = float(np.linalg.norm((q@(aa/scale_a)@z.conj().T)-(a/scale_a))) / norm_a
        db = float(np.linalg.norm((q@(bb/scale_b)@z.conj().T)-(b/scale_b))) / norm_b
        tri = max(np.linalg.norm(np.tril(aa/scale_a,-1)),
                  np.linalg.norm(np.tril(bb/scale_b,-1)))
        unitarity = max(np.linalg.norm(q.conj().T@q-np.eye(n)),
                        np.linalg.norm(z.conj().T@z-np.eye(n)))
        beta_phase = np.max(np.abs(np.imag(beta/scale_b)))
        if np.min(np.real(beta/scale_b)) < -1e-12 or beta_phase > 1e-12:
            raise AssertionError(f"case {k}: beta not real nonnegative")
        if tri > 1e-12 or max(da,db) > 5e-8 or unitarity > 1e-8:
            raise AssertionError(f"case {k}: bad factors: da={da}, db={db}, "
                                 f"unitary={unitarity}, triangular={tri}")
        sa, sb = a/scale_a, b/scale_b
        ref_aa, ref_bb, _, _ = linalg.qz(sa,sb,output='complex')
        our_alpha, our_beta = alpha/scale_a, beta/scale_b
        ref_alpha, ref_beta = np.diag(ref_aa), np.diag(ref_bb)
        ours_norm = np.hypot(np.abs(our_alpha),np.abs(our_beta))
        ref_norm = np.hypot(np.abs(ref_alpha),np.abs(ref_beta))
        if min(np.min(ours_norm),np.min(ref_norm)) <= 0:
            raise AssertionError(f"case {k}: indeterminate generalized eigenvalue")
        # Homogeneous/chordal comparison works also for beta=0 (infinite roots).
        chordal = np.abs(our_alpha[:,None]*ref_beta[None,:] -
                         our_beta[:,None]*ref_alpha[None,:]) / (
                         ours_norm[:,None]*ref_norm[None,:])
        rr, cc = linear_sum_assignment(chordal)
        projective = float(np.max(chordal[rr,cc]))
        if projective > 5e-8:
            # Forward homogeneous eigenvalues are not stable for nearly
            # defective pencils. Verify that the approximate pencil has a
            # singular value near zero at each computed homogeneous root.
            backward = 0.0
            for al, be in zip(our_alpha, our_beta):
                sigma = float(linalg.svdvals(be*sa-al*sb)[-1])
                backward = max(backward, sigma / max(1e-300,
                    abs(be)*linalg.norm(sa)+abs(al)*linalg.norm(sb)))
            if backward > 2e-8:
                raise AssertionError(f"case {k}: projective error={projective}, "
                                     f"pencil backward error={backward}")
            worst['ill_conditioned_spectra'] += 1
            worst['backward_eig'] = max(worst['backward_eig'],float(backward))
        eigen_error = projective
        if not np.any(np.abs(our_beta) < 1e-14) and not np.any(np.abs(ref_beta) < 1e-14):
            ref = linalg.eigvals(sa,sb)
            ours = our_alpha/our_beta
            costs = np.abs(ours[:,None]-ref[None,:])
            rr, cc = linear_sum_assignment(costs)
            eigen_forward = np.max(costs[rr,cc]) / max(1.,np.max(np.abs(ref)))
            if eigen_forward > 5e-8:
                # In nearly defective nonnormal matrices, two backward-stable
                # algorithms can disagree in forward eigenvalues. Test the
                # smallest singular value of every computed pencil instead.
                backward = 0.0
                for lam in ours:
                    sigma = float(linalg.svdvals(sa-lam*sb)[-1])
                    backward = max(backward, sigma / max(1.,linalg.norm(sa)+abs(lam)*linalg.norm(sb)))
                if backward > 2e-8:
                    raise AssertionError(f"case {k}: forward error={eigen_forward}, "
                                         f"eigenvalue backward error={backward}")
                worst['ill_conditioned_spectra'] += 1
                worst['backward_eig'] = max(worst['backward_eig'],float(backward))
        if np.max(np.abs(alpha-np.diag(aa))/scale_a) > 1e-9 or \
            np.max(np.abs(beta-np.diag(bb))/scale_b) > 1e-9:
            raise AssertionError(f"case {k}: optional diagonal outputs wrong")
        worst['reconstruct'] = max(worst['reconstruct'],da,db)
        worst['unitary'] = max(worst['unitary'],float(unitarity))
        worst['eigenvalue'] = max(worst['eigenvalue'],float(eigen_error))
    return worst


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['both','checked','optimized'],default='both')
    parser.add_argument('--cases',type=int,default=180)
    parser.add_argument('--method',choices=['auto','iterative','inverse'],default='auto')
    parser.add_argument('--max-order',type=int,default=32)
    parser.add_argument('--seed',type=int,default=73131)
    args = parser.parse_args()
    data = list(cases(args.cases,args.seed,args.max_order))
    with tempfile.TemporaryDirectory(prefix='scifort_qz_') as temp:
        work = Path(temp)
        for mode in (['checked','optimized'] if args.mode == 'both' else [args.mode]):
            exe = compile_exe(work,mode)
            run_native(work,mode)
            worst = run_compare(exe,data,args.method)
            print(f"{mode}: passed {len(data)} SciPy comparisons, worst={worst}")


if __name__ == '__main__':
    main()
