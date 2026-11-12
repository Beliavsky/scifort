#!/usr/bin/env python3
"""Compare SciFort generalized eig/eigvals with SciPy (homogeneous spectrum).

Compiles standalone GNU Fortran in checked and optimized modes.  Compares
projective eigenvalues, normalized left/right pencil residuals, and eigenvalue
representations, not vector phases or column ordering.
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
           "src/linalg_svd.f90", "src/linalg_qr.f90",
           "src/linalg_matrix_functions.f90", "src/linalg_schur.f90",
           "src/linalg_qz.f90", "src/linalg_generalized_eigen.f90"]
ENV = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
       "MKL_NUM_THREADS": "1"}


def generate_cases(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = [0, 1, 2, 3, 4, 6, 9, min(max_order, 16), max_order]
    for k in range(count):
        n = sizes[k % len(sizes)]
        kind = 1 + ((k // len(sizes)) % 4)
        style = (k // (4 * len(sizes))) % 12
        a = rng.normal(size=(n, n))
        b = np.eye(n) + 0.2 * rng.normal(size=(n, n)) / max(np.sqrt(n), 1.)
        if kind in (3, 4):
            a = a + 1j * rng.normal(size=(n, n))
        if kind in (2, 4):
            b = b + 0.2j * rng.normal(size=(n, n)) / max(np.sqrt(n), 1.)
        if style == 1:
            a = np.triu(a)
            b = np.triu(b)
        elif style == 2:
            a = np.tril(a)
            b = np.tril(b)
        elif style == 3:
            a = np.diag(np.linspace(-2, 2, n)).astype(a.dtype)
            b = np.diag(np.linspace(0.6, 2.0, n)).astype(b.dtype)
        elif style == 4:
            a = np.zeros_like(a)
        elif style == 5:
            a = 1e80 * a
            b = 1e-80 * b
        elif style == 6:
            a = 1e-80 * a
            b = 1e80 * b
        elif style == 7:
            a = 1e120 * a
            b = 1e120 * b
        elif style == 8:
            a = 1e-120 * a
            b = 1e-120 * b
        elif style == 9:
            # A regular pencil with some infinite eigenvalues.
            a = np.diag(np.arange(1, n+1, dtype=float)).astype(a.dtype)
            b = np.diag(np.arange(1, n+1, dtype=float)).astype(b.dtype)
            if n:
                b[-1, -1] = 0.
        elif style == 10:
            # Both A and B singular, yet alpha/beta well-defined for each root.
            a = np.diag(np.ones(n)).astype(a.dtype)
            b = np.diag(np.ones(n)).astype(b.dtype)
            if n >= 2:
                a[0, 0] = 0.
                b[-1, -1] = 0.
        elif style == 11:
            # Regular pencil with a repeated finite spectrum.
            a = 2 * b
        if kind in (1, 2):
            a = np.asarray(a.real)
        if kind in (1, 3):
            b = np.asarray(b.real)
        yield kind, np.asarray(a), np.asarray(b)


def compile_exe(folder: Path, mode: str, test: bool):
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace",
              "-ffpe-trap=invalid,zero,overflow"] if mode == "checked"
             else ["-O3"])
    src = ("test/test_linalg_generalized_eigen.f90" if test else
           "tools/linalg_reference/driver_generalized_eigen.f90")
    exe = folder / (mode + ("_native" if test else "_driver"))
    cmd = ["gfortran", "-std=f2018", *flags,
           *(str(ROOT / p) for p in SOURCES), str(ROOT / src), "-o", str(exe)]
    p = subprocess.run(cmd, cwd=folder, text=True, capture_output=True)
    if p.returncode:
        raise RuntimeError(f"Fortran compilation failed: {p.stderr}")
    return exe


def chordal(a, b, c, d):
    """Projective distance after separate scaling of pencil A and B."""
    an = np.hypot(np.abs(a), np.abs(b))
    cn = np.hypot(np.abs(c), np.abs(d))
    return np.abs((a / an)[:, None] * (d / cn)[None, :] -
                  (b / an)[:, None] * (c / cn)[None, :])


def compare(exe: Path, data, tol: float):
    lines = [str(len(data))]
    for kind, a, b in data:
        n = a.shape[0]
        lines.append(f"{kind} {n}")
        lines.extend(f"{x.real:.17e} {x.imag:.17e} {y.real:.17e} {y.imag:.17e}"
                     for x, y in zip(a.ravel(order='F'), b.ravel(order='F')))
    p = subprocess.run([str(exe)], input="\n".join(lines) + "\n", text=True,
                       capture_output=True, timeout=180, env=ENV)
    if p.returncode:
        raise RuntimeError(f"Fortran driver failed: {p.stderr[-3000:]}")
    output = iter(p.stdout.splitlines())
    worst = {"spectrum": 0., "right": 0., "left": 0., "eigvals": 0.}
    for index, (kind, a, b) in enumerate(data):
        n = a.shape[0]
        statuses = [int(x) for x in next(output).split()]
        if statuses != [0, 0]:
            raise AssertionError(f"case {index}: kind={kind} n={n} statuses={statuses}")
        rows = np.array([[float(t) for t in next(output).split()] for _ in range(n)])
        rows = rows.reshape((n, 12))
        vecrows = np.array([[float(t) for t in next(output).split()] for _ in range(n*n)])
        if not n:
            continue
        w = rows[:, 0] + 1j*rows[:, 1]
        alpha = rows[:, 2] + 1j*rows[:, 3]
        beta = rows[:, 4] + 1j*rows[:, 5]
        w2 = rows[:, 6] + 1j*rows[:, 7]
        alpha2 = rows[:, 8] + 1j*rows[:, 9]
        beta2 = rows[:, 10] + 1j*rows[:, 11]
        vr = (vecrows[:,0] + 1j*vecrows[:,1]).reshape(n,n,order='F')
        vl = (vecrows[:,2] + 1j*vecrows[:,3]).reshape(n,n,order='F')
        sa = max(np.max(np.abs(a.real)), np.max(np.abs(a.imag)), 1e-300)
        sb = max(np.max(np.abs(b.real)), np.max(np.abs(b.imag)), 1e-300)
        a0, b0 = a/sa, b/sb
        # SciPy computes generalized homogeneous eigenvalues and both bases.
        ref, _, _ = linalg.eig(a0,b0,left=True,right=True,homogeneous_eigvals=True)
        # The homogeneous coordinates must be scaled separately to make
        # the projective comparison meaningful on extreme A/B scale ratios.
        match = chordal(alpha/sa,beta/sb,ref[0],ref[1])
        rr, cc = linear_sum_assignment(match)
        spec = float(np.max(match[rr, cc]))
        same = float(np.max(chordal(alpha/sa,beta/sb,alpha2/sa,beta2/sb).diagonal()))
        worst["spectrum"] = max(worst["spectrum"], spec)
        worst["eigvals"] = max(worst["eigvals"], same)
        if spec > tol or same > tol:
            raise AssertionError(f"case {index}, kind={kind}, n={n}: "
                                 f"spectrum={spec:.3e}, eigvals={same:.3e}")
        for j in range(n):
            aa = alpha[j]/sa
            bb = beta[j]/sb
            scale = np.hypot(abs(aa),abs(bb))
            if scale == 0.:
                raise AssertionError(f"case {index}: zero homogeneous pair")
            aa, bb = aa/scale, bb/scale
            op = bb*a0-aa*b0
            denom = max(1e-300, abs(bb)*linalg.norm(a0) + abs(aa)*linalg.norm(b0))
            right = float(linalg.norm(op@vr[:,j]))/denom
            left = float(linalg.norm(op.conj().T@vl[:,j]))/denom
            worst["right"] = max(worst["right"], right)
            worst["left"] = max(worst["left"], left)
            if right > tol or left > tol:
                raise AssertionError(f"case {index}, n={n}, root={j}: "
                                     f"right={right:.3e} left={left:.3e}")
            if abs(linalg.norm(vr[:,j])-1.) > 1e-10 or abs(linalg.norm(vl[:,j])-1.) > 1e-10:
                raise AssertionError(f"case {index}: eigenvectors not normalized")
            if abs(beta[j]) == 0. and abs(alpha[j]) > 0.:
                if not np.isinf(w[j].real):
                    raise AssertionError(f"case {index}: infinite eigenvalue not represented")
            elif np.isfinite(w[j].real) and np.isfinite(w[j].imag):
                # Use projective comparison for high- or low-magnitude eigenvalues.
                wf = w[j] * sb/sa
                if np.isfinite(wf) and abs(wf)*abs(beta[j]/sb) < 1e280:
                    if abs((beta[j]/sb)*wf-(alpha[j]/sa)) > 1e-8:
                        raise AssertionError(f"case {index}: w and alpha/beta disagree")
            if not (np.isinf(w2[j]) and np.isinf(w[j])) and abs(w2[j] - w[j]) > 1e-8 * max(1., abs(w[j])):
                raise AssertionError(f"case {index}: eigvals differs from eig")
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['both','checked','optimized'], default='both')
    parser.add_argument('--cases', type=int, default=240)
    parser.add_argument('--max-order', type=int, default=20)
    parser.add_argument('--seed', type=int, default=82026)
    parser.add_argument('--tolerance', type=float, default=1e-7)
    args = parser.parse_args()
    data = list(generate_cases(args.cases,args.max_order,args.seed))
    with tempfile.TemporaryDirectory(prefix='scifort_generalized_eigen_') as temp:
        for mode in (['checked','optimized'] if args.mode == 'both' else [args.mode]):
            folder = Path(temp)
            native = compile_exe(folder,mode,test=True)
            p = subprocess.run([str(native)], cwd=folder, capture_output=True,
                               text=True, env=ENV, timeout=30)
            if p.returncode:
                raise RuntimeError(f"Native test failed: {p.stdout}\n{p.stderr}")
            exe = compile_exe(folder,mode,test=False)
            worst = compare(exe,data,args.tolerance)
            print(f"{mode}: native passed; {len(data)} SciPy comparisons passed; " +
                  ', '.join(f'{k} {v:.3e}' for k,v in worst.items()),flush=True)


if __name__ == '__main__':
    main()
