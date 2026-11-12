#!/usr/bin/env python3
"""Compare SciFort LDL-based indefinite solves with scipy.linalg.solve.

Tests solve_indefinite and reusable ldl + ldl_solve over real symmetric,
complex Hermitian, and complex symmetric matrices in both orientations.
The supplied matrix has deliberately contradictory unused-triangle entries.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
import scipy
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    "src/kinds.f90", "src/linalg.f90", "src/linalg_ldl.f90",
    "src/linalg_indefinite_solve.f90", "tools/linalg_reference/driver_indefinite_solve.f90",
)


@dataclass
class Case:
    case_id: int
    family: int
    lower: bool
    a: np.ndarray
    b: np.ndarray
    expected_x: np.ndarray
    full_a: np.ndarray


def compile_driver(directory: Path, mode: str) -> Path:
    flags = (["-std=f2018", "-O0", "-g", "-fcheck=all", "-fbacktrace",
              "-ffpe-trap=invalid,zero,overflow"] if mode == "checked"
             else ["-std=f2018", "-O3", "-march=x86-64"])
    binary = directory / f"driver_indefinite_solve_{mode}"
    cmd = ["gfortran", *flags, *(str(ROOT / p) for p in SOURCES), "-o", str(binary)]
    result = subprocess.run(cmd, cwd=directory, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"Compilation failed:\n{result.stderr}")
    return binary


def cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = (0, 1, 2, 3, 4, 5, 8, 12, 20, 32)
    for case_id in range(count):
        n = sizes[case_id % len(sizes)]
        family = 1 + (case_id // len(sizes)) % 3
        lower = bool((case_id // (3*len(sizes))) % 2)
        kind = (case_id // (6*len(sizes))) % 6
        nrhs = (0, 1, 2, 5)[(case_id // len(sizes)) % 4]
        scale = (1e-80, 1.0, 1e80)[case_id % 3]
        if family == 1:
            raw = rng.normal(size=(n,n))
        else:
            raw = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
        if family == 2:
            a = raw + raw.conj().T
        else:
            a = raw + raw.T
        if kind == 1:
            # Guaranteed nonsingular if all random noise were omitted.
            a = (rng.normal(size=(n,n)) +
                 (1j*rng.normal(size=(n,n)) if family != 1 else 0)) * 0.005
            a = a+a.conj().T if family == 2 else a+a.T
            np.fill_diagonal(a, 0)
            for j in range(0,n-1,2):
                v = 1.5+0.6j if family != 1 else 1.5
                a[j,j+1] = v
                a[j+1,j] = v.conjugate() if family == 2 else v
            if n % 2:
                a[-1,-1] = 1.0
        elif kind == 2 and n > 1:
            # A late pivot and a near-zero diagonal promote 2x2 blocks.
            a[-1,-1] = 0.0
            a[0,0] = 0.0
            a[0,-1] = 3.0+0.7j if family != 1 else 3.0
            a[-1,0] = a[0,-1].conjugate() if family == 2 else a[0,-1]
        elif kind == 3:
            a = np.diag([((-1)**j)*(n+2) for j in range(n)]).astype(a.dtype)
            a += (raw+raw.conj().T if family == 2 else raw+raw.T)*0.03
        elif kind == 4 and n > 0:
            a += np.eye(n)*((-1)**(case_id % 2))*0.3
        elif kind == 5 and n > 0:
            a += np.diag([((-1)**j)*(n+2) for j in range(n)])
        if n == 1 and abs(a[0,0]) == 0:
            a[0,0] = 1.0
        a *= scale
        true_x = rng.normal(size=(n,nrhs))
        if family != 1:
            true_x = true_x + 1j*rng.normal(size=(n,nrhs))
        b = a @ true_x
        supplied = a.copy()
        for j in range(n):
            for i in range(n):
                if (lower and i < j) or (not lower and i > j):
                    supplied[i,j] = (-12.25 + (7.5j if family != 1 else 0))*scale
        yield Case(case_id,family,lower,supplied,b,true_x,a)


def fortran_input(cases_to_run):
    lines = []
    for case in cases_to_run:
        n, nrhs = case.b.shape
        lines.append(f"{n} {nrhs} {case.family} {int(case.lower)}")
        for z in case.a.T.flat:
            lines.append(f"{z.real:.17e} {z.imag:.17e}")
        for z in case.b.T.flat:
            lines.append(f"{z.real:.17e} {z.imag:.17e}")
    return "\n".join(lines)+"\n"


def parse_output(out, cases_to_run):
    lines = iter(out.splitlines())
    for case in cases_to_run:
        n, nrhs = case.b.shape
        sd, sf = map(int, next(lines).split())
        data = np.array([[float(v) for v in next(lines).split()]
                         for _ in range(n*nrhs)], dtype=float)
        data = data.reshape((nrhs,n,4)).transpose((1,0,2))
        xd = data[:,:,0] + 1j*data[:,:,1]
        xf = data[:,:,2] + 1j*data[:,:,3]
        yield sd, sf, xd, xf


def verify(case: Case, result):
    sd, sf, xd, xf = result
    n, nrhs = case.b.shape
    if sd != 0 or sf != 0:
        raise AssertionError(f"case {case.case_id}, family={case.family}, "
                             f"lower={case.lower}, size={n}: statuses {sd} {sf}")
    if n == 0 or nrhs == 0:
        return 0.0
    symmetric = case.family != 2
    actual_a = (np.tril(case.a) +
                (np.tril(case.a,-1).T if symmetric else np.tril(case.a,-1).conj().T)
                if case.lower else np.triu(case.a) +
                (np.triu(case.a,1).T if symmetric else np.triu(case.a,1).conj().T))
    scale = max(1e-300,np.linalg.norm(actual_a,ord=np.inf))
    if not np.allclose(actual_a,case.full_a,rtol=0,atol=1e-14*scale):
        raise AssertionError(f"case {case.case_id}: incorrectly synthesized input")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", linalg.LinAlgWarning)
        oracle = linalg.solve(case.a,case.b,lower=case.lower,
                              assume_a="her" if case.family == 2 else "sym")
    # Separate numpy oracle validates the full reconstructed matrix.
    numpy_x = np.linalg.solve(actual_a,case.b)
    norm_true = max(1e-300,np.linalg.norm(numpy_x))
    error = max(np.linalg.norm(xd-numpy_x),np.linalg.norm(xf-numpy_x))/norm_true
    discrepancy = np.linalg.norm(xd-xf)/norm_true
    oracle_error = np.linalg.norm(oracle-numpy_x)/norm_true
    residual = max(np.linalg.norm(actual_a@xd-case.b),
                   np.linalg.norm(actual_a@xf-case.b)) / max(1e-300,np.linalg.norm(case.b))
    if not (np.all(np.isfinite(xd)) and np.all(np.isfinite(xf))):
        raise AssertionError(f"case {case.case_id}: nonfinite outputs")
    if error > 3e-9 or residual > 2e-9 or discrepancy > 2e-12 or oracle_error > 1e-8:
        raise AssertionError(f"case {case.case_id} family={case.family} n={n}, "
                             f"errors: solution={error:.3g} residual={residual:.3g} "
                             f"factored={discrepancy:.3g} SciPy={oracle_error:.3g}")
    return float(max(error,residual))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cases",type=int,default=900)
    p.add_argument("--seed",type=int,default=20261109)
    p.add_argument("--mode",choices=["checked","release"],default="checked")
    args = p.parse_args()
    batch = list(cases(args.cases,args.seed))
    with tempfile.TemporaryDirectory(prefix="scifort_ldl_solve_") as dirname:
        binary = compile_driver(Path(dirname),args.mode)
        result = subprocess.run([str(binary)],input=fortran_input(batch),
                                text=True,capture_output=True)
        if result.returncode:
            raise RuntimeError(f"Fortran execution failed: {result.stderr}")
        actual = list(parse_output(result.stdout,batch))
    errors = [verify(case,res) for case,res in zip(batch,actual,strict=True)]
    print(f"{len(errors)} SciPy {scipy.__version__} indefinite solve comparisons "
          f"passed; maximum relative error {max(errors,default=0):.3g}; "
          f"mode={args.mode}")


if __name__ == "__main__":
    main()
