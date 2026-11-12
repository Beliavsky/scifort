#!/usr/bin/env python3
"""Independent SciPy 1.17 LDL validation; compiles a Fortran stdin driver.

SciPy and SciFort can select different 1x1/2x2 pivots. Compare reconstructed
matrices, permitted block structure, and permutation semantics, not raw factors.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import scipy
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("src/kinds.f90", "src/linalg.f90", "src/linalg_ldl.f90",
           "tools/linalg_reference/driver_ldl.f90")


def compile_driver(directory: Path, mode: str) -> Path:
    flags = (["-std=f2018", "-O0", "-g", "-fcheck=all", "-fbacktrace",
              "-ffpe-trap=invalid,zero,overflow"] if mode == "checked"
             else ["-std=f2018", "-O3", "-march=x86-64"])
    binary = directory / f"driver_ldl_{mode}"
    cmd = ["gfortran", *flags, *(str(ROOT / path) for path in SOURCES),
           "-o", str(binary)]
    result = subprocess.run(cmd, cwd=directory, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"compile failed: {result.stderr}")
    return binary


def cases(count: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = (0, 1, 2, 3, 4, 5, 8, 12, 20, 32)
    for case in range(count):
        n = sizes[case % len(sizes)]
        family = 1 + (case // len(sizes)) % 3
        lower = (case // (3*len(sizes))) % 2 == 0
        kind = (case // (6*len(sizes))) % 8
        scale = (1e-80, 1.0, 1e80)[case % 3]
        a = rng.normal(size=(n,n))
        if family != 1:
            a = a + 1j*rng.normal(size=(n,n))
        if family == 2:
            a = a + a.conj().T
        else:
            a = a + a.T
        if kind == 1:
            a += np.eye(n)*n
        if kind == 2:
            np.fill_diagonal(a,0)
        if kind == 3:
            a[:] = 0
            for j in range(0,n-1,2):
                a[j,j+1] = 1.5+0.25j if family != 1 else 1.5
                a[j+1,j] = a[j,j+1].conjugate() if family == 2 else a[j,j+1]
        if kind == 4 and n:
            a[-1,:] = 0
            a[:,-1] = 0
        if kind == 5 and n > 1:
            a[-1,:] = a[0,:]
            a[:,-1] = a[:,0].conj() if family == 2 else a[:,0]
        if kind == 6:
            a = -a
        if kind == 7 and n:
            a += np.eye(n)*0.04
        a *= scale
        yield case, family, lower, a


def input_for_cases(batch):
    lines = []
    for _, family, lower, a in batch:
        n = a.shape[0]
        lines.append(f"{n} {family} {1 if lower else 0}")
        for z in a.T.flat:
            lines.append(f"{z.real:.17e} {z.imag:.17e}")
    return "\n".join(lines) + "\n"


def parse_output(text: str, batch):
    lines = iter(text.splitlines())
    for _, _, _, a in batch:
        n = len(a)
        status = int(next(lines))
        data = np.array([[float(v) for v in next(lines).split()]
                         for _ in range(n*n)], dtype=float).reshape((n,n,4))
        lu = (data[:,:,0] + 1j*data[:,:,1]).T
        d = (data[:,:,2] + 1j*data[:,:,3]).T
        perm = np.array([int(next(lines)) for _ in range(n)]) - 1
        yield status, lu, d, perm


def verify(case, result):
    case_id, family, lower, a = case
    status, lu, d, perm = result
    n = len(a)
    if status != 0:
        raise AssertionError(f"case {case_id}: failed status {status}")
    if n == 0:
        return 0.0
    # The unused input triangle must have no effect on the factorization.
    a = np.tril(a) + np.tril(a,-1).conj().T if lower and family == 2 else (
        np.tril(a) + np.tril(a,-1).T if lower else (
        np.triu(a) + np.triu(a,1).conj().T if family == 2 else
        np.triu(a) + np.triu(a,1).T))
    if sorted(perm.tolist()) != list(range(n)):
        raise AssertionError(f"case {case_id}: invalid permutation")
    tri = lu[perm,:]
    err_tri = np.linalg.norm(tri - (np.tril(tri) if lower else np.triu(tri)), ord=np.inf)
    if err_tri > 2e-12:
        raise AssertionError(f"case {case_id}: nontriangular factor: {err_tri}")
    if not np.allclose(np.diag(tri),1,atol=2e-12,rtol=0):
        raise AssertionError(f"case {case_id}: nonunit factor")
    if np.any(np.abs(d[np.abs(np.subtract.outer(np.arange(n),np.arange(n)))>1]) > 1e-12*max(1,np.max(np.abs(a)))):
        raise AssertionError(f"case {case_id}: nonblock D")
    if np.any((np.abs(np.diag(d,1)[:-1]) > 0) & (np.abs(np.diag(d,1)[1:]) > 0)):
        raise AssertionError(f"case {case_id}: overlapping 2x2 blocks")
    identity = lu @ d @ (lu.conj().T if family == 2 else lu.T)
    scale = max(np.linalg.norm(a,ord=np.inf),1e-300)
    error = np.linalg.norm(identity-a,ord=np.inf)/scale
    if error > 1e-10:
        raise AssertionError(f"case {case_id}: scaled reconstruction error {error}")
    py_lu, py_d, py_perm = linalg.ldl(a,lower=lower,hermitian=(family != 3))
    scipy_recon = py_lu @ py_d @ (py_lu.conj().T if family == 2 else py_lu.T)
    oracle_error = np.linalg.norm(scipy_recon-a,ord=np.inf)/scale
    if oracle_error > 1e-10:
        raise AssertionError(f"case {case_id}: SciPy reference residual {oracle_error}")
    return float(error)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=900)
    parser.add_argument("--seed", type=int, default=20261109)
    parser.add_argument("--mode", choices=["checked", "release"], default="checked")
    args = parser.parse_args()
    batch = list(cases(args.cases,args.seed))
    with tempfile.TemporaryDirectory(prefix="scifort_ldl_") as dirname:
        binary = compile_driver(Path(dirname),args.mode)
        proc = subprocess.run([str(binary)],input=input_for_cases(batch),
                              text=True,capture_output=True)
        if proc.returncode:
            raise RuntimeError(f"Fortran driver failure: {proc.stderr}")
        reference = list(parse_output(proc.stdout,batch))
    errors = [verify(case,result) for case,result in zip(batch,reference,strict=True)]
    print(f"{len(errors)} SciPy {scipy.__version__} LDL reference cases passed; "
          f"maximum relative residual {max(errors,default=0):.3g}; mode={args.mode}")


if __name__ == "__main__":
    main()
