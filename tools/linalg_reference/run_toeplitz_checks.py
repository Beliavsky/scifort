#!/usr/bin/env python3
"""Independent NumPy/SciPy checks of SciFort's Toeplitz solvers.

Run from any working directory. Requires gfortran, NumPy, and SciPy.
Builds only the relevant SciFort modules, not the whole FPM project.
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import solve_toeplitz as scipy_solve_toeplitz
from scipy.linalg import toeplitz

ROOT = Path(__file__).resolve().parents[2]


def compile_driver(build, mode):
    flags = ["-std=f2018", "-ffree-line-length-none"]
    if mode == "checked":
        flags += ["-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
    else:
        flags += ["-O3", "-march=native"]
    cmd = ["gfortran", *flags, "-J", str(build), "-I", str(build),
           "-o", str(build / "driver_toeplitz")]
    cmd += [str(ROOT / "src" / name) for name in
            ("kinds.f90", "linalg.f90", "linalg_dense.f90", "linalg_toeplitz.f90")]
    cmd += [str(ROOT / "tools/linalg_reference/driver_toeplitz.f90")]
    subprocess.run(cmd, check=True, cwd=build, capture_output=True, text=True)
    return build / "driver_toeplitz"


def numbers(v, is_complex):
    if is_complex:
        return " ".join(f"({z.real:.17e},{z.imag:.17e})" for z in
                        np.asarray(v).ravel(order="F"))
    return " ".join(f"{z:.17e}" for z in np.asarray(v).ravel(order="F"))


def make_case(rng, idx):
    n = int(rng.choice([1, 2, 3, 4, 5, 8, 12, 20, 32, 48]))
    nrhs = int(rng.choice([1, 1, 2, 4]))
    kind_id = idx % 4 + 1
    cplx_c = kind_id in (2, 4)
    cplx_b = kind_id in (2, 3)
    explicit_row = idx % 3 != 0
    scale = float([1e-80, 1e-24, 1, 1e24, 1e80][idx % 5])

    c = rng.normal(size=n)
    if cplx_c:
        c = c + 1j * rng.normal(size=n)
    c[0] = (n + 3.0) + (0.4j if cplx_c and idx % 2 else 0.0)
    if explicit_row:
        r = rng.normal(size=n)
        if cplx_c:
            r = r + 1j * rng.normal(size=n)
        r[0] = -123 + (0.27j if cplx_c else 0)
    else:
        r = None

    # These valid nonsingular cases force a zero leading principal pivot.
    force_fallback = idx % 23 == 7 and n >= 2
    if force_fallback:
        c[0] = 0.0
        if r is None:
            r = np.conjugate(c).copy()
            explicit_row = True
        c[1] += 0.5
        r[1] += 0.4
    c *= scale
    if r is not None:
        r *= scale
    t = toeplitz(c, r)
    x_ref = rng.normal(size=(n,nrhs))
    if cplx_b:
        x_ref = x_ref + 1j * rng.normal(size=(n,nrhs))
    b = t @ x_ref
    # Solve general complex systems with real b when requested.
    if kind_id == 4:
        b = b.real
    return kind_id, c, r, b, x_ref, force_fallback


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=1600)
    parser.add_argument("--seed", type=int, default=20261109)
    parser.add_argument("--mode", choices=["checked", "release"], default="checked")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    cases = [make_case(rng, idx) for idx in range(args.cases)]
    lines = []
    for kind_id, c, r, b, _, _ in cases:
        n, nrhs = b.shape
        lines.append(f"{n} {nrhs} {kind_id} {int(r is not None)}")
        lines.append(numbers(c, kind_id in (2,4)))
        if r is not None:
            lines.append(numbers(r, kind_id in (2,4)))
        lines.append(numbers(b, kind_id in (2,3)))
    with tempfile.TemporaryDirectory(prefix="scifort_toeplitz_") as directory:
        exe = compile_driver(Path(directory), args.mode)
        proc = subprocess.run([str(exe)], input="\n".join(lines) + "\n",
                              text=True, capture_output=True)
        if proc.returncode:
            raise RuntimeError(f"Fortran failed: {proc.stderr[-3000:]}")
    rows = iter(proc.stdout.splitlines())
    max_relative_solution = 0.0
    max_relative_residual = 0.0
    fallbacks = 0
    worst = None
    for idx, (kind_id, c, r, b, x_ref, forced) in enumerate(cases):
        status, dense = (int(z) for z in next(rows).split())
        if status:
            raise AssertionError(f"case {idx}: status={status}")
        if forced and not dense:
            raise AssertionError(f"case {idx}: expected dense fallback")
        fallbacks += dense
        n,nrhs = b.shape
        z = np.empty((n,nrhs), dtype=complex if kind_id != 1 else float)
        for j in range(nrhs):
            for i in range(n):
                value = next(rows).split()
                z[i,j] = complex(*map(float,value)) if len(value) == 2 else float(value[0])
        matrix = toeplitz(c,r)
        dense_ref = np.linalg.solve(matrix,b)
        if not forced:
            scipy_ref = scipy_solve_toeplitz((c,r) if r is not None else c, b)
            # Some near-degenerate cases can amplify differences, so compare
            # all three independently, not just against the constructed x_ref.
            if not np.allclose(scipy_ref,dense_ref,rtol=5e-8,atol=5e-10):
                raise AssertionError(f"SciPy and dense oracles disagree for {idx}")
        diff = np.max(np.abs(z-dense_ref))
        err = float(diff / max(1.0,np.max(np.abs(dense_ref))))
        residual = np.max(np.abs(matrix@z-b))
        den = max(1.0,float(np.max(np.abs(matrix))*np.max(np.abs(z))*n + np.max(np.abs(b))))
        rel_res = float(residual / den)
        if err > max_relative_solution:
            max_relative_solution, worst = err,idx
        max_relative_residual = max(max_relative_residual,rel_res)
        if err > 4e-9 or rel_res > 4e-9:
            raise AssertionError(f"case {idx}: rel err {err:.3e}, residual {rel_res:.3e}")
    try:
        extra = next(rows)
    except StopIteration:
        extra = None
    if extra is not None:
        raise AssertionError("Unexpected extra driver output")
    print(f"PASS {args.cases}/{args.cases} comparisons ({args.mode}, SciPy {scipy.__version__})")
    print(f"Pivoted-LU fallbacks: {fallbacks}")
    print(f"Max relative solution error: {max_relative_solution:.3e} (case {worst})")
    print(f"Max normalized residual: {max_relative_residual:.3e}")


if __name__ == "__main__":
    main()
