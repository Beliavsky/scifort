#!/usr/bin/env python3
"""Compare SciFort half-open eigenvalue intervals with scipy.linalg.eigh.

Tests standard and generalized problems, real/complex matrices, repeated spectra,
empty intervals, infinite bounds, upper/lower triangles, and eigenpair residuals.
Use OPENBLAS_NUM_THREADS=1 to keep small LAPACK calls reproducibly fast.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    'src/kinds.f90',
    'src/linalg.f90',
    'src/linalg_banded.f90',
    'src/linalg_hessenberg.f90',
    'src/linalg_banded_tridiagonal.f90',
    'src/linalg_banded_eigen.f90',
    'src/linalg_dense.f90',
    'src/linalg_dense_eigen.f90',
]
ENV = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1',
       'MKL_NUM_THREADS': '1'}


def build(work: Path, mode: str, native: bool) -> Path:
    flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace',
              '-ffpe-trap=invalid,zero,overflow'] if mode == 'checked' else ['-O3'])
    source = ('test/test_linalg_eigen_by_value.f90' if native
              else 'tools/linalg_reference/driver_eigen_by_value.f90')
    executable = work / f'{"native" if native else "driver"}_{mode}'
    command = ['gfortran', '-std=f2018', *flags, *[str(ROOT / x) for x in SOURCES],
               str(ROOT / source), '-o', str(executable)]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True,
                            env=ENV, timeout=90)
    if result.returncode:
        raise RuntimeError('Compilation failed:\n' + result.stderr[-5000:])
    return executable


def stable_norm(x):
    """Stable Frobenius/Euclidean norm for extremely scaled arrays."""
    maximum = float(np.max(np.abs(x), initial=0.0))
    if maximum == 0.0:
        return 0.0
    return float(maximum * np.linalg.norm(np.asarray(x) / maximum))


def generate(count: int, max_order: int, seed: int):
    rng = np.random.default_rng(seed)
    sizes = [n for n in sorted({1, 2, 3, 4, 5, 7, 10, 16, 24, max_order})
             if n <= max_order]
    cases = []
    for k in range(count):
        n = int(rng.choice(sizes))
        is_complex = k % 2 == 1
        lower = (k // 2) % 2 == 1
        typ = (k // 4) % 4  # 0=standard; 1,2,3=generalized
        vectors = (k // 3) % 4 != 0
        family = (k // 7) % 12
        a = rng.standard_normal((n, n))
        x = rng.standard_normal((n, n))
        if is_complex:
            a = a + 1j*rng.standard_normal((n,n))
            x = x + 1j*rng.standard_normal((n,n))
        a = (a + a.conj().T)/2
        b = (x.conj().T @ x)/max(n,1) + np.eye(n)
        if family == 1:
            a = np.diag(rng.choice([-3., 0., 2., 5.], size=n)).astype(a.dtype)
        elif family == 2:
            a = np.eye(n, dtype=a.dtype)
        elif family == 3:
            a = np.zeros_like(a)
        elif family == 4:
            b = np.diag(np.geomspace(.1,10.,n)).astype(a.dtype)
        elif family == 5:
            a *= 1e-80
            b *= 1e-80 if typ else 1.0
        elif family == 6:
            a *= 1e80
            b *= 1e80 if typ else 1.0
        elif family == 7:
            q, _ = np.linalg.qr(x)
            a = (q * np.linspace(1.0, 1.0+1e-5, n)) @ q.conj().T
        elif family == 8:
            a = np.diag(np.linspace(-4.0, 5.0, n)).astype(a.dtype)
        elif family == 9:
            b = np.eye(n)
            a *= .001
        elif family == 10:
            a = a @ a.conj().T / max(n,1)
        elif family == 11:
            b = np.diag(np.geomspace(.01,100.,n)).astype(a.dtype)
        a = np.asarray(a, dtype=np.complex128 if is_complex else np.float64)
        b = np.asarray(b, dtype=a.dtype)
        scipy_kw = {'b': b, 'type': typ} if typ else {}
        full = linalg.eigvalsh(a, **scipy_kw)
        # Construct boundaries away from floating-point-sensitive eigenvalues.
        scale = max(float(np.max(np.abs(full))), 1.0e-280)
        # A mathematically zero spectrum may contain subnormal LAPACK noise.
        # Avoid selecting on that noise with infinitesimal artificial bounds.
        if not np.any(np.abs(full) > 1.0e-280):
            scale = 1.0
        edges = [-np.inf]
        for left, right in zip(full[:-1], full[1:]):
            if right-left > 1e-7*scale:
                midpoint = float(left/2.0 + right/2.0)
                if left < midpoint < right:
                    edges.append(midpoint)
        edges.append(np.inf)
        style = k % 8
        if style == 0:
            lo, hi = -np.inf, np.inf
        elif style == 1:
            lo, hi = float(full[-1]) + scale, float(full[-1]) + 2.0*scale
        elif style == 2:
            lo, hi = float(full[0]) - 2.0*scale, float(full[0]) - scale
        elif style == 3 and family in (1, 2, 3, 8) and typ == 0:
            # Exact equality: lower bound excludes repeated diagonal eigenvalues.
            lo = float(np.min(full))
            hi = float(np.max(full)) + max(scale, 1.0)
        else:
            first = int(rng.integers(0,len(edges)-1))
            last = int(rng.integers(first+1,len(edges)))
            lo, hi = edges[first], edges[last]
        raw_a = a.copy()
        raw_b = b.copy()
        if k % 5 == 1:
            ignored = np.triu_indices(n,1) if lower else np.tril_indices(n,-1)
            raw_a[ignored] = np.nan
            raw_b[ignored] = np.nan
        cases.append(dict(a=a,b=b,raw_a=raw_a,raw_b=raw_b,is_complex=is_complex,
                          lower=lower,typ=typ,vectors=vectors,lo=lo,hi=hi))
    return cases


def compare(executable: Path, cases: list, tolerance: float):
    lines = [str(len(cases))]
    for case in cases:
        n = len(case['a'])
        lines.append(f'{2 if case["is_complex"] else 1} {n} {int(case["lower"])} '
                     f'{case["typ"]} {int(case["vectors"])} {n} {n}')
        lines.append(f'{case["lo"]:.17e} {case["hi"]:.17e}')
        for matrix in (case['raw_a'],case['raw_b']):
            for value in matrix.ravel(order='F'):
                lines.append(f'{value.real:.17e} {value.imag:.17e}')
    result = subprocess.run([str(executable)], input='\n'.join(lines)+'\n',
                            text=True, capture_output=True, env=ENV, timeout=150)
    if result.returncode:
        raise RuntimeError('Fortran driver failed:\n' + result.stderr[-2500:])
    output = iter(result.stdout.splitlines())
    worst = {'eigenvalues':0.0, 'residual':0.0, 'orthogonality':0.0}
    for k, case in enumerate(cases):
        n = len(case['a'])
        a, b, typ = case['a'],case['b'],case['typ']
        kwargs = {'b': b, 'type': typ} if typ else {}
        expected = linalg.eigh(a, subset_by_value=(case['lo'],case['hi']),
                               eigvals_only=not case['vectors'], **kwargs)
        reference = expected[0] if case['vectors'] else expected
        status, count = map(int,next(output).split())
        if status != 0 or count != len(reference):
            raise AssertionError(f'case {k}: status={status}, selected={count}, '
                                 f'SciPy selected={len(reference)}, '
                                 f'kind={typ}, interval={(case["lo"],case["hi"])}')
        obtained = np.array([float(next(output)) for _ in range(count)])
        # Near-zero eigenvalues require a matrix-scale denominator.  A
        # relative error against just the selected eigenvalues is unstable.
        matrix_scale = (stable_norm(a) if typ == 0 else
                        stable_norm(a) / stable_norm(b) if typ == 1 else
                        stable_norm(a) * stable_norm(b))
        denominator = max(float(stable_norm(reference)), float(matrix_scale),
                          np.finfo(float).tiny)
        error = float(stable_norm(obtained-reference))/denominator
        worst['eigenvalues'] = max(worst['eigenvalues'],error)
        if error > tolerance:
            raise AssertionError(f'case {k}: eigenvalues difference {error:.3g}, '
                                 f'SciPy {reference}, SciFort {obtained}')
        if not case['vectors'] or count == 0:
            continue
        entries = np.array([[float(x) for x in next(output).split()]
                            for _ in range(n*count)])
        z = (entries[:,0] + 1j*entries[:,1]).reshape((n,count), order='F')
        if not case['is_complex']:
            z = z.real
        if typ == 0:
            residual = a @ z - z*obtained
            norm_mat = stable_norm(a)
            gram = z.conj().T @ z
        elif typ == 1:
            residual = a @ z - (b @ z)*obtained
            norm_mat = stable_norm(a) + np.max(abs(obtained))*stable_norm(b)
            gram = z.conj().T @ b @ z
        elif typ == 2:
            residual = a @ (b @ z) - z*obtained
            norm_mat = stable_norm(a)*stable_norm(b) + np.max(abs(obtained))
            gram = z.conj().T @ b @ z
        else:
            residual = b @ (a @ z) - z*obtained
            norm_mat = stable_norm(a)*stable_norm(b) + np.max(abs(obtained))
            gram = z.conj().T @ np.linalg.solve(b,z)
        res = stable_norm(residual) / max(
            norm_mat*stable_norm(z), np.finfo(float).tiny)
        ortho = stable_norm(gram - np.eye(count))/max(np.sqrt(count), 1)
        worst['residual'] = max(worst['residual'],float(res))
        worst['orthogonality'] = max(worst['orthogonality'],float(ortho))
        if res > 20*tolerance or ortho > 20*tolerance:
            raise AssertionError(f'case {k}: residual={res:.3g} orthogonality={ortho:.3g}')
    if next(output,None) is not None:
        raise AssertionError('Unexpected trailing lines in Fortran driver output')
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['checked','optimized','both'], default='both')
    parser.add_argument('--cases', type=int, default=240)
    parser.add_argument('--max-order', type=int, default=48)
    parser.add_argument('--seed', type=int, default=20261108)
    parser.add_argument('--tolerance', type=float, default=5e-9)
    args = parser.parse_args()
    modes = ['checked','optimized'] if args.mode == 'both' else [args.mode]
    cases = generate(args.cases,args.max_order,args.seed)
    with tempfile.TemporaryDirectory(prefix='scifort-eigen-value-') as temp:
        work = Path(temp)
        for mode in modes:
            native = build(work,mode,True)
            native_run = subprocess.run([str(native)],capture_output=True,
                                        text=True,env=ENV,timeout=120)
            if native_run.returncode:
                raise RuntimeError(native_run.stdout + native_run.stderr)
            driver = build(work,mode,False)
            metrics = compare(driver,cases,args.tolerance)
            print(f'{mode}: {len(cases)} SciPy comparisons passed; '
                  f'eigenvalue={metrics["eigenvalues"]:.3e}, '
                  f'residual={metrics["residual"]:.3e}, '
                  f'orthogonality={metrics["orthogonality"]:.3e}; native passed')


if __name__ == '__main__':
    main()
