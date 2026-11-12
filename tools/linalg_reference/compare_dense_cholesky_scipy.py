#!/usr/bin/env python3
"""SciPy comparisons for reusable dense real/complex Cholesky factorization."""
from __future__ import annotations
from dataclasses import dataclass
import subprocess
from pathlib import Path

import numpy as np
from scipy import linalg


@dataclass
class Case:
    variant: int
    a: np.ndarray  # upper or lower triangle describes the SPD matrix
    positive_definite: np.ndarray
    b: np.ndarray
    b2: np.ndarray
    lower: bool
    f_reference: np.ndarray
    x_reference: np.ndarray
    x2_reference: np.ndarray


def build_cases(count: int, seed: int, max_order: int, stress: bool = False) -> list[Case]:
    rng = np.random.default_rng(seed)
    sizes = [1, 2, 3, 4, 7, 12, 20, max_order]
    if stress:
        sizes += [36, 64, 96, 128, 192]
    cases = []
    for idx in range(count):
        variant = 1 + idx % 2
        lower = bool((idx // 2) % 2)
        n = int(rng.choice(sizes if idx % 3 else [int(rng.integers(1, max_order + 1))]))
        nrhs = int(rng.choice([1, 2, 3, 5]))
        z = rng.standard_normal((n, n)) / np.sqrt(n)
        if variant == 2:
            z = z + 1j*rng.standard_normal((n, n)) / np.sqrt(n)
        positive_definite = z.conj().T @ z
        positive_definite += (0.3 + rng.random()) * np.eye(n)
        # Deliberate scale variation while keeping the condition number moderate.
        if idx % 17 == 0:
            positive_definite *= 1e-50
        elif idx % 23 == 0:
            positive_definite *= 1e50
        elif stress and idx % 9 == 0 and n > 1:
            # A more demanding spectrum, generated with a unitary similarity.
            q, _ = np.linalg.qr(z)
            eigenvalues = np.geomspace(1., 1e-7, n)
            positive_definite = (q * eigenvalues) @ q.conj().T
            positive_definite = (positive_definite + positive_definite.conj().T) / 2
        # Only the selected triangle is authoritative; corrupt the ignored half.
        a = positive_definite.copy()
        unwanted = np.triu_indices(n, 1) if lower else np.tril_indices(n, -1)
        corrupt = 10*rng.standard_normal(len(unwanted[0]))
        if variant == 2:
            corrupt = corrupt + 7j*rng.standard_normal(len(unwanted[0]))
        a[unwanted] = corrupt
        if variant == 2 and idx % 15 == 0:
            # Consistent with LAPACK POTRF: imaginary diagonal components are ignored.
            a[np.diag_indices(n)] += 1j*rng.standard_normal(n)
        truth = rng.standard_normal((n, nrhs))
        truth2 = rng.standard_normal((n, nrhs))
        if variant == 2:
            truth = truth + 1j*rng.standard_normal((n, nrhs))
            truth2 = truth2 + 1j*rng.standard_normal((n, nrhs))
        b = positive_definite @ truth
        b2 = positive_definite @ truth2
        factor, was_lower = linalg.cho_factor(a, lower=lower, check_finite=True)
        assert was_lower == lower
        ref = linalg.cho_solve((factor, lower), b)
        ref2 = linalg.cho_solve((factor, lower), b2)
        cases.append(Case(variant, a, positive_definite, b, b2, lower, factor, ref, ref2))
    return cases


def run(executable: Path, cases: list[Case], tolerance: float) -> dict[str,float]:
    lines = [str(len(cases))]
    for case in cases:
        n, nrhs = case.b.shape
        lines.append(f'{case.variant} {n} {nrhs} {int(case.lower)}')
        for array in (case.a, case.b, case.b2):
            for v in np.asarray(array).ravel(order='F'):
                lines.append(f'{np.real(v):.17e} {np.imag(v):.17e}')
    result = subprocess.run([str(executable)], input='\n'.join(lines) + '\n',
                            text=True, capture_output=True, check=True, timeout=180)
    linesout = result.stdout.splitlines()
    pos = 0
    worst_solution = 0.
    worst_factor = 0.
    worst_backward = 0.
    worst_reconstruction = 0.
    for idx, case in enumerate(cases):
        n, nrhs = case.b.shape
        if pos >= len(linesout):
            raise AssertionError(f'No output for case {idx}')
        status = int(linesout[pos])
        pos += 1
        if status:
            raise AssertionError(f'Nonzero status {status} at case {idx} '
                                 f'(n={n}, complex={case.variant==2}, lower={case.lower})')
        nfactor = n*n
        nsol = n*nrhs
        if pos + nfactor + 2*nsol > len(linesout):
            raise AssertionError(f'Truncated output at case {idx}')
        allvalues = np.array([complex(*map(float, ln.split())) for ln in
                              linesout[pos:pos+nfactor+2*nsol]])
        pos += nfactor+2*nsol
        f = allvalues[:nfactor].reshape((n,n),order='F')
        x = allvalues[nfactor:nfactor+nsol].reshape((n,nrhs),order='F')
        x2 = allvalues[nfactor+nsol:].reshape((n,nrhs),order='F')
        if case.variant == 1:
            if max(np.max(np.abs(v.imag)) for v in (f,x,x2)):
                raise AssertionError(f'Unexpected imaginary output in real case {idx}')
            f = f.real
            x = x.real
            x2 = x2.real
        mask = np.tril(np.ones((n,n),dtype=bool)) if case.lower else np.triu(np.ones((n,n),dtype=bool))
        if np.any(f[~mask]):
            raise AssertionError(f'Unused output triangle not zero at case {idx}')
        fref = np.tril(case.f_reference) if case.lower else np.triu(case.f_reference)
        ferr = np.linalg.norm(f-fref)/max(np.linalg.norm(fref), np.finfo(float).tiny)
        worst_factor = max(worst_factor, float(ferr))
        rec = f @ f.conj().T if case.lower else f.conj().T @ f
        recerr = np.linalg.norm(rec-case.positive_definite)/np.linalg.norm(case.positive_definite)
        worst_reconstruction = max(worst_reconstruction,float(recerr))
        for obtained, ref, b in ((x, case.x_reference, case.b),
                                 (x2, case.x2_reference, case.b2)):
            err = np.linalg.norm(obtained-ref) / max(1.,np.linalg.norm(ref))
            residual = np.linalg.norm(case.positive_definite @ obtained-b)
            denominator = (np.linalg.norm(case.positive_definite)*np.linalg.norm(obtained) +
                           np.linalg.norm(b))
            backward = residual/max(np.finfo(float).tiny,denominator)
            worst_solution=max(worst_solution,float(err))
            worst_backward=max(worst_backward,float(backward))
            if err > tolerance or backward > max(1e-12, tolerance):
                raise AssertionError(f'Case {idx}, n={n} gave error={err:.3e}, '
                                     f'backward={backward:.3e}, tolerance={tolerance:.3e}')
        if ferr > max(2e-11,tolerance) or recerr > max(2e-11,tolerance):
            raise AssertionError(f'Factor mismatch case {idx}: ferr={ferr:.3e}, rec={recerr:.3e}')
    if pos != len(linesout):
        raise AssertionError('Unexpected surplus Fortran output')
    print(f'passed {len(cases)}/{len(cases)} SciPy cases; '
          f'worst solution={worst_solution:.3e}, factor={worst_factor:.3e}, '
          f'reconstruction={worst_reconstruction:.3e}, backward={worst_backward:.3e}')
    return {'solution':worst_solution, 'factor':worst_factor,
            'reconstruction':worst_reconstruction, 'backward':worst_backward}
