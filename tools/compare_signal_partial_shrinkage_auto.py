#!/usr/bin/env python3
"""Compare complex-Gaussian diagonal shrinkage and partial coherency with NumPy.

The reference independently implements the complex Wishart moment plug-in,
regularized precision inverse, and conditional coherency. Deterministic random
spectral matrices are generated from independent complex or real observations;
real cases only exercise algebraic compatibility, not proper-complex statistics.
"""
from __future__ import annotations
import argparse
import io
import subprocess
import numpy as np


def alpha_numpy(s: np.ndarray, n: float) -> float:
    if n <= 1 + np.finfo(float).eps:
        return 1.0
    scaled = s / np.max(np.abs(s))
    d = np.real(np.diag(scaled))
    mask = ~np.eye(s.shape[0], dtype=bool)
    diag_prod = float(np.sum((d[:,None]*d[None,:])[mask]))
    off = float(np.sum((np.abs(scaled)**2)[mask]))
    if off <= 1.e-30:
        return 0.0
    return float(np.clip((n*diag_prod-off) / ((n-1)*(n+1)*off), 0, 1))


def create_cases() -> list[tuple[list[np.ndarray], float, float, int, int, int]]:
    rng = np.random.default_rng(84137)
    cases = []
    for idx in range(196):
        nc = 2 + idx % 5
        nf = 1 + (idx//3) % 4
        nt = 1 + (idx//7) % 3
        n_eff = [1., 2., 3., 4., 6., 9., 15., 25., 4.5][idx % 9]
        ridge = [0.015, 0.08, 0.0][idx % 3]
        # Avoid numerical near-singularity in no-ridge cases: independent
        # samples >= channels + 3, and positive-definite covariance.
        if ridge == 0.0 and n_eff < nc+3:
            ridge = 0.025
        mats = []
        for step in range(nf*nt):
            if idx % 16 == 0 and step == 0:
                s = np.diag(rng.uniform(.4, 2.0, size=nc)).astype(complex)
            elif idx % 23 == 0 and step == nf*nt-1:
                s = np.zeros((nc, nc), dtype=complex)
            else:
                n = max(1, int(np.ceil(n_eff)))
                x = rng.standard_normal((n, nc))
                if idx % 4 != 0:
                    x = (x + 1j * rng.standard_normal((n, nc))) / np.sqrt(2)
                mix = np.eye(nc) + (rng.standard_normal((nc, nc)) +
                                         1j*rng.standard_normal((nc,nc))) * 0.2
                if idx % 4 == 0:
                    mix = mix.real
                x = x @ mix
                s = (x.conj().T @ x)/n
                scale = 10.0**rng.uniform(-1.0, 1.0)
                s *= scale
                s = (s + s.conj().T)/2
            mats.append(s)
        cases.append((mats, n_eff, ridge, nf, nt, nc))
    return cases


def main(exe: str) -> None:
    cases = create_cases()
    stream = io.StringIO()
    stream.write(f'{len(cases)}\n')
    for mats, n_eff, ridge, nf, nt, nc in cases:
        stream.write(f'{nf} {nt} {nc} {n_eff:.17g} {ridge:.17g}\n')
        for mat in mats:
            for j in range(nc):
                for i in range(nc):
                    z = mat[i, j]
                    stream.write(f'{z.real:.18e} {z.imag:.18e}\n')
    cmd = subprocess.run([exe], input=stream.getvalue(), text=True,
                         capture_output=True, check=True)
    tokens = iter(cmd.stdout.split())
    nreturned = int(next(tokens))
    assert nreturned == len(cases)
    max_error = 0.0
    total_valid = total_singular = total_elements = 0
    for idx, (mats, n_eff, ridge, nf, nt, nc) in enumerate(cases):
        status, nf2, nt2, nc2 = [int(next(tokens)) for _ in range(4)]
        assert (nf2, nt2, nc2) == (nf, nt, nc)
        count_singular = 0
        for m, mat in enumerate(mats):
            valid = bool(int(next(tokens)))
            alpha = float(next(tokens))
            is_zero = not np.any(mat)
            if is_zero:
                assert not valid, f'case {idx} matrix {m}: zero spectrum accepted'
                assert np.isnan(alpha)
                count_singular += 1
                total_singular += 1
                continue
            alpha_expected = alpha_numpy(mat, n_eff)
            max_error = max(max_error, abs(alpha-alpha_expected))
            assert abs(alpha-alpha_expected) <= 2.0e-12, (idx,m,alpha,alpha_expected)
            reg = (1-alpha_expected)*mat + alpha_expected*np.diag(np.diag(mat))
            reg += ridge*np.trace(mat).real / nc * np.eye(nc)
            # All constructed nonzero matrices must be invertible under the
            # selected ridge or alpha, except a rank-deficient clipped-alpha=0.
            try:
                eig = np.linalg.eigvalsh(reg)
                expected_valid = min(eig)>64*np.finfo(float).eps*np.max(np.abs(reg))
            except np.linalg.LinAlgError:
                expected_valid = False
            assert valid == expected_valid, (idx, m, n_eff, alpha, eig, ridge)
            if not valid:
                count_singular += 1
                total_singular += 1
                continue
            p = np.linalg.inv(reg)
            coh = -p/np.sqrt(np.real(np.diag(p))[:,None]*np.real(np.diag(p))[None,:])
            np.fill_diagonal(coh, 1.)
            power = np.abs(coh)**2
            for j in range(nc):
                for i in range(nc):
                    reported = [float(next(tokens)) for _ in range(5)]
                    expected = [power[i,j], p[i,j].real, p[i,j].imag,
                                coh[i,j].real, coh[i,j].imag]
                    for a, b in zip(reported, expected):
                        max_error = max(max_error, abs(a-b)/(1+abs(b)))
                        assert abs(a-b) <= 5e-9*(1+abs(b)), (idx,m,i,j,a,b)
                    total_elements += 5
            total_valid += 1
        assert status == (2 if count_singular else 0), (idx, status, count_singular)
    assert next(tokens, None) is None, 'unexpected trailing output'
    print(f'PASS {len(cases)} NumPy cases; {total_valid} invertible bins; '
          f'{total_singular} singular bins; {total_elements} precision/coherency values; '
          f'max scaled error {max_error:.3g}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    main(parser.parse_args().driver)
