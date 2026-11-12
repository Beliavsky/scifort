#!/usr/bin/env python3
"""Independently verify Fortran group CSD coherence using NumPy canonical SVD.

This oracle does not use Wilks' determinant formula: it forms conditional
block covariances using Schur complements, whitens the blocks by Cholesky
factors, and computes squared singular values of the cross-block coherence.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import numpy as np


def expected(s: np.ndarray, aa: np.ndarray, bb: np.ndarray, cc: np.ndarray,
             ridge: float, shrinkage: float) -> tuple[float, float]:
    n = s.shape[0]
    scale = np.trace(s).real / n
    reg = (1.0 - shrinkage) * s
    reg[np.diag_indices(n)] = np.diag(s) + ridge * scale

    def sub(rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
        return reg[np.ix_(rows, cols)]

    block_a = sub(aa, aa)
    block_b = sub(bb, bb)
    cross = sub(aa, bb)
    if cc.size:
        cond = sub(cc, cc)
        block_a = block_a - sub(aa, cc) @ np.linalg.solve(cond, sub(cc, aa))
        block_b = block_b - sub(bb, cc) @ np.linalg.solve(cond, sub(cc, bb))
        cross = cross - sub(aa, cc) @ np.linalg.solve(cond, sub(cc, bb))
    la = np.linalg.cholesky(block_a)
    lb = np.linalg.cholesky(block_b)
    # la^{-1} cross lb^{-H} has singular values equal to canonical coherences.
    whitened = np.linalg.solve(la, cross)
    whitened = np.linalg.solve(lb, whitened.conj().T).conj().T
    squared = np.linalg.svd(whitened, compute_uv=False)**2
    if np.any(squared > 1.0 + 3e-12):
        raise AssertionError(f'nonphysical canonical coherence: {squared}')
    squared = np.clip(squared, 0, 1-1e-15)
    log_information = -np.log1p(-squared).sum()
    coherence = -np.expm1(-log_information)
    return float(coherence), float(log_information)


def main(binary: str) -> None:
    rng = np.random.default_rng(291047)
    cases: list[tuple] = []
    chunks: list[str] = []
    maxerr = 0.0
    n_cells = 0
    num_cases = 196
    chunks.append(str(num_cases))
    for case in range(num_cases):
        nc = int(rng.integers(3, 9))
        nf = int(rng.integers(1, 6))
        nt = int(rng.integers(1, 5))
        perm = rng.permutation(nc)
        n_a = int(rng.integers(1, min(4, nc)))
        n_b = int(rng.integers(1, min(4, nc-n_a)+1))
        n_c = int(rng.integers(0, min(3, nc-n_a-n_b)+1))
        a = perm[:n_a]
        b = perm[n_a:n_a+n_b]
        c = perm[n_a+n_b:n_a+n_b+n_c]
        ridge = [0, 0.00001, 0.005, 0.05, 0.3][case % 5]
        shrink = [0, 0.1, 0.2, 0.6, 0.95][(case // 5) % 5]
        complex_data = case % 2 == 0
        matrices = np.empty((nt, nf, nc, nc), dtype=np.complex128)
        ref = []
        chunks.append(f'{nc} {nf} {nt} {n_a} {n_b} {n_c} {ridge:.17e} {shrink:.17e}')
        chunks.append(' '.join(str(x+1) for x in a))
        chunks.append(' '.join(str(x+1) for x in b))
        if n_c:
            chunks.append(' '.join(str(x+1) for x in c))
        for t in range(nt):
            for f in range(nf):
                z = rng.normal(size=(nc, nc+2))
                if complex_data:
                    z = z + 1j*rng.normal(size=z.shape)
                # Correlated channels with nontrivial cross-spectrum phases.
                scales = np.exp(rng.uniform(-2.5, 2.5, nc))
                z = scales[:, None] * z
                s = z @ z.conj().T + np.diag((scales*0.07)**2)
                if case % 11 == 0:
                    # Nearly deterministic shared latent component.
                    latent = rng.normal(size=(nc, 1))
                    if complex_data:
                        latent = latent + 1j*rng.normal(size=(nc, 1))
                    s = s + 25 * (scales[:,None]*latent) @ (scales[:,None]*latent).conj().T
                matrices[t, f] = s
                ref.append(expected(s, a, b, c, ridge, shrink))
                for j in range(nc):
                    for i in range(nc):
                        val = s[i,j]
                        chunks.append(f'{val.real:.17e} {val.imag:.17e}')
        cases.append((nc,nf,nt,a,b,c,ridge,shrink,ref))
    stdin = '\n'.join(chunks) + '\n'
    process = subprocess.run([str(binary)], input=stdin, text=True, capture_output=True,
                             env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'),
                             check=True)
    rows = iter(process.stdout.splitlines())
    for case, (nc,nf,nt,a,b,c,ridge,shrink,ref) in enumerate(cases):
        status, got_nf, got_nt = map(int, next(rows).split())
        assert status == 0 and got_nf == nf and got_nt == nt, (case,status,got_nf,got_nt)
        for index, want in enumerate(ref):
            values = next(rows).split()
            got = np.array([float(values[0]),float(values[1])])
            assert values[2] == 'T', (case,index,values)
            diff = float(max(abs(got-np.asarray(want))))
            maxerr = max(maxerr,diff)
            n_cells += 1
            if not np.isclose(got[0],want[0],rtol=4e-11,atol=8e-12) or not np.isclose(
                    got[1],want[1],rtol=6e-11,atol=8e-12):
                raise AssertionError(f'case {case} bin {index}: Fortran={got} NumPy={want}; '
                                     f'groups {a},{b}|{c}, ridge={ridge}, shrink={shrink}')
    assert next(rows, None) is None, 'unexpected extra driver output'
    print(f'PASS {num_cases} independent NumPy canonical-SVD scenarios; '
          f'{n_cells} time-frequency cells; max absolute error {maxerr:.3g}',flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
