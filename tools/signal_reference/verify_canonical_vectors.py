#!/usr/bin/env python3
"""Independent NumPy SVD validation for conditional canonical spectral directions.

Compare variance normalization and cross-covariance diagonalization for the
actual Fortran coefficients, and compare isolated mode directions (up to phase)
against NumPy SVD. The implementation of SciFort's Jacobi solver is not reused.
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np


def conditional(s, aa, bb, cc, ridge, shrink):
    n = len(s)
    reg = (1 - shrink) * s.copy()
    reg[np.diag_indices(n)] = np.diag(s) + ridge * np.trace(s).real / n

    def sub(i, j):
        return reg[np.ix_(i, j)]

    sa, sb, sab = sub(aa, aa), sub(bb, bb), sub(aa, bb)
    if len(cc):
        saa = sub(aa, cc)
        sbb = sub(bb, cc)
        sc = sub(cc, cc)
        sa = sa - saa @ np.linalg.solve(sc, saa.conj().T)
        sb = sb - sbb @ np.linalg.solve(sc, sbb.conj().T)
        sab = sab - saa @ np.linalg.solve(sc, sbb.conj().T)
    la, lb = np.linalg.cholesky(sa), np.linalg.cholesky(sb)
    w = np.linalg.solve(la, sab)
    w = np.linalg.solve(lb, w.conj().T).conj().T
    ua, sing, vh = np.linalg.svd(w, full_matrices=False)
    refa = np.linalg.solve(la.conj().T, ua)
    refb = np.linalg.solve(lb.conj().T, vh.conj().T)
    return sa, sb, sab, sing, refa, refb


def main(binary):
    rng = np.random.default_rng(291049)
    n_cases = 150
    chunks = [str(n_cases)]
    cases = []
    for case in range(n_cases):
        nc = int(rng.integers(3, 10))
        nf = int(rng.integers(1, 5))
        nt = int(rng.integers(1, 4))
        perm = rng.permutation(nc)
        na = int(rng.integers(1, min(5, nc)))
        nb = int(rng.integers(1, min(5, nc-na)+1))
        nd = int(rng.integers(0, min(3, nc-na-nb)+1))
        aa = perm[:na]
        bb = perm[na:na+nb]
        cc = perm[na+nb:na+nb+nd]
        ridge = [0, 1e-5, 0.01, 0.1][case % 4]
        shrink = [0, 0.03, 0.25, 0.65, 0.95][case % 5]
        chunks.append(f'{nc} {nf} {nt} {na} {nb} {nd} {ridge:.17e} {shrink:.17e}')
        chunks.append(' '.join(str(q+1) for q in aa))
        chunks.append(' '.join(str(q+1) for q in bb))
        if nd:
            chunks.append(' '.join(str(q+1) for q in cc))
        expected = []
        for t in range(nt):
            for f in range(nf):
                z = rng.normal(size=(nc, nc+2))
                if case % 2 == 0:
                    z = z + 1j*rng.normal(size=z.shape)
                z *= np.exp(rng.uniform(-2, 2, nc))[:, None]
                s = z @ z.conj().T + 0.05*np.eye(nc)
                if case % 7 == 0:
                    s += 20 * np.outer(z[:, 0], z[:, 0].conj())
                expected.append(conditional(s, aa, bb, cc, ridge, shrink))
                for j in range(nc):
                    for i in range(nc):
                        v = s[i,j]
                        chunks.append(f'{v.real:.17e} {v.imag:.17e}')
        cases.append((nf, nt, na, nb, expected))
    process = subprocess.run([binary], input='\n'.join(chunks)+'\n', text=True,
                             capture_output=True, check=True,
                             env=dict(os.environ, OPENBLAS_NUM_THREADS='1',
                                      OMP_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    rows = iter(process.stdout.splitlines())
    worst = 0.0
    n_cells = 0
    n_isolated = 0
    n_nonisolated = 0
    for case, (nf, nt, na, nb, expected) in enumerate(cases):
        status, got_nf, got_nt = map(int, next(rows).split())
        assert (status, got_nf, got_nt) == (0, nf, nt), (case, status, nf, nt)
        for index, (sa, sb, sab, singular, refa, refb) in enumerate(expected):
            valid = next(rows).strip()
            assert valid == 'T', (case, index, valid)
            nm = min(na, nb)
            a = np.empty((na, nm),dtype=complex)
            b = np.empty((nb, nm),dtype=complex)
            modes = np.empty(nm)
            stable = np.empty(nm,dtype=bool)
            for k in range(nm):
                parts = next(rows).split()
                stable[k] = parts[0] == 'T'
                vals = np.array([float(v) for v in parts[1:]])
                assert len(vals) == 1 + 2*(na+nb)
                modes[k] = vals[0]
                a[:,k] = vals[1:1+2*na:2] + 1j*vals[2:1+2*na:2]
                b[:,k] = vals[1+2*na::2] + 1j*vals[2+2*na::2]
            eye = np.eye(nm)
            errs = (a.conj().T @ sa @ a - eye,
                    b.conj().T @ sb @ b - eye,
                    a.conj().T @ sab @ b - np.diag(singular),
                    modes - singular**2)
            err = max(float(np.max(abs(e))) for e in errs)
            worst = max(err, worst)
            if err > 3e-8:
                raise AssertionError(f'case={case}, bin={index}: invariant error {err}; '
                                     f'singular={singular}, modes={modes}')
            for k in range(nm):
                if stable[k]:
                    n_isolated += 1
                    overlap_a = abs(np.vdot(refa[:,k], sa @ a[:,k]))
                    overlap_b = abs(np.vdot(refb[:,k], sb @ b[:,k]))
                    if abs(overlap_a-1) > 1e-6 or abs(overlap_b-1) > 1e-6:
                        raise AssertionError(f'case={case}, bin={index}, mode={k}, '
                                             f'NumPy-SVD overlap {overlap_a} {overlap_b}')
                else:
                    n_nonisolated += 1
            n_cells += 1
    assert next(rows,None) is None, 'unexpected extra driver output'
    print(f'PASS {n_cases} NumPy-SVD canonical-vector scenarios, {n_cells} cells, '
          f'{n_isolated} isolated and {n_nonisolated} degenerate directions; '
          f'max invariant error {worst:.3g}',flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
