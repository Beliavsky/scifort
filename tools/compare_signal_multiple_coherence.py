#!/usr/bin/env python3
"""Independent NumPy inversion and Schur complement reference for spectral R^2.

Tests 700 frequency/time/channel/configuration output rows against an
independently constructed Hermitian PSD spectral matrix. NumPy is test-only.
"""
from __future__ import annotations
import math
import subprocess
import sys
import numpy as np

NF, NT, NC = 7, 4, 5

def spectrum(f: int, t: int) -> np.ndarray:
    z = np.empty((NC, NC), np.complex128)
    for j in range(1, NC + 1):
        for i in range(1, NC + 1):
            a = 0.13*i + 0.17*j + 0.11*f + 0.07*t
            z[i-1, j-1] = (math.sin(a) + (0.17*i if i == j else 0.0)
                            + 1j*(math.cos(2*a) + 0.11*j))
    s = z.conj().T @ z
    s[np.diag_indices(NC)] += [0.15 + 0.035*i for i in range(1, NC+1)]
    return s


def moment_shrinkage(s: np.ndarray, n: float) -> float:
    # Fixed-weight Gaussian-Wishart plug-in used in the signal module.
    scale = np.abs(s).max()
    d = sum((s[i,i].real/scale)*(s[j,j].real/scale)
            for i in range(NC) for j in range(NC) if i != j)
    o = sum(abs(s[i,j]/scale)**2
            for i in range(NC) for j in range(NC) if i != j)
    if o == 0:
        return 0.0
    return float(np.clip((n*d-o)/((n*n-1)*o), 0.0, 1.0))


def check(exe: str) -> None:
    out = subprocess.run([exe], check=True, text=True, capture_output=True).stdout
    worst = 0.0
    count = 0
    for line in out.splitlines():
        mode, f, t, i, coh, residual, explained, shrink = line.split()
        mode, f, t, i = map(int, (mode, f, t, i))
        coh, residual, explained, shrink = map(float,(coh,residual,explained,shrink))
        raw = spectrum(f,t)
        ridge, alpha = {1:(0.0,0.0),2:(0.04,0.2),3:(0.10,0.8),
                        4:(0.05,0.04*f+0.03*t),5:(0.02,moment_shrinkage(raw,5.4))}[mode]
        reg = raw.copy()
        for k in range(NC):
            for j in range(NC):
                if k != j:
                    reg[k,j] *= 1-alpha
        reg += ridge * np.trace(raw).real/NC * np.eye(NC)
        # Schur complement is independent of the Fortran precision algorithm.
        idx = i-1
        keep = np.array([j for j in range(NC) if j != idx])
        cond = (reg[idx,idx]-reg[idx,keep] @
                np.linalg.solve(reg[np.ix_(keep,keep)],reg[keep,idx])).real
        ref = np.array([1-cond/reg[idx,idx].real,cond,reg[idx,idx].real-cond,alpha])
        actual = np.array([coh,residual,explained,shrink])
        errors = abs(ref-actual)/(1+abs(ref))
        worst = max(worst, float(errors.max()))
        if not np.allclose(actual,ref,atol=2e-11,rtol=2e-11):
            raise AssertionError(f'Mismatch mode={mode} f={f} t={t} ch={i}: {actual} != {ref}')
        count += 1
    assert count == 5*NF*NT*NC, count
    print(f'PASS independent NumPy Schur complements: {count} cases, maximum normalized error {worst:.3g}')


if __name__ == '__main__':
    check(sys.argv[1])
