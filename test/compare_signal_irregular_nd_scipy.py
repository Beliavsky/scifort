"""Compare rank-3 multichannel Lomb-Scargle power/amplitude with SciPy.

    python test/compare_signal_irregular_nd_scipy.py checked
    python test/compare_signal_irregular_nd_scipy.py optimized

The rank-15 and noncontiguous-array cases are checked separately by the
native Fortran regression.  This driver intentionally uses the public
mathematics of scipy.signal.lombscargle as an independent numeric reference.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.signal import lombscargle

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "reference" / "signal_irregular_nd_reference_driver.f90"


def test_cases():
    return [(n, nf, nc, variant)
            for n in (6, 12, 23, 47, 86)
            for nf in (1, 5, 13, 39)
            for nc in (1, 3)
            for variant in range(6)]


def inputs(n, nf, nc):
    i = np.arange(1, n+1, dtype=float)
    x = .37*(i-1) + .009*i*i + .13*np.sin(.17*i)
    w = .2 + .15*((i.astype(int)*7) % 13)
    j = np.arange(1, nf+1, dtype=float)
    freqs = .13 + .31*j + .004*j*j
    y = np.empty((2,n,nc))
    for ch in range(1, nc+1):
        for k in range(1,3):
            y[k-1,:,ch-1] = (np.cos(.85*x+.23*k*ch)
                            + .37*np.sin(1.74*x-.13*ch)
                            + .05*np.cos(.19*i*k) + .4*ch+.13*k)
    return x, y, freqs, w


def expected(n, nf, nc, variant):
    x, y, omega, w = inputs(n, nf, nc)
    parameters = {}
    if variant in (1,2,3,4):
        parameters["weights"] = w
    if variant in (1,2,4):
        parameters["floating_mean"] = True
    if variant == 3:
        parameters["precenter"] = True
    parameters["normalize"] = ("amplitude" if variant in (2,5)
                                else "normalize" if variant == 1 else "power")
    result = np.empty((2,nf,nc), dtype=complex if variant in (2,5) else float)
    for ch in range(nc):
        for k in range(2):
            result[k,:,ch] = lombscargle(x,y[k,:,ch],omega,**parameters)
    return result


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'checked'
    if mode not in ('checked','optimized'):
        raise SystemExit('Choose checked or optimized')
    flags = ['-std=f2018','-Werror=implicit-interface']
    flags += (['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3'])
    cases = test_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-lomb-nd-') as tmp:
        work = Path(tmp)
        binary = work / 'lomb_reference'
        command = [os.environ.get('FC','gfortran'), *flags, 
                   str(ROOT/'src/kinds.f90'),
                   str(ROOT/'src/signal_irregular.f90'),
                   str(ROOT/'src/signal_irregular_nd.f90'),
                   str(DRIVER), '-o', str(binary)]
        subprocess.run(command, cwd=work, check=True, capture_output=True)
        proc = subprocess.run([str(binary)], cwd=work,
                              input=''.join(' '.join(map(str,c))+'\n' for c in cases),
                              text=True, capture_output=True, check=True)
    data = iter(proc.stdout.splitlines())
    worst_absolute = worst_scaled = 0.
    value_count = 0
    for case in cases:
        n,nf,nc,variant = case
        header = next(data).split()
        assert header == ['CASE',str(n),str(nf),str(nc),str(variant),'0'], (case,header)
        nvals = 2*nf*nc
        if variant in (2,5):
            vals = np.array([complex(*map(float,next(data).split())) for _ in range(nvals)])
        else:
            vals = np.array([float(next(data)) for _ in range(nvals)])
        got = vals.reshape((2,nf,nc), order='F')
        target = expected(*case)
        absolute = np.max(np.abs(got-target))
        scale = max(1., float(np.max(np.abs(target))))
        relative = absolute / scale
        if not np.isfinite(relative) or relative > 3.e-10:
            raise AssertionError((case, absolute, relative, target.ravel()[:4], got.ravel()[:4]))
        worst_absolute = max(worst_absolute,absolute)
        worst_scaled = max(worst_scaled,relative)
        value_count += nvals
    assert next(data, None) is None
    print(f'PASS {mode}: {len(cases)} SciPy 1.17 comparisons, {value_count} spectral values; '
          f'max absolute error={worst_absolute:.3g}, max scaled error={worst_scaled:.3g}')


if __name__ == '__main__':
    main()
