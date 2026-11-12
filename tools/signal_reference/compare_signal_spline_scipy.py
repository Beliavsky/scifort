#!/usr/bin/env python3
"""Seeded cubic/quadratic spline comparisons against SciPy 1.17."""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261021)


def main():
    cases = []
    for degree in (2, 3):
        for complex_input in (False, True):
            for n in (1, 2, 3, 4, 8, 15, 31, 65):
                for _ in range(15):
                    lam = 0.0
                    if degree == 3 and n >= 3 and RNG.random() < 0.35:
                        lam = float(RNG.choice([0.05, 0.2, 0.6, 3.0]))
                    sig = RNG.standard_normal(n)
                    if complex_input:
                        sig = sig + 1j * RNG.standard_normal(n)
                    dx = float(RNG.choice([0.3, 0.8, 1.0, 1.7, 3.0]))
                    x0 = float(RNG.uniform(-4, 4))
                    if n == 1:
                        coords = np.array([0.0])
                    else:
                        coords = np.r_[np.arange(n).astype(float),
                                       RNG.uniform(-2*(n-1), 3*(n-1), size=17),
                                       [-2*(n-1), 3*(n-1), n-1.0]]
                    cases.append((degree, complex_input, sig, lam, dx, x0, x0 + dx*coords))
    with tempfile.TemporaryDirectory(prefix='scifort-spline-') as tmp:
        exe = Path(tmp) / 'spline_driver'
        subprocess.run(['gfortran', '-std=f2018', '-O2', '-fcheck=bounds', '-J',tmp,
                        str(ROOT/'src/kinds.f90'), str(ROOT/'src/signal_spline.f90'),
                        str(ROOT/'tools/signal_reference/spline_reference_driver.f90'),
                        '-o',str(exe)], cwd=tmp, check=True)
        lines = [str(len(cases))]
        for degree, complex_input, sig, lam, dx, x0, xx in cases:
            lines.append(f'{len(sig)} {len(xx)} {degree} {int(complex_input)} {lam:.17g} {dx:.17g} {x0:.17g}')
            lines += [f'{np.real(s):.17g} {np.imag(s):.17g}' for s in sig]
            lines += [f'{x:.17g}' for x in xx]
        run = subprocess.run([str(exe)], input='\n'.join(lines)+'\n', text=True,
                             capture_output=True, check=True, cwd=tmp)
    output = iter(run.stdout.splitlines())
    max_coeff_error = 0.0
    max_eval_error = 0.0
    for index, (degree, complex_input, sig, lam, dx, x0, xx) in enumerate(cases):
        status, n, m = map(int, next(output).split())
        assert (status,n,m) == (0,len(sig),len(xx)), (index,status,n,m)
        coeff = np.array([complex(*map(float,next(output).split())) for _ in range(n)])
        vals = np.array([complex(*map(float,next(output).split())) for _ in range(m)])
        if degree == 3:
            refc = signal.cspline1d(sig,lamb=lam)
            refv = signal.cspline1d_eval(refc, xx, dx=dx,x0=x0)
        else:
            refc = signal.qspline1d(sig,lamb=lam)
            if complex_input:
                # SciPy 1.17 currently allocates a real temporary for complex qspline evaluation.
                refv = (signal.qspline1d_eval(refc.real, xx, dx=dx, x0=x0) +
                        1j*signal.qspline1d_eval(refc.imag, xx, dx=dx, x0=x0))
            else:
                refv = signal.qspline1d_eval(refc, xx, dx=dx,x0=x0)
        errc = float(np.max(np.abs(refc-coeff)))
        errv = float(np.max(np.abs(refv-vals)))
        max_coeff_error = max(max_coeff_error,errc)
        max_eval_error = max(max_eval_error,errv)
        if errc > 3e-10 or errv > 5e-10:
            raise AssertionError(f'case {index}, degree={degree}, n={n}, complex={complex_input}, '
                                 f'lam={lam}, coeff_err={errc:.3g}, eval_err={errv:.3g}')
    print(f'PASS: {len(cases)} spline coefficient/evaluation cases; '
          f'max coefficient error {max_coeff_error:.3e}; max evaluation error {max_eval_error:.3e}')


if __name__ == '__main__':
    main()
