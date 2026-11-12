#!/usr/bin/env python3
"""Compile SciFort's 2-D spline routines and compare with SciPy.

Usage: python tools/signal_reference/compare_signal_spline2d_scipy.py
Requires gfortran, NumPy and SciPy; uses no FPM dependencies.
"""
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
    for degree, lam, precision in [
        (2, 0.0, -1.0), (3, 0.0, -1.0), (3, 0.0, 1e-10),
        (3, 0.1, -1.0), (3, 1.0, -1.0), (3, 5.0, -1.0),
        (3, 0.1, 1e-9), (3, 5.0, 1e-10),
        (4, 0.0, -1.0), (4, 0.1, -1.0), (4, 5.0, -1.0),
    ]:
        for m, n in [(74, 79), (86, 88)]:
            for complex_input in (False, True):
                x = RNG.normal(size=(m,n))
                # Use correlated images and nonzero-mean images as well.
                if len(cases) % 2:
                    x = x.cumsum(axis=0) * 0.04 + 1
                if complex_input:
                    x = x + 1j*RNG.normal(size=(m,n))
                if degree == 2:
                    expected = signal.qspline2d(x, lamb=lam, precision=precision)
                elif degree == 3:
                    if complex_input and lam > 1 / 144:
                        # SciPy symiirorder2 rejects complex input; the Fortran
                        # implementation extends it by linearity.
                        expected = (signal.cspline2d(x.real,lamb=lam,precision=precision)
                                    + 1j*signal.cspline2d(x.imag,lamb=lam,precision=precision))
                    else:
                        expected = signal.cspline2d(x,lamb=lam,precision=precision)
                else:
                    if complex_input:
                        # SciPy demotes complex128 to complex64 internally
                        # (and rejects complex smoothed coefficients). Compare
                        # against componentwise double-precision reference.
                        expected = (signal.spline_filter(x.real,lmbda=lam)
                                    + 1j*signal.spline_filter(x.imag,lmbda=lam))
                    else:
                        expected = signal.spline_filter(x,lmbda=lam)
                cases.append((degree,lam,precision,x,expected))

    sources = [ROOT/'src'/f'{p}.f90' for p in
               ('kinds','signal_spatial','signal_spline','signal_spline2d')]
    driver = ROOT/'tools/signal_reference/spline2d_reference_driver.f90'
    with tempfile.TemporaryDirectory(prefix='scifort-spline2d-') as temp:
        temp = Path(temp)
        objs = []
        for source in sources:
            obj = temp/(source.stem+'.o')
            subprocess.run(['gfortran','-std=f2018','-O2','-fcheck=bounds',
                            '-J',str(temp),'-I',str(temp),'-c',str(source),
                            '-o',str(obj)],check=True,cwd=temp)
            objs.append(str(obj))
        exe = temp/'spline2d_reference_driver'
        subprocess.run(['gfortran','-std=f2018','-O2','-J',str(temp),'-I',str(temp),
                        str(driver),*objs,'-o',str(exe)],check=True,cwd=temp)
        lines = [str(len(cases))]
        for degree,lam,precision,x,_ in cases:
            m,n=x.shape
            lines.append(f'{m} {n} {degree} {int(np.iscomplexobj(x))} {lam:.17g} {precision:.17g}')
            lines.extend(f'{z.real:.17e} {z.imag:.17e}' for z in x.T.flat)
        text = '\n'.join(lines)+'\n'
        proc = subprocess.run([str(exe)],input=text,check=True,text=True,
                              capture_output=True,cwd=temp)
        outputs = iter(proc.stdout.splitlines())
        worst = 0.0
        count = 0
        for degree,lam,precision,x,ref in cases:
            status, m, n = map(int,next(outputs).split())
            if status != 0 or (m,n) != x.shape:
                raise AssertionError(f'status={status}, shape={(m,n)}, expected={x.shape}')
            raw = np.array([complex(*map(float,next(outputs).split())) for _ in range(m*n)])
            got = raw.reshape((n,m)).T
            err = float(np.max(abs(got-ref)))
            tol = 2e-5 if degree == 4 else 5e-7
            if err > tol:
                raise AssertionError(f'{degree=} {lam=} {precision=} {m=} {n=}, complex={np.iscomplexobj(x)}, max input={abs(x).max():.4g}: {err=}, {tol=}')
            worst = max(worst,err)
            count += 1
        print(f'PASS: {count}/{count} SciPy spline2d reference cases, max abs error {worst:.3g}')


if __name__ == '__main__':
    main()
