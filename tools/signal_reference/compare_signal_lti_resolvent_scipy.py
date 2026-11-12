#!/usr/bin/env python3
"""Independent SciPy linear-solve checks for SciFort LTI state-space responses.

Pass --driver /path/to/executable to reuse a compiled test probe. Without
--driver, this script compiles SciFort's complete standalone signal module set.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg, signal

from run_signal_suite import MODULES, ROOT


def generate_cases():
    rng = np.random.default_rng(20261018)
    cases = []
    for mode in (1, 2, 3, 4):
        digital = mode in (2, 4)
        is_complex = mode in (3, 4)
        for index in range(140):
            n = int(rng.choice([0, 1, 2, 3, 5, 8, 12, 20, 36, 60]))
            if n == 0:
                a = np.empty((0, 0), complex)
                b = np.empty((0, 1), complex)
                c = np.empty((1, 0), complex)
            else:
                if digital:
                    diagonal = rng.uniform(-0.83, 0.83, n)
                else:
                    diagonal = -rng.uniform(0.3, 5, n)
                if index % 29 == 0:
                    # Coefficients of high-order transfer polynomials overflow.
                    diagonal = (-rng.uniform(1e8, 1e10, n) if not digital
                                else rng.uniform(0.80, 0.90, n))
                a = np.diag(diagonal).astype(complex)
                a += rng.normal(scale=0.007, size=(n, n))
                if is_complex:
                    a += 1j*rng.normal(scale=0.007, size=(n, n))
                b = rng.normal(size=(n, 1)).astype(complex)
                c = rng.normal(size=(1, n)).astype(complex)
                if is_complex:
                    b += 1j*rng.normal(size=(n, 1))
                    c += 1j*rng.normal(size=(1, n))
            d = rng.normal(size=(1, 1)).astype(complex)
            if is_complex:
                d += 1j*rng.normal(size=(1, 1))
            w = (rng.uniform(-1, 1, 10) if digital else
                 rng.uniform(-7, 7, 10))
            if not digital and index % 29 == 0:
                w = rng.uniform(1e7, 1e10, 10)
            cases.append((mode, a, b, c, d, w))
    return cases


def create_input(cases):
    data = []
    for mode, a, b, c, d, w in cases:
        n = a.shape[0]
        data.append(f'{mode} {n} {len(w)}')
        for array in (a.ravel(order='F'), b.ravel(order='F'),
                      c.ravel(order='F'), d.ravel(order='F')):
            data.extend(f'{value.real:.17g} {value.imag:.17g}' for value in array)
        data.extend(f'{value:.17g}' for value in w)
    return '\n'.join(data)+'\n'


def check(driver):
    cases = generate_cases()
    process = subprocess.run([str(driver)],input=create_input(cases),
                             text=True,capture_output=True,check=True,timeout=180)
    output = iter(process.stdout.splitlines())
    max_scaled = 0.0
    nvalues = 0
    for number, (mode, a, b, c, d, w) in enumerate(cases):
        status, count = map(int,next(output).split())
        if status != 0 or count != len(w):
            raise AssertionError((number, 'bad Fortran status', status, count))
        values = np.array([list(map(float,next(output).split()))
                           for _ in range(count)])
        actual = values[:, 1]+1j*values[:, 2]
        n = a.shape[0]
        reference = []
        for wi in w:
            z = np.exp(1j*wi) if mode in (2,4) else 1j*wi
            if n:
                x = linalg.solve(z*np.eye(n)-a, b, check_finite=False)
                reference.append((c @ x+d)[0,0])
            else:
                reference.append(d[0,0])
        reference = np.array(reference)
        scaled = np.max(np.abs(actual-reference)/(1.0+np.abs(reference)))
        max_scaled = max(max_scaled, float(scaled))
        nvalues += len(w)
        if scaled > 3e-11:
            raise AssertionError((number,mode,n,scaled))
        # Cross-check SciPy's signal API for selected low-order real systems.
        if n <= 5 and n and mode in (1,2) and number % 4 == 0:
            system = signal.StateSpace(a.real,b.real,c.real,d.real)
            with np.errstate(all='ignore'):
                if mode == 1:
                    _, scipy_freq = signal.freqresp(system,w=w)
                else:
                    discrete = signal.dlti(a.real,b.real,c.real,d.real,dt=1)
                    _, scipy_freq = signal.dfreqresp(discrete,w=w)
            if np.max(np.abs(scipy_freq-actual)/(1+np.abs(actual))) > 1e-8:
                raise AssertionError((number,'SciPy signal API disagrees'))
    if list(output):
        raise AssertionError('Unexpected extra Fortran output')
    print(f'PASS: {len(cases)} state-space designs, {nvalues} frequency points; '
          f'maximum scaled difference {max_scaled:.3e}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--driver',type=Path)
    args = parser.parse_args()
    if args.driver:
        check(args.driver.resolve())
        return
    with tempfile.TemporaryDirectory(prefix='scifort-lti-resolvent-') as temp:
        work = Path(temp)
        objs = []
        for module in MODULES:
            src = ROOT/'src'/f'{module}.f90'
            obj = work/f'{module}.o'
            subprocess.run(['gfortran','-O2','-std=f2018','-Wimplicit-interface',
                            '-Werror=implicit-interface','-J',str(work),'-I',str(work),
                            '-c',str(src),'-o',str(obj)],check=True,cwd=work)
            objs.append(str(obj))
        exe = work/'lti_resolvent_driver'
        subprocess.run(['gfortran','-O2',str(ROOT/'tools/signal_reference/lti_resolvent_driver.f90'),
                        '-J',str(work),'-I',str(work),*objs,'-o',str(exe)],check=True,cwd=work)
        check(exe)


if __name__ == '__main__':
    main()
