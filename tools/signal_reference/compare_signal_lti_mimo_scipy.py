#!/usr/bin/env python3
"""Compare independent SciPy linear solves and MIMO Bode responses.

Use --driver to reuse a compiled Fortran probe. Otherwise compile the
standalone SciFort signal module sources with gfortran, without FPM.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import linalg

from run_signal_suite import MODULES, ROOT


def generate_cases():
    rng = np.random.default_rng(20261018)
    cases = []
    for mode in (1, 2, 3, 4):
        digital = mode % 2 == 0
        is_complex = mode > 2
        for index in range(120):
            ns = int(rng.choice([0, 1, 2, 3, 5, 8, 16, 30, 60]))
            ni = int(rng.integers(1, 5))
            no = int(rng.integers(1, 5))
            if ns:
                diagonal = (rng.uniform(-0.85, 0.85, ns) if digital
                            else -rng.uniform(0.3, 6.0, ns))
                if index % 17 == 0:
                    diagonal = (rng.uniform(0.7, 0.9, ns) if digital else
                                -rng.uniform(1e8, 1e10, ns))
                a = np.diag(diagonal).astype(np.complex128)
                a += rng.normal(scale=0.004, size=(ns, ns))
                b = rng.normal(size=(ns, ni)).astype(np.complex128)
                c = rng.normal(size=(no, ns)).astype(np.complex128)
                if is_complex:
                    a += 0.004j*rng.normal(size=(ns, ns))
                    b += 1j*rng.normal(size=(ns, ni))
                    c += 1j*rng.normal(size=(no, ns))
            else:
                a = np.empty((0, 0), np.complex128)
                b = np.empty((0, ni), np.complex128)
                c = np.empty((no, 0), np.complex128)
            d = rng.normal(size=(no, ni)).astype(np.complex128)
            if is_complex:
                d += 1j*rng.normal(size=(no, ni))
            w = rng.uniform(-3.14, 3.14, 8) if digital else rng.uniform(-8.0, 8.0, 8)
            if not digital and ns and index % 17 == 0:
                w = rng.uniform(1e7, 1e10, 8)
            cases.append((mode, a, b, c, d, w))
    return cases


def make_input(cases):
    lines = []
    for mode, a, b, c, d, w in cases:
        ns, ni, no = a.shape[0], b.shape[1], c.shape[0]
        lines.append(f'{mode} {ns} {ni} {no} {len(w)}')
        for array in (a, b, c, d):
            lines.extend(f'{val.real:.17g} {val.imag:.17g}'
                         for val in array.ravel(order='F'))
        lines.extend(f'{val:.17g}' for val in w)
    return '\n'.join(lines)+'\n'


def check(driver):
    cases = generate_cases()
    proc = subprocess.run([str(driver)],input=make_input(cases),text=True,
                          capture_output=True,check=True,timeout=210)
    output = iter(proc.stdout.splitlines())
    worst_response = 0.0
    worst_magnitude = 0.0
    worst_phase = 0.0
    samples = 0
    channels = 0
    for number, (mode, a, b, c, d, w) in enumerate(cases):
        s1, s2, count = map(int,next(output).split())
        if (s1, s2, count) != (0, 0, len(w)):
            raise AssertionError((number,'bad Fortran status',s1,s2,count))
        ni, no, ns = b.shape[1], c.shape[0], a.shape[0]
        rows = np.array([list(map(float,next(output).split()))
                         for _ in range(len(w)*ni*no)])
        rows = rows.reshape((len(w),ni,no,5)).transpose(0,2,1,3)
        axis = rows[:,:,:,0]
        resp = rows[:,:,:,1]+1j*rows[:,:,:,2]
        magnitude = rows[:,:,:,3]
        phase = rows[:,:,:,4]
        reference = []
        for omega in w:
            z = np.exp(1j*omega) if mode in (2,4) else 1j*omega
            if ns:
                x = linalg.solve(z*np.eye(ns)-a,b,check_finite=False)
                reference.append(c@x+d)
            else:
                reference.append(d)
        expected = np.array(reference)
        scaled = float(np.max(np.abs(resp-expected)/(1.0+np.abs(expected))))
        worst_response = max(worst_response,scaled)
        if scaled > 4e-11:
            raise AssertionError((number,mode,ns,ni,no,'response error',scaled))
        physical_w = w/(0.2 if mode in (2,4) else 1.0)
        if np.max(np.abs(axis-physical_w[:,None,None])) > 5e-9:
            raise AssertionError((number,'frequency axis'))
        expected_mag = 20*np.log10(np.abs(expected))
        expected_phase = np.unwrap(np.angle(expected),axis=0)*180/np.pi
        magerr = float(np.max(np.abs(magnitude-expected_mag)))
        phaseerr = float(np.max(np.abs(phase-expected_phase)))
        worst_magnitude = max(worst_magnitude,magerr)
        worst_phase = max(worst_phase,phaseerr)
        if magerr > 5e-8 or phaseerr > 1e-8:
            raise AssertionError((number,mode,ns,'Bode error',magerr,phaseerr))
        samples += len(w)
        channels += no*ni*len(w)
    if list(output):
        raise AssertionError('extra output in Fortran probe')
    print(f'PASS: {len(cases)} MIMO systems, {samples} frequencies, '
          f'{channels} channel evaluations; maximum scaled error '
          f'{worst_response:.3e}; max dB error {worst_magnitude:.3e}; '
          f'max phase error {worst_phase:.3e} degrees')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--driver',type=Path)
    args = parser.parse_args()
    if args.driver:
        check(args.driver.resolve())
        return
    with tempfile.TemporaryDirectory(prefix='scifort-mimo-') as temp:
        work = Path(temp)
        objs = []
        for name in MODULES:
            src = ROOT/'src'/f'{name}.f90'
            obj = work/f'{name}.o'
            subprocess.run(['gfortran','-O2','-std=f2018',
                            '-Wimplicit-interface','-Werror=implicit-interface',
                            '-J',str(work),'-I',str(work),'-c',str(src),
                            '-o',str(obj)],check=True,cwd=work)
            objs.append(str(obj))
        exe = work/'lti_mimo_driver'
        subprocess.run(['gfortran','-O2','-J',str(work),'-I',str(work),
                        str(ROOT/'tools/signal_reference/lti_mimo_driver.f90'),
                        *objs,'-o',str(exe)],check=True,cwd=work)
        check(exe)


if __name__ == '__main__':
    main()
