#!/usr/bin/env python3
"""Compare SciFort analog/digital partial fractions and round trips to SciPy.

Requires numpy, scipy and gfortran, but not fpm.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import numpy as np
import scipy
from scipy import signal
from run_signal_suite import MODULES, ROOT


def cases():
    rng = np.random.default_rng(1726)
    samples = []
    for digital in (False, True):
        for order in range(1, 7):
            for i in range(45):
                if i % 13 == 0 and order >= 2:
                    base = rng.uniform(0.2, 0.65)
                    pole_values = [base, base]
                    pole_values += list(rng.uniform(-0.8, 0.8, order-2))
                    poles = np.asarray(pole_values)
                elif i % 5 == 0 and order >= 2:
                    pole_values = [0.4+0.25j, 0.4-0.25j]
                    pole_values += list(rng.uniform(-0.8, 0.8, order-2))
                    poles = np.asarray(pole_values)
                else:
                    poles = rng.uniform(-0.85, 0.85, order)
                a = np.poly(poles)
                # Proper, equal-order, and improper polynomials.
                length = max(1, order + [0, 1, 2, -1][i % 4])
                if i % 9 == 0:
                    b = rng.normal(size=length) + 1j*rng.normal(size=length)
                else:
                    b = rng.normal(size=length)
                if i % 17 == 0:
                    b *= 0
                if i % 8 == 0:
                    a *= rng.uniform(0.2,2.0)
                samples.append((digital, np.asarray(b,complex), np.asarray(a,complex)))
    return samples


def render_stdin(samples):
    rows = [str(len(samples))]
    for digital, b, a in samples:
        rows.append(f'{2 if digital else 1} {len(b)} {len(a)}')
        for v in np.concatenate([b,a]):
            rows.append(f'{v.real:.17e} {v.imag:.17e}')
    return '\n'.join(rows)+'\n'


def parse_result(lines, pos):
    status, nr, nk, npoles = map(int, lines[pos].split())
    pos += 1
    if status:
        return (status, np.array([]), np.array([]), np.array([]), None), pos

    def consume(count):
        nonlocal pos
        out = np.empty(count,dtype=complex)
        for i in range(count):
            re,im = map(float,lines[pos].split())
            out[i] = re+1j*im
            pos += 1
        return out

    r = consume(nr)
    p = consume(npoles)
    k = consume(nk)
    tf_status, nb, na = map(int, lines[pos].split())
    pos += 1
    bout = consume(nb)
    aout = consume(na)
    return (status, r,p,k, (tf_status,bout,aout)), pos


def grouped(r,p):
    # Allow magnitude-tied conjugate roots to be emitted in either order.
    remaining = set(range(len(p)))
    groups = []
    while remaining:
        first = min(remaining)
        same = [j for j in sorted(remaining) if abs(p[j]-p[first]) < 1e-3]
        groups.append((np.mean(p[same]),r[same]))
        remaining.difference_update(same)
    return groups


def response(coeff_b, coeff_a, digital, x):
    if digital:
        return np.polynomial.polynomial.polyval(x,coeff_b)/np.polynomial.polynomial.polyval(x,coeff_a)
    return np.polyval(coeff_b,x)/np.polyval(coeff_a,x)


def check_case(case, got, index):
    digital,b,a = case
    operation = signal.residuez if digital else signal.residue
    expected_r, expected_p, expected_k = operation(b,a)
    status,r,p,k,tf = got
    if status:
        raise AssertionError(f'case {index}: Fortran status {status}')
    if len(p) != len(expected_p):
        raise AssertionError(f'case {index}: pole count {len(p)} != {len(expected_p)}')
    expected_groups = grouped(expected_r,expected_p)
    for gp,gr in grouped(r,p):
        candidates = [(abs(ep-gp), j,er) for j,(ep,er) in enumerate(expected_groups)]
        delta,j,er = min(candidates,key=lambda x:x[0])
        if delta > 3e-4 or len(gr) != len(er):
            raise AssertionError(f'case {index}: pole mismatch {gp} or group length {len(gr)}')
        if not np.allclose(gr,er,rtol=3e-3,atol=3e-4):
            raise AssertionError(f'case {index}: residue mismatch {gr} vs {er}')
        expected_groups.pop(j)
    if expected_groups:
        raise AssertionError(f'case {index}: missing groups')
    if not np.allclose(k,expected_k,atol=3e-7,rtol=3e-7):
        raise AssertionError(f'case {index}: direct {k} vs {expected_k}')
    tf_status,bout,aout = tf
    if tf_status:
        raise AssertionError(f'case {index}: inverse status {tf_status}')
    x = np.array([0.21+0.63j,0.87+0.19j,1.23+0.43j])
    if digital:
        x = np.exp(-1j*np.array([0.29,0.83,1.74]))
    yr = response(b,a,digital,x)
    yf = response(bout,aout,digital,x)
    if not np.allclose(yr,yf,rtol=2e-5,atol=2e-5):
        raise AssertionError(f'case {index}: transfer response mismatch {yr} vs {yf}')


def main():
    samples = cases()
    with tempfile.TemporaryDirectory(prefix='scifort-pfrac-') as directory:
        work = Path(directory)
        flags = ['-std=f2018','-O2','-Wimplicit-interface','-Werror=implicit-interface']
        objects=[]
        for name in MODULES:
            src = ROOT/'src'/f'{name}.f90'
            obj = work/f'{name}.o'
            subprocess.run(['gfortran',*flags,'-I',str(work),'-J',str(work),
                '-c',str(src),'-o',str(obj)],check=True,cwd=work)
            objects.append(str(obj))
        exe = work/'probe'
        subprocess.run(['gfortran',*flags,'-I',str(work),'-J',str(work),
            str(Path(__file__).with_name('probe_partial_fraction.f90')),*objects,
            '-o',str(exe)],check=True,cwd=work)
        process = subprocess.run([str(exe)],input=render_stdin(samples),text=True,
            capture_output=True,check=True)
        lines=process.stdout.splitlines()
        pos=0
        for i,case in enumerate(samples):
            result,pos=parse_result(lines,pos)
            check_case(case,result,i)
        if pos != len(lines):
            raise AssertionError('unconsumed probe output')
    print(f'PASS: {len(samples)} SciPy {scipy.__version__} partial-fraction comparisons and inverse round trips')


if __name__ == '__main__':
    main()
