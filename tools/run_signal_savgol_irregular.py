#!/usr/bin/env python3
"""Compare SciFort's irregular-grid local-polynomial filter with NumPy WLS.

Also check uniform odd-window special cases against scipy.signal.savgol_filter.
From project root, run: python tools/run_signal_savgol_irregular.py
"""
from __future__ import annotations

import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.signal import savgol_filter

ROOT = Path(__file__).resolve().parents[1]


def make_cases():
    """Deterministic real, complex, weighted, derivative and edge cases."""
    rng = np.random.default_rng(120197)
    cases = []
    for case in range(220):
        n = int(rng.integers(14, 65))
        if case % 4 == 0:
            x = np.arange(n, dtype=float) * float(rng.uniform(0.10, 0.95))
        else:
            x = np.cumsum(rng.uniform(0.10, 0.60, n))
        wl = int(rng.integers(5, min(n, 15) + 1))
        order = int(rng.integers(1, min(wl, 6)))
        deriv = int(rng.integers(0, order + 2))
        complex_signal = case % 2 == 0
        weighted = case % 3 != 0
        real_signal = rng.standard_normal(n) + np.sin(x * 0.47)
        imag_signal = rng.standard_normal(n) - np.cos(x * 0.71) if complex_signal else np.zeros(n)
        weights = np.exp(rng.uniform(-0.8, 0.8, n))
        if case % 7 == 0 and weighted and wl > order + 1:
            weights[int(rng.integers(0, n))] = 0.0
        query = np.sort(rng.uniform(x[0] - 0.02, x[-1] + 0.02, int(rng.integers(5, 27))))
        if case % 5 == 0:
            query[::3] = x[np.linspace(0, n-1, len(query[::3]), dtype=int)]
        cases.append((x, real_signal+1j*imag_signal, weights, query, wl, order, deriv,
                      complex_signal, weighted))
    return cases


def expected_for(case):
    x, y, weights, query, wl, order, deriv, complex_signal, weighted = case
    if not weighted:
        weights = np.ones_like(x)
    result = []
    for at in query:
        # Select the closest consecutive observations; a stable tie selects
        # the earlier timestamp (the rule used by the Fortran implementation).
        candidates = np.argsort(np.abs(x-at),kind='stable')[:wl]
        indices = np.sort(candidates)
        assert indices[-1] - indices[0] + 1 == wl
        time = x[indices]
        scale = max(abs(time-at))
        basis = np.column_stack([((time-at)/scale)**j for j in range(order+1)])
        root_weight = np.sqrt(weights[indices] / max(weights[indices]))
        fitted = np.linalg.lstsq(root_weight[:,None]*basis,
                                 root_weight*y[indices],rcond=None)[0]
        value = math.factorial(deriv)*fitted[deriv]/scale**deriv if deriv <= order else 0j
        result.append(value)
    return np.asarray(result)


def input_stream(cases):
    lines = [str(len(cases))]
    for x, y, w, q, wl, order, deriv, complex_signal, weighted in cases:
        lines.append(f'{len(x)} {len(q)} {wl} {order} {deriv} {int(complex_signal)} {int(weighted)}')
        for row in (x,y.real,y.imag,w,q):
            lines.append(' '.join(f'{float(v):.17g}' for v in row))
    return '\n'.join(lines)+'\n'


def compare(exe, cases, reference):
    output = subprocess.run([str(exe)],input=input_stream(cases),text=True,
                            capture_output=True,check=True).stdout.splitlines()
    reader = iter(output)
    worst = 0.0
    compared = 0
    for k,(case, expected) in enumerate(zip(cases,reference)):
        status, count = (int(value) for value in next(reader).split())
        assert (status,count) == (0,len(expected)), (k,status,count,len(expected))
        got = np.asarray([complex(*map(float,next(reader).split())) for _ in expected])
        np.testing.assert_allclose(got,expected,rtol=8e-9,atol=3e-8,
                                   err_msg=f'case {k}: window={case[4]}, order={case[5]}, derivative={case[6]}')
        worst = max(worst,float(np.max(np.abs(got-expected))))
        compared += len(got)
    assert next(reader,None) is None
    return compared,worst


def check_uniform_scipy(exe):
    rng = np.random.default_rng(511)
    cases=[]
    for k in range(60):
        n = int(rng.integers(20,60))
        delta = float(rng.uniform(0.15,1.3))
        x = delta*np.arange(n)
        signal = rng.normal(size=n)
        wl = int(rng.choice([5,7,9,11,13]))
        order = int(rng.integers(1,min(wl,5)))
        deriv = int(rng.integers(0,order+1))
        cases.append((x,signal.astype(complex),np.ones(n),x,wl,order,deriv,False,False))
    output = subprocess.run([str(exe)],input=input_stream(cases),text=True,
                            capture_output=True,check=True).stdout.splitlines()
    reader = iter(output)
    worst = 0.0
    for i,case in enumerate(cases):
        x,y,_,_,wl,order,deriv,_,_ = case
        assert tuple(map(int,next(reader).split())) == (0,len(x))
        got=np.asarray([float(next(reader).split()[0]) for _ in x])
        ref=savgol_filter(y.real,wl,order,deriv=deriv,delta=x[1]-x[0],mode='interp')
        np.testing.assert_allclose(got,ref,rtol=2e-9,atol=2e-9,err_msg=f'uniform test {i}')
        worst=max(worst,float(np.max(np.abs(got-ref))))
    assert next(reader,None) is None
    return len(cases), worst


def build_and_run(directory, optimized, cases, reference):
    directory.mkdir()
    flags = ['-O3','-march=native'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    common = ['gfortran','-std=f2018',*flags,'-Wall','-Wextra',
              '-Werror=implicit-interface','-J',str(directory),'-I',str(directory)]
    sources = [ROOT/'src/kinds.f90',ROOT/'src/signal_savgol.f90']
    native = directory/'test_signal_savgol_irregular'
    subprocess.run([*common,*(str(f) for f in sources),
                    str(ROOT/'test/test_signal_savgol_irregular.f90'),'-o',str(native)],
                   cwd=directory,check=True,capture_output=True,text=True)
    subprocess.run([str(native)],check=True,capture_output=True,text=True)
    reference_program = directory/'reference_signal_savgol_irregular'
    subprocess.run([*common,*(str(f) for f in sources),
                    str(ROOT/'tools/reference_signal_savgol_irregular.f90'),'-o',str(reference_program)],
                   cwd=directory,check=True,capture_output=True,text=True)
    compared,worst = compare(reference_program,cases,reference)
    n_uniform,uniform_error = check_uniform_scipy(reference_program)
    print(f"  {len(cases)} NumPy WLS comparisons / {compared} outputs; maximum error {worst:.3e}")
    print(f"  {n_uniform} uniform-grid SciPy comparisons; maximum error {uniform_error:.3e}")
    print('  Native Fortran regression passed')


def main():
    cases = make_cases()
    reference = [expected_for(case) for case in cases]
    with tempfile.TemporaryDirectory(prefix='scifort_savgol_irregular_') as scratch:
        for optimized in (False, True):
            print('optimized -O3:' if optimized else 'checked -fcheck=all:')
            build_and_run(Path(scratch)/('optimized' if optimized else 'checked'),
                          optimized,cases,reference)
    print('All irregular Savitzky-Golay checks passed.')


if __name__ == '__main__':
    main()
