#!/usr/bin/env python3
"""Compare the SciFort advanced windows against SciPy 1.17 (including invalid cases)."""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
import numpy as np
import warnings
from scipy.signal import windows

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261021)
CASES = []

for n in [0, 1, 2, 3, 4, 5, 8, 15, 16, 23, 31, 32, 65, 96, 127, 128, 191]:
    for sy in [0, 1]:
        for att in [45., 60., 80., 100., 120., 150., 200.]:
            CASES.append((1, n, sy, att, 0., 1))
            CASES.append((4, n, sy, att, 0., 1))
        for nbar, sl in [(1,20.), (2,20.), (4,30.), (4,65.), (5,75.), (8,110.), (12,50.)]:
            for norm in [0, 1]:
                CASES.append((2, n, sy, float(nbar), sl, norm))
            CASES.append((5, n, sy, float(nbar), sl, 1))
        for beta in [0., 2., 8., 14., 30., 80.]:
            CASES.append((3, n, sy, beta, 0., 1))
            CASES.append((6, n, sy, beta, 0., 1))
for _ in range(120):
    n = int(RNG.integers(4, 180))
    sy = int(RNG.integers(0,2))
    fam = int(RNG.choice([1, 2, 3, 4, 5, 6]))
    if fam in [3, 6]: n = 2*(n//2)
    p1 = float(RNG.choice([0,3,5,7,9,11,17,20,35,45,60,80,100]))
    p2 = float(RNG.choice([20,30,40,55,75,110]))
    norm = int(RNG.integers(0,2))
    if fam in [2,5]: p1 = float(RNG.integers(1,13))
    CASES.append((fam,n,sy,p1,p2,norm))


def expected(case):
    f, n, sy, p1, p2, normal = case
    try:
        if f in [1,4]:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=UserWarning,
                    message='This window is not suitable for spectral analysis')
                return windows.chebwin(n, p1, sym=bool(sy))
        if f in [2,5]: return windows.taylor(n,int(p1),p2,norm=(bool(normal) if f==2 else True),sym=bool(sy))
        return windows.kaiser_bessel_derived(n,p1,sym=bool(sy))
    except ValueError: return None


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-window-check-') as temp:
        tmp = Path(temp)
        src = ROOT/'src'
        subprocess.run(['gfortran', '-O2', '-std=f2018', '-fcheck=bounds',
            '-J',str(tmp),'-I',str(tmp),str(src/'kinds.f90'),
            str(src/'signal_dpss.f90'),str(src/'signal_windows_extra.f90'),
            str(src/'signal_spectral.f90'),
            str(ROOT/'tools/signal_reference/window_advanced_driver.f90'),
            '-o',str(tmp/'driver')],cwd=tmp,check=True)
        inp = ''.join(' '.join(str(v) for v in case)+'\n' for case in CASES)
        output = subprocess.run([str(tmp/'driver')],cwd=tmp,input=inp,
            text=True,capture_output=True,check=True).stdout.splitlines()
    assert len(output) == len(CASES), (len(output), len(CASES))
    count, max_error = 0, 0.0
    for case,line in zip(CASES,output):
        toks = line.split()
        st, size = int(toks[0]),int(toks[1])
        actual = np.asarray([float(t.replace('D','E')) for t in toks[2:]])
        want = expected(case)
        if want is None:
            assert st==1 and size==0, (case,st,size)
        else:
            assert st==0 and size==len(want)==len(actual), (case,st,size)
            delta = np.max(abs(actual-want),initial=0.)
            max_error = max(max_error,delta)
            if not np.allclose(actual,want,atol=3e-10,rtol=3e-10):
                i = np.argmax(abs(actual-want))
                raise AssertionError((case,float(delta),i,float(actual[i]),float(want[i])))
        count += 1
    print(f'PASS: {count}/{len(CASES)} SciPy 1.17 advanced-window comparisons, max abs err {max_error:.3e}')


if __name__ == '__main__':
    main()
