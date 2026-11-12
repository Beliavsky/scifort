#!/usr/bin/env python3
"""Verify SciFort YT/KNV0 method controls and output against SciPy 1.17.

Run at SciFort root. Requires gfortran, NumPy and SciPy. Pole-placement
feedback gains are not unique, so the tests compare closed-loop poles and
iteration/determinant properties, not particular values of K.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import sys

import numpy as np
from scipy.signal import place_poles

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).parent))
from verify_mimo_place_poles_scipy import test_cases  # noqa: E402


def build(work):
    exe = work / 'place_poles_methods'
    sources = ('kinds', 'linalg', 'signal_place_poles_mimo',
               'signal_place_poles')
    cmd = ['gfortran', '-std=f2018', '-O2', '-fcheck=all',
           '-J', str(work), '-I', str(work)]
    cmd.extend(str(ROOT / 'src' / f'{name}.f90') for name in sources)
    cmd += [str(ROOT / 'tools/signal_reference/place_poles_methods_driver.f90'),
            '-o', str(exe)]
    subprocess.run(cmd, cwd=work, check=True, capture_output=True)
    return exe


def main():
    cases = test_cases()[:200]
    records = []
    for a, b, poles in cases:
        methods = ['YT']
        if np.all(np.abs(np.imag(poles)) < 1e-12):
            methods.append('KNV0')
        for method in methods:
            # Explicit negative rtol forces all iterations, as in SciPy.
            records.append((a, b, poles, method, 1))
            records.append((a, b, poles, method, 6))
    payload = []
    for a,b,poles,method,maxiter in records:
        payload.append(f'{a.shape[0]} {b.shape[1]} {method} {maxiter} -1')
        payload.append(' '.join(f'{v:.17g}' for v in a.ravel(order='F')))
        payload.append(' '.join(f'{v:.17g}' for v in b.ravel(order='F')))
        payload.append(' '.join(f'({v.real:.17g},{v.imag:.17g})' for v in poles))
    with TemporaryDirectory(prefix='scifort_pole_iterative_') as tmp:
        exe = build(Path(tmp))
        output = subprocess.run([str(exe)], input='\n'.join(payload)+'\n',
                                text=True, capture_output=True, check=True).stdout
    lines = output.splitlines()
    if len(lines) != len(records):
        raise AssertionError(f'{len(lines)} responses, expected {len(records)}')
    worst = 0.0
    determinants = []
    for i, ((a,b,poles,method,maxiter), line) in enumerate(zip(records,lines)):
        vals = np.fromstring(line,sep=' ')
        n,m = b.shape
        if vals.size != 4+n*m or int(vals[0]) != 0:
            raise AssertionError(f'case {i}: unexpected status or size: {vals[:5]}')
        count, rtol, determinant = int(vals[1]), vals[2], vals[3]
        if count != maxiter:
            raise AssertionError(f'case {i}: iteration count {count}, expected {maxiter}')
        if not np.isfinite(determinant) or determinant <= 0 or not 0 <= rtol <= 1:
            raise AssertionError(f'case {i}: invalid diagnostics {vals[:5]}')
        gain = vals[4:].reshape((m,n),order='F')
        got = np.linalg.eigvals(a-b@gain)
        error = max(abs(v-got[np.argmin(abs(got-v))]) for v in poles)
        worst = max(worst,error)
        if error > 3.e-6:
            raise AssertionError(f'case {i}: pole eigenvalue error {error}')
        if maxiter == 1:
            determinants.append(determinant)
        else:
            if determinant < determinants[-1]*(1-1.e-8):
                raise AssertionError(f'case {i}: determinant fell in later iteration')
            # SciPy is an independent feasibility and API reference.
            fit = place_poles(a,b,poles,method=method,rtol=-1,maxiter=6)
            assert fit.gain_matrix.shape == gain.shape
    print(f'Iterative YT/KNV0: {len(records)}/{len(records)} comparisons, '
          f'maximum closed-loop pole error {worst:.4g}')


if __name__ == '__main__':
    main()
