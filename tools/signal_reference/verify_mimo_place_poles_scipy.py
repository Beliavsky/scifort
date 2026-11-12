#!/usr/bin/env python3
"""Verify SciFort's partially actuated MIMO place_poles against SciPy.

Run from the SciFort root. Requires gfortran, NumPy and SciPy. The gain
matrix is not unique for MIMO, so the reference is the requested closed-loop
spectrum, not SciPy's particular optimized gain matrix.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import warnings

import numpy as np
from scipy.signal import place_poles

ROOT = Path(__file__).resolve().parents[2]


def compile_driver(work):
    exe = work / 'place_poles_mimo_compare'
    subprocess.run([
        'gfortran', '-std=f2018', '-O2', '-fcheck=all', '-fbacktrace',
        '-J', str(work), '-I', str(work),
        str(ROOT / 'src/kinds.f90'),
        str(ROOT / 'src/linalg.f90'),
        str(ROOT / 'src/signal_place_poles_mimo.f90'),
        str(ROOT / 'src/signal_place_poles.f90'),
        str(ROOT / 'tools/signal_reference/place_poles_driver.f90'),
        '-o', str(exe)], check=True, capture_output=True, cwd=work)
    return exe


def test_cases():
    rng = np.random.default_rng(229011)
    cases = []
    for i in range(480):
        n = int(rng.integers(3, 11))
        m = int(rng.integers(2, n))
        a = rng.normal(size=(n, n))*.6
        b = rng.normal(size=(n, m))*.8
        if i % 7 == 0:
            b[:, 0] *= .01
        if i % 17 == 0:
            a *= 3
        if i % 41 == 0:
            b[:, 1] = b[:, 0]
        poles = -(0.3+rng.random(n)*4.0)
        if i % 3 == 0:
            poles = poles.astype(complex)
            poles[0] = -0.5-rng.random()+1j*rng.uniform(.25, 3)
            poles[1] = poles[0].conjugate()
        if i % 5 == 0:
            poles[-1] = poles[-2] if not (np.iscomplexobj(poles) and poles[-2].imag != 0) \
                else -(rng.random()+.3)
        if max(np.count_nonzero(abs(poles-v) < 1.e-10) for v in poles) > \
                np.linalg.matrix_rank(b):
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                place_poles(a, b, poles)
        except ValueError:
            continue
        cases.append((a, b, poles))

    # KNV example from SciPy's public documentation (real desired poles).
    a = np.array([[1.380, -.2077, 6.715, -5.676],
                  [-.5814, -4.290, 0., .675],
                  [1.067, 4.273, -6.654, 5.893],
                  [.0480, 4.273, 1.343, -2.104]])
    b = np.array([[0.,5.679],[1.136,1.136],[0.,0.],[-3.146,0.]])
    cases.append((a,b,np.array([-.2,-.5,-5.0566,-8.6659])))
    # The complex-pole example from the same documentation.
    a = np.array([[0., 7./3., 0., 0.], [0.,0.,0.,7./9.],
                  [0.,0.,0.,0.], [0.,0.,0.,0.]])
    b = np.array([[0.,0.],[0.,0.],[1.,0.],[0.,1.]])
    poles = np.array([-3, -1, -2-1j, -2+1j])/3.
    cases.append((a,b,poles))
    return cases


def main():
    cases = test_cases()
    payload = []
    for a,b,poles in cases:
        payload.append(f'{a.shape[0]} {b.shape[1]}')
        payload.append(' '.join(f'{v:.17g}' for v in a.ravel(order='F')))
        payload.append(' '.join(f'{v:.17g}' for v in b.ravel(order='F')))
        payload.append(' '.join(f'({v.real:.17g},{v.imag:.17g})' for v in poles))
    with TemporaryDirectory(prefix='scifort_mimo_poles_') as tmp:
        exe = compile_driver(Path(tmp))
        process = subprocess.run([str(exe)], input='\n'.join(payload)+'\n',
                                 text=True, capture_output=True, check=True)
    lines = process.stdout.splitlines()
    if len(lines) != len(cases):
        raise AssertionError(f'{len(lines)} lines; expected {len(cases)}')
    worst_error = 0.
    for i, (line, (a,b,poles)) in enumerate(zip(lines,cases)):
        result = np.fromstring(line,sep=' ')
        n,m = b.shape
        if len(result) != 1+n*m or result[0] != 0:
            raise AssertionError(f'case {i}: invalid status {result[:1]}')
        gain = result[1:].reshape((m,n),order='F')
        computed = np.linalg.eigvals(a-b@gain)
        error = np.max(np.abs(np.sort_complex(computed)-np.sort_complex(poles)))
        worst_error = max(worst_error,float(error))
        if error > 2.e-6:
            raise AssertionError(f'case {i}: eigenvalue error {error}')
    print(f'MIMO pole placement: {len(cases)}/{len(cases)} SciPy-validated cases; '
          f'maximum eigenvalue error {worst_error:.4g}')


if __name__ == '__main__':
    main()
