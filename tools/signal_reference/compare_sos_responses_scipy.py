#!/usr/bin/env python3
"""Numerically validate SOS-native impulse/step responses against SciPy 1.17.

Unlike discrete-time TF response routines this test never converts a large
SOS cascade to high-degree polynomial coefficients. It uses scipy.signal.sosfilt
as the independent time-domain reference and scipy.signal.dimpulse/dstep for
low-order transfer-function agreement.
"""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
import numpy as np
from scipy import signal


@dataclass
class Case:
    sos: np.ndarray
    n: int
    dt: float
    impulse: bool
    zi: np.ndarray | None
    stride: int


def cases():
    rng = np.random.default_rng(20261028)
    result = []
    for i in range(288):
        complex_filter = i % 3 == 0
        nsections = int(rng.integers(1, 8))
        if i % 13 == 0:
            nsections = int(rng.choice([12, 24, 40, 64]))
            order = 2 * nsections
            cutoff = float(rng.uniform(.09, .68))
            kind = 'lowpass' if i % 2 else 'highpass'
            sos = signal.butter(order, cutoff, output='sos', btype=kind)
        elif complex_filter:
            sos = np.empty((nsections,6), dtype=np.complex128)
            for j in range(nsections):
                poles = rng.uniform(.05,.8,2) * np.exp(1j*rng.uniform(-3,3,2))
                zeros = rng.uniform(.05,.6,2) * np.exp(1j*rng.uniform(-3,3,2))
                b = np.poly(zeros) * rng.uniform(.25,.7)
                a = np.poly(poles)
                sos[j] = np.r_[b,a]
        else:
            # Conjugate-pole biquads, strictly within unit disk.
            sos = np.empty((nsections,6), dtype=float)
            for j in range(nsections):
                p = rng.uniform(.05,.75)*np.exp(1j*rng.uniform(.05,2.7))
                z = rng.uniform(.05,.7)*np.exp(1j*rng.uniform(.05,2.7))
                sos[j] = np.r_[rng.uniform(.35,.85)*np.poly([z,z.conjugate()]),
                               np.poly([p,p.conjugate()])].real
        complex_filter = bool(np.iscomplexobj(sos))
        n = int(rng.choice([1,2,3,7,15,39,80,121,257]))
        dt = float(rng.choice([.01, .1, .2, .5, 1., 2.]))
        zi = None
        if i % 4 == 0:
            zi = rng.normal(0., .015, (nsections,2))
            if complex_filter:
                zi = zi + 1j*rng.normal(0., .015, (nsections,2))
        result.append(Case(sos, n, dt, bool(i%2), zi, (2 if i % 5 == 0 else 1)))
    return result


def encode(test_cases: list[Case]) -> str:
    lines = [str(len(test_cases))]
    for c in test_cases:
        s = c.sos
        lines.append(f'{len(s)} {c.n} {c.dt:.17g} {int(np.iscomplexobj(s))} '
                     f'{int(c.impulse)} {int(c.zi is not None)} {c.stride}')
        for sec in s:
            lines.append(' '.join(f'{v:.17g}' for v in [*sec.real,*sec.imag]))
        if c.zi is not None:
            for state in c.zi:
                lines.append(' '.join(f'{v:.17g}' for v in [*state.real,*state.imag]))
    return '\n'.join(lines)+'\n'


def compare(executable: str) -> None:
    test_cases = cases()
    completed = subprocess.run([executable], input=encode(test_cases), text=True,
                               capture_output=True, check=True)
    rows = iter(completed.stdout.splitlines())
    worst = 0.
    total = 0
    for case_id, c in enumerate(test_cases):
        status, n = map(int,next(rows).split())
        if status != 0 or n != c.n:
            raise AssertionError(f'case {case_id}: invalid result status={status}, n={n}')
        got = np.array([list(map(float,next(rows).split())) for _ in range(n)])
        x = np.zeros(n, dtype=c.sos.dtype) if c.impulse else np.ones(n, dtype=c.sos.dtype)
        if c.impulse:
            x[0] = 1
        expected = signal.sosfilt(c.sos,x,zi=c.zi)[0] if c.zi is not None else signal.sosfilt(c.sos,x)
        error = np.max(np.abs(got[:,1]+1j*got[:,2]-expected))
        # Condition-aware tolerance on long cascades and moderate-polynomial cases.
        bound = 7.e-12 * max(1.,float(np.max(np.abs(expected))))
        if not (np.isfinite(error) and error <= bound):
            raise AssertionError(f'case {case_id}: max_error={error:g} > {bound:g}')
        clock = np.arange(n)*c.dt
        if not np.allclose(got[:,0],clock,atol=3.e-13,rtol=5.e-14):
            raise AssertionError(f'case {case_id}: time grid mismatch')
        # Compare low-order zero-state models with SciPy's public LTI response.
        if len(c.sos) <= 3 and c.zi is None and not np.iscomplexobj(c.sos) and n > 1:
            b, a = signal.sos2tf(c.sos)
            if c.impulse:
                tt, response = signal.dimpulse((b,a,c.dt),n=n)
            else:
                tt, response = signal.dstep((b,a,c.dt),n=n)
            if not np.allclose(np.ravel(response[0]),expected,atol=1.e-9,rtol=1.e-9):
                raise AssertionError(f'case {case_id}: SciPy LTI reference mismatch')
        worst = max(worst,float(error))
        total += n
    if next(rows,None) is not None:
        raise AssertionError('unexpected rows in Fortran output')
    print(f'PASS {len(test_cases)} SOS response cases, {total} samples, '
          f'max abs error {worst:.3e}; SciPy {signal.__version__ if hasattr(signal,"__version__") else "1.17"}',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', help='compiled sos_response_driver path')
    args = parser.parse_args()
    compare(args.executable)
