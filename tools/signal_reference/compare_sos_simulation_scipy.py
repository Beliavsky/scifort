#!/usr/bin/env python3
"""Independent SciPy 1.17 comparison for arbitrary-input SOS simulation.

SciFort's dlsim_sos is a native extension, not a public scipy.signal routine.
SciPy's SOS filtering provides coefficient, sample, and DF-II state references.
"""
from __future__ import annotations

import argparse
import itertools
import subprocess

import numpy as np
from scipy import signal


def construct_cases() -> list[dict]:
    rng = np.random.default_rng(20261028)
    cases = []
    for variant, nsec, nin, interpolation, use_state in itertools.product(
            range(4), (1, 2, 5, 13, 64), (9, 31, 73), (0, 1), (0, 1)):
        dt = 0.04 + 0.0125 * variant
        sos = np.zeros((nsec, 6), dtype=np.complex128)
        for i in range(nsec):
            a1 = -0.40 + 0.13 * np.cos(0.45 * (i+1))
            a2 = 0.16 + 0.035 * np.sin(0.13 * (i+1))
            gain = 1 + a1 + a2
            sos[i] = [0.42*gain, 0.38*gain, 0.2*gain, 1, a1, a2]
            if variant in (2, 3):
                sos[i, 0] += 0.02j * gain
                sos[i, 1] -= 0.02j * gain
                sos[i, 4] += 0.005j * np.cos(i+1)
        xx = np.arange(nin, dtype=np.float64)
        u = (0.5*np.sin(0.21*xx) + np.cos(0.48*xx)
             + rng.normal(scale=0.2, size=nin)).astype(complex)
        if variant in (1, 3):
            u += 1j*(np.sin(0.39*xx) + rng.normal(scale=0.1, size=nin))
        t = (dt * xx * (0.7 + 0.006 * xx) if interpolation else None)
        zi = ((rng.normal(size=(nsec, 2)) * 0.025).astype(complex)
              if use_state else np.zeros((nsec, 2), dtype=complex))
        if variant in (1, 2, 3) and use_state:
            zi += 1j * rng.normal(size=(nsec, 2)) * 0.012
        if variant in (0, 1):
            sos.imag[:] = 0
        if variant in (0, 2):
            u.imag[:] = 0
        if variant == 0:
            zi.imag[:] = 0
        cases.append(dict(variant=variant, sos=sos, u=u, dt=dt, t=t,
                          zi=zi, state_flag=use_state, time_flag=interpolation))
    return cases


def feed_data(cases: list[dict]) -> str:
    lines = [str(len(cases))]
    for c in cases:
        sos, u, t, zi = c['sos'], c['u'], c['t'], c['zi']
        lines.append(f"{c['variant']} {len(sos)} {len(u)} "
                     f"{c['time_flag']} {c['state_flag']} {c['dt']:.17g}")
        for v in sos.ravel(order='C'):
            lines.append(f'{v.real:.17g} {v.imag:.17g}')
        for v in u:
            lines.append(f'{v.real:.17g} {v.imag:.17g}')
        if t is not None:
            lines += [f'{v:.17g}' for v in t]
        if c['state_flag']:
            for v in zi.ravel(order='C'):
                lines.append(f'{v.real:.17g} {v.imag:.17g}')
    return '\n'.join(lines) + '\n'


def complex_line(v: str) -> complex:
    parts = v.split()
    assert len(parts) == 2, v
    return complex(float(parts[0]), float(parts[1]))


def compare(executable: str) -> None:
    cases = construct_cases()
    run = subprocess.run([executable], input=feed_data(cases), text=True,
                         capture_output=True, check=True)
    lines = iter(run.stdout.splitlines())
    max_error = 0.0
    total = 0
    total_state = 0
    for i, case in enumerate(cases, start=1):
        label, number, status = next(lines).split()
        assert (label, int(number), int(status)) == ('CASE', i, 0), (i, label, number, status)
        length_label, nt_str = next(lines).split()
        assert length_label == 'LEN', i
        nt = int(nt_str)
        got_t = np.array([float(next(lines)) for _ in range(nt)])
        got_y = np.array([complex_line(next(lines)) for _ in range(nt)])
        nsec = len(case['sos'])
        got_x = np.array([complex_line(next(lines)) for _ in range(nt*2*nsec)]).reshape(nt,2*nsec)
        got_zf = np.array([complex_line(next(lines)) for _ in range(2*nsec)]).reshape(nsec,2)
        if case['t'] is None:
            expected_t = np.arange(len(case['u'])) * case['dt']
            given_u = case['u']
        else:
            expected_t = np.linspace(0, case['t'][-1], int(np.floor(case['t'][-1]/case['dt']))+1)
            given_u = (np.interp(expected_t, case['t'], case['u'].real)
                       + 1j*np.interp(expected_t, case['t'], case['u'].imag))
        if case['variant'] in (0, 2):
            given_u = given_u.real
        sos = case['sos'] if case['variant'] in (2, 3) else case['sos'].real
        zi = case['zi'] if case['variant'] != 0 else case['zi'].real
        sos = np.ascontiguousarray(sos)
        zi = np.ascontiguousarray(zi)
        given_u = np.ascontiguousarray(given_u)
        expected_y, expected_zf = signal.sosfilt(sos, given_u, zi=zi)
        expected_x = np.empty_like(got_x)
        state = zi.copy()
        for k, val in enumerate(given_u):
            expected_x[k,:] = state.reshape(-1)
            _, state = signal.sosfilt(sos,np.array([val]),zi=state)
        for name, expected, got in (
                ('time',expected_t,got_t),('output',expected_y,got_y),
                ('pre-input states',expected_x,got_x),('final states',expected_zf,got_zf)):
            error = np.max(np.abs(np.asarray(expected)-got))
            max_error = max(max_error, error)
            if not np.allclose(got,expected,rtol=8e-11,atol=8e-12):
                raise AssertionError(f'Case {i} {name}: max error {error:.4g} '
                                     f'variant={case["variant"]} sections={nsec}')
        total += len(got_y)
        total_state += len(got_x.ravel())
    try:
        extra = next(lines)
    except StopIteration:
        extra = None
    assert extra is None, f'Unexpected driver output: {extra}'
    print(f'PASS {len(cases)} SOS simulation cases, {total} output samples, '
          f'{total_state} state samples; max abs error {max_error:.4e}',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    args = parser.parse_args()
    compare(args.executable)
