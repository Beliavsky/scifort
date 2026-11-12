#!/usr/bin/env python3
"""Compare parallel SOS bank simulation against independent SciPy sosfilt calls.

SciFort's multichannel dlsim_sos_bank is an extension rather than a public
SciPy function. SciPy's sosfilt supplies the numerical and DF-II references.
"""
from __future__ import annotations

import argparse
import itertools
import subprocess

import numpy as np
from scipy import signal


def make_cases() -> list[dict]:
    rng = np.random.default_rng(20261028)
    cases = []
    for variant, nsec, nch, time_flag, state_flag in itertools.product(
        range(8), (1, 3, 9, 24), (1, 3), (0, 1), (0, 1)
    ):
        nt = 13 + (nsec % 3) * 5 + (2 * nch) + (variant % 4)
        dt = 0.05 + 0.006 * (variant % 4)
        sos = np.zeros((nsec, 6, nch), dtype=complex)
        for channel in range(nch):
            for section in range(nsec):
                a1 = -0.35 + 0.06 * np.cos(0.3 * (section + 1))
                a2 = 0.09 + 0.03 * np.sin(0.7 * (section + 1))
                gain = 1 + a1 + a2
                scale = 0.95 + 0.035 * channel if variant >= 4 else 0.95
                sos[section, :, channel] = [
                    0.45 * gain * scale, 0.35 * gain * scale,
                    0.20 * gain * scale, 1, a1, a2,
                ]
                if variant % 4 in (2, 3):
                    sos[section, 0, channel] += 1j * gain * 0.004
                    sos[section, 2, channel] -= 1j * gain * 0.004
                    sos[section, 4, channel] += 1j * 0.003
        xx = np.arange(nt, dtype=float)[:, None]
        yy = np.arange(nch, dtype=float)[None, :]
        u = (np.sin(xx * 0.14 + yy * 0.37)
             + 0.4 * np.cos(xx * 0.39 - yy)
             + rng.normal(0, 0.05, size=(nt, nch))).astype(complex)
        if variant % 4 in (1, 3):
            u += 1j * (0.2 * np.cos(xx * 0.18 - yy)
                       + rng.normal(0, 0.02, size=(nt, nch)))
        t = dt * np.arange(nt) * (0.8 + 0.003 * np.arange(nt)) if time_flag else None
        zi = (rng.normal(0, 0.008, (nsec, 2, nch))
              + 1j * rng.normal(0, 0.005, (nsec, 2, nch))).astype(complex)
        if variant % 4 == 0:
            zi = zi.real.astype(complex)
        if not state_flag:
            zi[:] = 0
        cases.append(dict(variant=variant, sos=sos, u=u, t=t, zi=zi,
                          dt=dt, time_flag=time_flag, state_flag=state_flag))
    return cases


def feed(cases: list[dict]) -> str:
    lines = [str(len(cases))]
    for case in cases:
        variant = case['variant']
        sos, u, zi, t = case['sos'], case['u'], case['zi'], case['t']
        nsec, _, nch = sos.shape
        nin = len(u)
        lines.append(f"{variant} {nsec} {nin} {nch} "
                     f"{case['time_flag']} {case['state_flag']} {case['dt']:.17g}")
        for k in range(nch):
            for section in range(nsec):
                for coef in range(6):
                    z = sos[section, coef, k]
                    lines.append(f"{z.real:.17g} {z.imag:.17g}")
        for k in range(nch):
            for i in range(nin):
                z = u[i, k]
                lines.append(f"{z.real:.17g} {z.imag:.17g}")
        if t is not None:
            lines.extend(f"{v:.17g}" for v in t)
        if case['state_flag']:
            for k in range(nch):
                for section in range(nsec):
                    for q in range(2):
                        z = zi[section, q, k]
                        lines.append(f"{z.real:.17g} {z.imag:.17g}")
    return '\n'.join(lines) + '\n'


def complex_line(value: str) -> complex:
    vals = value.split()
    assert len(vals) == 2, value
    return complex(float(vals[0]), float(vals[1]))


def compare(executable: str) -> None:
    cases = make_cases()
    run = subprocess.run([executable], input=feed(cases), text=True,
                         capture_output=True, check=True)
    lines = iter(run.stdout.splitlines())
    total_output = 0
    total_state = 0
    max_error = 0.0
    for num, case in enumerate(cases, start=1):
        label, row, status = next(lines).split()
        assert (label, int(row), int(status)) == ('CASE', num, 0), (num, label, status)
        length_label, nt_str = next(lines).split()
        assert length_label == 'LEN', (num, length_label)
        nt = int(nt_str)
        nsec, _, nch = case['sos'].shape
        got_t = np.array([float(next(lines)) for _ in range(nt)])
        got_y = np.empty((nt, nch), dtype=complex)
        got_x = np.empty((nt, 2 * nsec, nch), dtype=complex)
        got_zf = np.empty((nsec, 2, nch), dtype=complex)
        for k in range(nch):
            got_y[:, k] = [complex_line(next(lines)) for _ in range(nt)]
        for k in range(nch):
            for i in range(nt):
                got_x[i, :, k] = [complex_line(next(lines)) for _ in range(2 * nsec)]
        for k in range(nch):
            for section in range(nsec):
                got_zf[section, :, k] = [complex_line(next(lines)) for _ in range(2)]
        if case['t'] is None:
            expected_t = np.arange(len(case['u'])) * case['dt']
            inputs = case['u']
        else:
            expected_t = np.linspace(0, case['t'][-1],
                                     int(np.floor(case['t'][-1]/case['dt'])) + 1)
            inputs = np.column_stack([
                np.interp(expected_t, case['t'], case['u'][:, j].real)
                + 1j * np.interp(expected_t, case['t'], case['u'][:, j].imag)
                for j in range(nch)
            ])
        assert len(expected_t) == nt
        expected_y = np.empty_like(got_y)
        expected_x = np.empty_like(got_x)
        expected_zf = np.empty_like(got_zf)
        for channel in range(nch):
            sos = case['sos'][:, :, 0 if case['variant'] < 4 else channel]
            if case['variant'] % 4 in (0, 1):
                sos = sos.real
            u = inputs[:, channel]
            if case['variant'] % 4 in (0, 2):
                u = u.real
            zi = case['zi'][:, :, channel]
            if case['variant'] % 4 == 0:
                zi = zi.real
            sos = np.ascontiguousarray(sos)
            zi = np.ascontiguousarray(zi)
            u = np.ascontiguousarray(u)
            expected_y[:, channel], expected_zf[:, :, channel] = signal.sosfilt(sos, u, zi=zi)
            state = zi.copy()
            for i, sample in enumerate(u):
                expected_x[i, :, channel] = state.ravel(order='C')
                _, state = signal.sosfilt(sos, np.array([sample]), zi=state)
        for name, expected, got in (
            ('time', expected_t, got_t), ('output', expected_y, got_y),
            ('pre-input states', expected_x, got_x), ('final states', expected_zf, got_zf)
        ):
            err = np.max(np.abs(expected - got))
            max_error = max(max_error, err)
            if not np.allclose(expected, got, rtol=8e-11, atol=8e-12):
                raise AssertionError(f"case {num} {name}: error {err:.3e}, "
                                     f"variant={case['variant']} sections={nsec}")
        total_output += got_y.size
        total_state += got_x.size
    try:
        extra = next(lines)
    except StopIteration:
        extra = None
    assert extra is None, f"unexpected driver line: {extra}"
    print(f"PASS {len(cases)} parallel SOS cases; {total_output} outputs, "
          f"{total_state} state samples; max abs error {max_error:.4e}", flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    args = parser.parse_args()
    compare(args.executable)
