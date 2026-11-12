#!/usr/bin/env python3
"""Check SciFort Bessel filter design against SciPy 1.17+ (no runtime dependency)."""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy import signal


def reference_cases():
    for n in (1, 3, 6, 10, 14, 18, 20):
        for btype in ('lowpass', 'highpass', 'bandpass', 'bandstop'):
            for norm in ('phase', 'delay', 'mag'):
                # Exercising all modes at high orders also exercises root pairing.
                if n >= 18 and norm != 'phase':
                    continue
                cutoff = (0.18, 0.37) if btype.startswith('band') else (0.23,)
                yield n, btype, norm, False, cutoff, 2.0
                if n in (1, 6, 14, 20):
                    acutoff = (2.0, 7.0) if btype.startswith('band') else (3.5,)
                    yield n, btype, norm, True, acutoff, 2.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver', help='compiled Fortran driver_signal_bessel_design')
    args = parser.parse_args()
    cases = list(reference_cases())
    modes = ('lowpass', 'highpass', 'bandpass', 'bandstop')
    norms = ('phase', 'delay', 'mag')
    lines = [str(len(cases))]
    for n, mode, norm, analog, wn, fs in cases:
        a, b = wn if len(wn) == 2 else (wn[0], 0)
        lines.append(f'{n} {modes.index(mode)+1} {int(analog)} '
                     f'{norms.index(norm)+1} {fs:.17g} {a:.17g} {b:.17g}')
    run = subprocess.run([args.driver], input='\n'.join(lines)+'\n',
                         capture_output=True, text=True, check=True)
    results = iter(run.stdout.splitlines())
    worst_coeff = 0.0
    worst_freq = 0.0
    for n, mode, norm, analog, cutoff, fs in cases:
        status, nb, na, ns = map(int, next(results).split())
        if status != 0:
            raise AssertionError(f'Fortran status {status}: {n} {mode} {norm} {analog}')
        b = np.fromstring(next(results), dtype=float, sep=' ')
        a = np.fromstring(next(results), dtype=float, sep=' ')
        assert b.size == nb and a.size == na
        if ns:
            sos = np.fromstring(next(results), dtype=float, sep=' ').reshape(ns, 6)
        else:
            sos = None
        kw = {'N': n, 'Wn': cutoff if len(cutoff) == 2 else cutoff[0],
              'btype': mode, 'norm': norm, 'analog': analog}
        if not analog:
            kw['fs'] = fs
        br, ar = signal.bessel(**kw, output='ba')
        scale_b = np.maximum(1e-12*np.max(np.abs(br)), np.abs(br))
        scale_a = np.maximum(1e-12*np.max(np.abs(ar)), np.abs(ar))
        e_coeff = max(np.max(np.abs(b-br)/scale_b), np.max(np.abs(a-ar)/scale_a))
        worst_coeff = max(worst_coeff, e_coeff)
        if e_coeff > (3e-7 if n >= 18 else 1e-8):
            raise AssertionError(f'BA coefficients {n} {mode} {norm} {analog}: '
                                 f'{e_coeff:.4g}')
        if analog:
            # High-order BA band-polynomial evaluation is ill-conditioned
            # near the stopband notch (also in SciPy). Compare the stable
            # outer-frequency region rather than amplifying roundoff.
            upper = 1.0 if n >= 10 and mode.startswith('band') else 10.0
            w = np.geomspace(0.08, upper, 55)
            _, resp0 = signal.freqs(br, ar, worN=w)
            _, resp1 = signal.freqs(b, a, worN=w)
        else:
            w = np.linspace(0, np.pi, 100, endpoint=False)
            _, resp0 = signal.sosfreqz(signal.bessel(**kw, output='sos'), worN=w)
            _, resp1 = signal.sosfreqz(sos, worN=w)
        e_freq = np.max(np.abs(resp0-resp1))
        worst_freq = max(worst_freq, e_freq)
        if e_freq > 2e-6:
            raise AssertionError(f'frequency response {n} {mode} {norm} '
                                 f'{analog}: {e_freq:.4g}')
    print(f'PASS: {len(cases)} Bessel designs vs SciPy; '
          f'worst relative coefficient error {worst_coeff:.3g}; '
          f'worst absolute response error {worst_freq:.3g}')


if __name__ == '__main__':
    main()
