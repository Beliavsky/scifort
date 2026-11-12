#!/usr/bin/env python3
"""Independent SOS group-delay checks against SciPy and an analytic reference.

The SciPy 1.17 group_delay() API accepts transfer-function polynomials, not
second-order-section matrices. For low-order cascades we compare directly
against sos2tf -> group_delay; for high-order cascades we sum independent
per-section SciPy group_delay results and cross-check the analytic derivative.
"""
from __future__ import annotations

import argparse
import subprocess
import warnings
from pathlib import Path

import numpy as np
from scipy import signal


def make_cascade(rng: np.random.Generator, sections: int, complex_filter: bool):
    sos = np.empty((sections, 6), dtype=np.complex128)
    for j in range(sections):
        # Cover widely varying gains without ever forming their product.
        g = 10.0 ** rng.uniform(-95, 95)
        zr, pr = rng.uniform(0.1, 0.97, 2)
        theta_z, theta_p = rng.uniform(0.15, 2.97, 2)
        if complex_filter:
            z1 = zr * np.exp(1j * theta_z)
            z2 = zr * np.exp(1j * rng.uniform(-3.0, 3.0))
            p1 = pr * np.exp(1j * theta_p)
            p2 = pr * np.exp(1j * rng.uniform(-3.0, 3.0))
        else:
            z1 = zr * np.exp(1j * theta_z)
            z2 = np.conj(z1)
            p1 = pr * np.exp(1j * theta_p)
            p2 = np.conj(p1)
        sos[j] = [g, -g*(z1+z2), g*z1*z2, 1, -(p1+p2), p1*p2]
    if not complex_filter:
        return sos.real.astype(np.float64)
    return sos


def analytic(sos: np.ndarray, w: np.ndarray, fs: float) -> np.ndarray:
    q = np.exp(-2j*np.pi*w/fs)
    result = np.zeros(w.size, dtype=np.float64)
    for section in sos:
        b = section[:3] / np.max(abs(section[:3]))
        a = section[3:] / np.max(abs(section[3:]))
        numerator = b[0] + b[1]*q + b[2]*q*q
        denominator = a[0] + a[1]*q + a[2]*q*q
        result += np.real((b[1]*q+2*b[2]*q*q)/numerator -
                          (a[1]*q+2*a[2]*q*q)/denominator)
    return result


def scipy_reference(sos: np.ndarray, w: np.ndarray, fs: float) -> np.ndarray:
    # Each section's SciPy evaluation is independent, avoiding unstable
    # high-order polynomial convolution in the reference calculation.
    total = np.zeros(w.size)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        for section in sos:
            b = section[:3] / np.max(abs(section[:3]))
            a = section[3:] / np.max(abs(section[3:]))
            freq, gd = signal.group_delay((b,a),w=w,fs=fs)
            np.testing.assert_allclose(freq,w,rtol=1e-12,atol=1e-10)
            total += gd
    return total


def build_cases():
    rng = np.random.default_rng(20261028)
    cases = []
    sizes = (1,2,3,4,8,24,64,128)
    rates = (2*np.pi, 2.0, 48000.0)
    for is_complex in (False, True):
        for ns in sizes:
            for _ in range(7):
                rate = float(rates[len(cases) % len(rates)])
                # Avoid an exact zero/pole at a band endpoint.
                w = np.sort(rng.uniform(0.09, 0.91, 20))*rate
                sos = make_cascade(rng, ns, is_complex)
                cases.append((sos,w,rate))
    # Realistic high-order IIR designs with zeros exactly on the unit circle;
    # the tested frequencies deliberately exclude those frequencies.
    for order in (4, 8, 16, 32):
        for kind in ('lowpass','highpass','bandpass'):
            cutoff = [0.21,0.47] if kind == 'bandpass' else 0.32
            sos = signal.butter(order, cutoff, btype=kind, output='sos')
            rate = 2.0*np.pi
            w = np.linspace(0.047,0.94,31)*np.pi
            cases.append((sos,w,rate))
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    cases = build_cases()
    messages = [str(len(cases))]
    for sos, w, rate in cases:
        messages.append(f'{len(sos)} {len(w)} {int(np.iscomplexobj(sos))} {rate:.17e}')
        for row in sos:
            messages.append(' '.join(f'{component:.17e}' for z in row
                                     for component in (np.real(z),np.imag(z))))
        messages.append(' '.join(f'{v:.17e}' for v in w))
    result = subprocess.run([str(args.executable)],
                            input='\n'.join(messages)+'\n', text=True,
                            capture_output=True, check=True)
    lines = iter(result.stdout.splitlines())
    count_values = 0
    biggest_error = 0.0
    scipy_count = 0
    for ix,(sos,w,rate) in enumerate(cases):
        status_str,n_str = next(lines).split()
        assert int(status_str) == 0, (ix, 'status',status_str)
        assert int(n_str) == len(w), (ix, 'length',n_str)
        output = np.array([[float(x) for x in next(lines).split()]
                           for _ in range(len(w))])
        np.testing.assert_allclose(output[:,0],w,rtol=1e-15,atol=1e-12)
        predicted = analytic(sos,w,rate)
        # An absolute tolerance accounts for near-unit-circle pole/zero
        # cancellations and longer cascades with group delay O(100).
        np.testing.assert_allclose(output[:,1],predicted,rtol=4e-11,atol=2e-10)
        biggest_error = max(biggest_error,float(np.max(abs(output[:,1]-predicted))))
        if len(sos) <= 4 and ix < 112:
            normalized = np.array(sos, copy=True)
            for row in normalized:
                row[:3] /= np.max(abs(row[:3]))
            b,a = signal.sos2tf(normalized)
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                _, reference = signal.group_delay((b,a),w=w,fs=rate)
            np.testing.assert_allclose(output[:,1],reference,rtol=4e-10,atol=1e-9)
            scipy_count += 1
        else:
            independent = scipy_reference(sos,w,rate)
            np.testing.assert_allclose(output[:,1],independent,rtol=4e-11,atol=2e-10)
            scipy_count += 1
        count_values += len(w)
    assert next(lines,None) is None, 'extra driver output'
    print(f'PASS {len(cases)} independent SOS group-delay cases, '
          f'{count_values} delays; SciPy checked {scipy_count} cases; '
          f'max abs diff to analytic {biggest_error:.3e}')


if __name__ == '__main__':
    main()
