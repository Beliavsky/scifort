#!/usr/bin/env python3
"""Compare independent Fortran Remez exchange designs with SciPy 1.17+.

Uses an equiripple-response error check as well as tap comparisons. Different
exchange grids can produce slightly different, but equally valid, minimax taps.
Requires NumPy, SciPy, and gfortran; no FPM or external Fortran libraries.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal
import scipy

ROOT = Path(__file__).resolve().parents[2]


def make_cases() -> list[tuple[int, np.ndarray, np.ndarray, np.ndarray, str, float]]:
    rng = np.random.default_rng(20261017)
    specs = []
    for i in range(240):
        n = int(rng.integers(9, 95))
        rate = float(rng.choice([1., 2., 48., 800.]))
        kind = ['bandpass', 'hilbert', 'differentiator'][i % 3]
        if kind == 'bandpass':
            shape = i % 4
            if shape == 0:
                p = rng.uniform(.18, .48)
                s = p + rng.uniform(.09, .18)
                edges = [0., p, s, 1.]
                desired = [1., 0.]
            elif shape == 1:
                n += (n+1) % 2  # highpass, Nyquist must not vanish
                s = rng.uniform(.18, .48)
                p = s + rng.uniform(.09, .18)
                edges = [0., s, p, 1.]
                desired = [0., 1.]
            elif shape == 2:
                edges = [0., .19, .29, .60, .70, 1.]
                desired = [0., 1., 0.]
            else:
                edges = [0., .17, .27, .45, .55, .75, .85, 1.]
                desired = [1., 0., .6, 0.]
            b = np.array(edges)*rate/2
        else:
            low = rng.uniform(.075, .18)
            high = rng.uniform(.83, .96)
            b = np.array([low, high])*rate/2
            desired = [1.]
        d = np.array(desired)*rng.uniform(.75, 1.2)
        w = rng.uniform(.8, 2.5, len(d))
        specs.append((n, b, d, w, kind, rate))
    return specs


def measure_error(h: np.ndarray, n: int, bands: np.ndarray, desired: np.ndarray,
                  weights: np.ndarray, kind: str, rate: float) -> float:
    errors = []
    for j, (lo, hi) in enumerate(bands.reshape(-1, 2)):
        frequencies = np.linspace(lo, hi, 512)
        w = 2*np.pi*frequencies/rate
        _, complex_h = signal.freqz(h, worN=w)
        amp = complex_h * np.exp(1j * w * (n-1)/2)
        amp = amp.real if kind == 'bandpass' else amp.imag
        d = desired[j] if kind != 'differentiator' else desired[j]*frequencies/rate
        weight = weights[j]
        if kind == 'differentiator' and desired[j] != 0.0:
            weight = weight*rate/np.maximum(frequencies, 1e-14*rate)
        errors.append(np.max(np.abs(amp-d)*weight))
    return float(max(errors))


def main() -> None:
    cases = make_cases()
    lines = [str(len(cases))]
    for n, bands, desired, weights, kind, rate in cases:
        codes = {'bandpass': 1, 'hilbert': 2, 'differentiator': 3}
        lines.append(f'{n} {len(desired)} {codes[kind]} 16 {rate:.17g}')
        for arr in (bands, desired, weights):
            lines.append(' '.join(f'{v:.17g}' for v in arr))
    with tempfile.TemporaryDirectory(prefix='scifort-remez-') as tmpdir:
        tmp = Path(tmpdir)
        cmd = ['gfortran', '-std=f2018', '-O2', '-Werror=implicit-interface', '-fcheck=all',
               '-J', str(tmp), '-I', str(tmp)]
        cmd += [str(ROOT / 'src' / f'{name}.f90') for name in ('kinds', 'signal_remez')]
        cmd += [str(ROOT / 'tools/signal_reference/remez_driver.f90'), '-o', str(tmp/'remez')]
        subprocess.run(cmd, cwd=tmp, check=True)
        process = subprocess.run([str(tmp/'remez')], input='\n'.join(lines)+'\n', text=True,
                                 capture_output=True, cwd=tmp, check=True, timeout=180)
    results = iter(process.stdout.splitlines())
    checked = 0
    nonconverged = 0
    reference_failed = 0
    worst_coeff_scaled = 0.0
    worst_ratio = 0.0
    max_excess = 0.0
    for i, (n, bands, desired, weights, kind, rate) in enumerate(cases):
        st, it, size, ripple = next(results).split()
        st, it, size = int(st), int(it), int(size)
        if st:
            if st != 2 or size != 0:
                raise AssertionError(f'case {i}: unexpected status={st}, size={size}')
            nonconverged += 1
            continue
        actual = np.fromstring(next(results), sep=' ')
        if size != n or len(actual) != n or not np.all(np.isfinite(actual)):
            raise AssertionError(f'case {i}: invalid returned coefficient array')
        try:
            expected = signal.remez(n, bands, desired, weight=weights,
                                    type=kind, fs=rate, maxiter=30)
        except ValueError:
            # Even SciPy can fail to converge on highly demanding filters.
            reference_failed += 1
            continue
        diff = float(np.max(np.abs(expected-actual)))
        aerr = measure_error(actual, n, bands, desired, weights, kind, rate)
        eerr = measure_error(expected, n, bands, desired, weights, kind, rate)
        ratio = aerr / max(eerr, 1e-12)
        if ratio > 1.035 and aerr > eerr + 1e-6:
            raise AssertionError(f'case {i}: weighted ripple {aerr:.6g} '
                                 f'exceeds SciPy {eerr:.6g} (ratio={ratio:.3f})')
        if diff > .03*(1+np.max(np.abs(expected))):
            raise AssertionError(f'case {i}: max tap diff={diff:.4g}')
        worst_coeff_scaled = max(worst_coeff_scaled, diff/(1+np.max(np.abs(expected))))
        if eerr >= 1e-6:
            worst_ratio = max(worst_ratio, ratio)
        max_excess = max(max_excess, aerr-eerr)
        checked += 1
    if next(results, None) is not None:
        raise AssertionError('Unexpected extra driver output')
    if checked < 230:
        raise AssertionError(f'Only {checked} comparisons passed; {nonconverged} nonconverged')
    print(f'PASS: {len(cases) - nonconverged}/{len(cases)} SciFort equiripple designs, '
          f'{checked} compared with SciPy {scipy.__version__}')
    print(f'  {reference_failed} SciPy reference nonconvergence, '
          f'{nonconverged} SciFort nonconvergence statuses')
    print(f'  maximum scale-normalized tap difference: {worst_coeff_scaled:.3e}')
    print(f'  maximum meaningful weighted-ripple ratio vs SciPy: {worst_ratio:.5f}')
    print(f'  maximum absolute weighted-ripple excess vs SciPy: {max_excess:.3e}')


if __name__ == '__main__':
    main()
