#!/usr/bin/env python3
"""Randomized minimum-order analog/digital IIR comparisons to SciPy.

Compiles a self-contained signal driver with GNU Fortran. Run:
    python tools/signal_reference/verify_iir_order_bands.py
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
import warnings

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
MODULES = [
    'kinds', 'signal_sequences', 'signal_convolve', 'signal_convolve2d', 'signal_detrend',
    'signal_peaks', 'signal_savgol', 'signal_lfilter', 'signal_sos',
    'signal_filtfilt', 'signal_filtfilt_nd', 'signal_dpss', 'signal_windows_extra', 'signal_spectral', 'signal_fir_advanced', 'signal_spectral_complex',
    'signal_analog_design', 'signal_design', 'signal_analog', 'signal_zpk_conversion', 'signal_state_space', 'signal_cont2discrete', 'signal_discrete_simulation', 'signal_continuous_simulation', 'signal_partial_fraction',
    'signal_elliptic', 'signal_analog_sos', 'signal_iir_extra', 'signal_fftconvolve', 'signal_fourier_resample',
    'signal_hilbert', 'signal_waveforms', 'signal_periodic',
    'signal_irregular', 'signal_resample', 'signal_remez', 'signal_czt', 'signal',
]
RNG = np.random.default_rng(20261016)


def make_case(family, analog, mode):
    """Sample valid stop/pass edge arrangements and attenuation targets."""
    fs = 2.0 if analog else float(RNG.choice([2., 100., 1000.]))
    scale = float(10**RNG.uniform(0.0, 2.0)) if analog else fs / 2
    gp = float(RNG.uniform(0.5, 4.))
    gs = float(RNG.uniform(17., 58.))
    if mode == 1:  # low
        wp = np.array([scale * RNG.uniform(0.10, 0.3)])
        ws = wp * float(RNG.uniform(1.4, 3.2))
    elif mode == 2:  # high
        wp = np.array([scale * RNG.uniform(0.50, 0.78)])
        ws = wp * float(RNG.uniform(0.35, 0.75))
    elif mode == 3:  # bandstop
        wp = np.array([0.11, 0.83])*scale
        ws = np.array([float(RNG.uniform(.22, .36)),
                       float(RNG.uniform(.51, .72))])*scale
    else:  # bandpass
        wp = np.array([float(RNG.uniform(.24, .36)),
                       float(RNG.uniform(.52, .68))])*scale
        ws = np.array([float(RNG.uniform(.08, .19)),
                       float(RNG.uniform(.75, .9))])*scale
    return gp, gs, fs, wp, ws


def main():
    cases = []
    for family in range(1, 4):
        for analog in (False, True):
            for mode in range(1, 5):
                for _ in range(45):
                    gp, gs, fs, wp, ws = make_case(family, analog, mode)
                    func = [signal.buttord, signal.cheb1ord, signal.cheb2ord][family-1]
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        n, cut = func(wp if mode > 2 else float(wp[0]),
                                      ws if mode > 2 else float(ws[0]),
                                      gp, gs, analog=analog,
                                      **({} if analog else {'fs':fs}))
                    if 1 <= n <= 10 and np.all(np.isfinite(cut)):
                        cases.append((family, analog, mode, gp, gs, fs, wp, ws,
                                      n, np.atleast_1d(cut)))
    rows = []
    for family, analog, mode, gp, gs, fs, wp, ws, *_ in cases:
        wp0 = float(wp[0]); wp1 = float(wp[1]) if mode > 2 else 0.0
        ws0 = float(ws[0]); ws1 = float(ws[1]) if mode > 2 else 0.0
        rows.append(f'{family} {int(analog)} {mode} {gp:.17g} {gs:.17g} '
                    f'{fs:.17g} {wp0:.17g} {wp1:.17g} {ws0:.17g} {ws1:.17g}')
    worst_cutoff = 0.
    worst_coeff = 0.
    designs_checked = 0
    with tempfile.TemporaryDirectory(prefix='scifort-order-bands-') as temp:
        work = Path(temp)
        flags = ['-std=f2018', '-Wimplicit-interface',
                 '-Werror=implicit-interface', '-O0', '-fcheck=all']
        objects = []
        for name in MODULES:
            obj = work / f'{name}.o'
            subprocess.run(['gfortran', *flags, '-I', str(work), '-J', str(work),
                            '-c', str(ROOT / 'src' / f'{name}.f90'),
                            '-o', str(obj)], cwd=work, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            objects.append(str(obj))
        exe = work / 'order_band_driver'
        subprocess.run(['gfortran', *flags, '-I', str(work),
                        str(ROOT / 'tools/signal_reference/order_band_driver.f90'),
                        *objects, '-o', str(exe)], check=True, cwd=work)
        proc = subprocess.run([str(exe)], input='\n'.join(rows)+'\n', text=True,
                              capture_output=True, check=True)
        lineiter = iter(proc.stdout.splitlines())
        for case_idx, (family, analog, mode, gp, gs, fs, wp, ws, nref, cref) in enumerate(cases):
            record = next(lineiter).split()
            status, order, nf = map(int, record[:3])
            measured = np.array([float(s) for s in record[3:3+nf]])
            dstatus, nb, na = map(int, record[5:8])
            if status != 0 or order != nref or nf != cref.size:
                raise AssertionError(f'case {case_idx}: order mismatch; family={family}, mode={mode}, '
                                     f'analog={analog}, wp={wp}, ws={ws}, result={record}')
            delta = max(abs(measured-cref))
            worst_cutoff = max(worst_cutoff, float(delta))
            scale = float(max(np.max(abs(cref)), 1e-6))
            if mode == 3:
                # SciPy's fminbound has a default 1e-5 absolute tolerance on
                # the analog-warped passband optimization interval.
                tolerance = 4e-5 * max(1.0, scale)
            else:
                tolerance = 2e-10 * scale
            if delta > tolerance:
                raise AssertionError(f'case {case_idx}: cutoff {measured} vs {cref}, '
                                     f'diff {delta}, tolerance {tolerance}; input {wp}, {ws}, '
                                     f'analog={analog}, fam={family}, mode={mode}')
            if dstatus != 0 or nb == 0 or na == 0:
                raise AssertionError(f'case {case_idx}: invalid design {record}')
            actual_b = np.fromstring(next(lineiter), sep=' ')
            actual_a = np.fromstring(next(lineiter), sep=' ')
            assert actual_b.size == nb and actual_a.size == na
            # Compare BA coefficients in well-conditioned low-order problems.
            if nref <= 5:
                family_name = ['butter', 'cheby1', 'cheby2'][family-1]
                kwargs = {'analog': analog, 'ftype': family_name, 'output': 'ba'}
                if not analog:
                    kwargs['fs'] = fs
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    br, ar = signal.iirdesign(wp if mode > 2 else float(wp[0]),
                                             ws if mode > 2 else float(ws[0]),
                                             gp, gs, **kwargs)
                # SciPy bandstop order optimization produces slightly
                # different edges under its default bounded-search tolerance.
                if mode != 3:
                    discrepancy = max(np.max(abs(actual_b-br)) / max(1, np.max(abs(br))),
                                      np.max(abs(actual_a-ar)) / max(1, np.max(abs(ar))))
                    worst_coeff = max(worst_coeff, float(discrepancy))
                    if discrepancy > 5e-7:
                        raise AssertionError(f'case {case_idx}: design BA coeff err {discrepancy}')
                    designs_checked += 1
        if next(lineiter,None) is not None:
            raise AssertionError('Unexpected trailing driver output')
    print(f'PASS: {len(cases)} order/cutoff cases, including analog/digital and all four bands')
    print(f'PASS: {designs_checked} low-order iirdesign coefficient comparisons')
    print(f'Max cutoff absolute discrepancy: {worst_cutoff:.3e} (bandstop uses optimizer)')
    print(f'Max scaled BA coefficient discrepancy: {worst_coeff:.3e}')

if __name__ == '__main__':
    main()
