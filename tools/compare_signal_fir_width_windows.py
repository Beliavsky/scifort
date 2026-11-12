#!/usr/bin/env python3
"""Verify SciFort's FIR transition width and shaped windows against SciPy.

Run with the Fortran binary produced by tools/run_signal_fir_width_windows.py.
The independent reference uses scipy.signal.firwin and firwin2, with no
Fortran implementation source copied from SciPy.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
import scipy
from scipy import signal


def cases():
    fs_values = [2., 40., 480.]
    ratios = [0.004, 0.02, 0.10, 0.25]
    for mode in ('lowpass', 'bandpass', 'highpass', 'bandstop'):
        for i in (1, 2):
            n = 15 if i == 1 else 26
            if mode in ('highpass', 'bandstop') and n % 2 == 0:
                n += 1
            for fs in fs_values:
                nyq = fs / 2
                cutoff = ([0.28*nyq] if mode in ('lowpass', 'highpass')
                          else [0.25*nyq, 0.62*nyq])
                for ratio in ratios:
                    for scale in (False, True):
                        yield f'firwin {mode} n={n} fs={fs} width={ratio} scale={scale}', (
                            signal.firwin(n, cutoff, pass_zero=mode, fs=fs,
                                          width=ratio*nyq, window='ignored-name',
                                          scale=scale))
    gains_base = np.array([0.5, 0.90, 0.25, 0.60, 0.15])
    for typ in (1, 2, 3, 4):
        anti = typ in (3, 4)
        for i in (1, 2):
            n = (13 if i == 1 else 29) if typ in (1, 3) else (14 if i == 1 else 30)
            fs = 2. if i == 1 else 500.
            freq = np.array([0., 0.18, 0.5, 0.83, 1.])*fs/2
            gain = gains_base.copy()
            if anti:
                gain[0] = 0.
            if typ in (2, 3):
                gain[-1] = 0.
            for grid in (1, 2):
                nf = (33 if grid == 1 else 41) if n < 29 else (65 if grid == 1 else 73)
                for family in range(1, 9):
                    window = {
                        1: ('kaiser', 5.2),
                        2: ('gaussian', 3.7),
                        3: ('tukey', 0.37),
                        4: ('general_hamming', 0.62),
                        5: ('general_gaussian', 1.6, 4.4),
                        6: ('general_cosine', [0.4, 0.3, 0.2, 0.1]),
                        7: 'blackman',
                        8: None,
                    }[family]
                    if family == 6:
                        # SciPy 1.17 firwin2 forwards xp/device to get_window;
                        # the general_cosine family rejects those extra keywords.
                        # Compare the same mathematical operation with a directly
                        # constructed SciPy symmetric general-cosine taper.
                        expected = signal.firwin2(n, freq, gain, nfreqs=nf,
                                                  window=None, antisymmetric=anti, fs=fs)
                        expected *= signal.windows.general_cosine(n, [0.4,0.3,0.2,0.1], sym=True)
                    else:
                        expected = signal.firwin2(n, freq, gain, nfreqs=nf,
                                                  window=window, antisymmetric=anti, fs=fs)
                    if family == 8:
                        k = np.arange(n, dtype=float)
                        weights = 0.2 + 0.8*np.sin(0.5*np.pi*k/(n-1))**2
                        expected *= weights
                    yield f'firwin2 type={typ} n={n} grid={nf} win={family}', expected


def main(executable: str) -> None:
    out = subprocess.check_output([executable], text=True)
    tokens = iter(out.split())
    max_error = 0.
    count = 0
    coefficients = 0
    for label, expected in cases():
        try:
            st, n = int(next(tokens)), int(next(tokens))
            actual = np.array([float(next(tokens)) for _ in range(n)])
        except (StopIteration, ValueError) as exc:
            raise AssertionError(f'bad Fortran output for {label}') from exc
        if st != 0 or n != len(expected):
            raise AssertionError(f'{label}: status={st}, len={n} expected {len(expected)}')
        error = np.max(np.abs(actual-expected))
        tolerance = 2e-11*(1+np.max(np.abs(expected)))
        if not np.isfinite(actual).all() or error > tolerance:
            raise AssertionError(f'{label}: max error {error:.3e} > {tolerance:.3e}\n'
                                 f'SciFort={actual}\nSciPy={expected}')
        count += 1
        coefficients += n
        max_error = max(max_error, error)
    if next(tokens, None) is not None:
        raise AssertionError('Unexpected trailing data from Fortran reference driver')
    print(f'PASS: {count} FIR designs, {coefficients} coefficients compared '
          f'to SciPy {scipy.__version__}; max absolute error {max_error:.3e}')


if __name__ == '__main__':
    main(sys.argv[1])
