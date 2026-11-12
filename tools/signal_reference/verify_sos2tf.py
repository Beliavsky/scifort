#!/usr/bin/env python3
"""Compare the real second-order-section conversion with scipy.signal.sos2tf.

Run from any directory with gfortran, NumPy and SciPy installed:
    python tools/signal_reference/verify_sos2tf.py
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
RNG = np.random.default_rng(20261015)


def cases():
    yield np.empty((0, 6), dtype=np.float64)
    for _ in range(400):
        count = int(RNG.integers(1, 11))
        sos = RNG.normal(0, 0.4, size=(count, 6))
        sos[:, 0] = RNG.uniform(-1.5, 1.5, size=count)
        sos[:, 3] = RNG.uniform(0.3, 1.8, size=count)
        if count and RNG.random() < .3:
            sos[-1, 2] = 0.0
            sos[-1, 5] = 0.0
        yield sos
    for family in (signal.butter, signal.cheby1, signal.cheby2):
        for order in range(1, 13):
            for mode in ('lowpass', 'highpass', 'bandpass', 'bandstop'):
                cutoff = .37 if mode in ('lowpass', 'highpass') else [.23, .61]
                args = ((order, 1.0, cutoff) if family == signal.cheby1
                        else (order, 40.0, cutoff) if family == signal.cheby2
                        else (order, cutoff))
                yield family(*args, btype=mode, output='sos')


def main() -> None:
    inputs = []
    expected = []
    for sos in cases():
        expected.append(signal.sos2tf(sos))
        inputs.append(str(len(sos)))
        for section in sos:
            inputs.append(' '.join(f'{v:.17g}' for v in section))

    with tempfile.TemporaryDirectory(prefix='scifort-sos2tf-') as temporary:
        path = Path(temporary)
        exe = path / 'sos2tf_driver'
        command = [
            'gfortran', '-std=f2018', '-O2', '-Werror=implicit-interface',
            '-J', str(path), '-I', str(path),
            str(ROOT / 'src/kinds.f90'),
            str(ROOT / 'src/signal_spectral.f90'),
            str(ROOT / 'src/signal_analog_design.f90'),
            str(ROOT / 'src/signal_design.f90'),
            str(ROOT / 'tools/signal_reference/sos2tf_driver.f90'),
            '-o', str(exe),
        ]
        subprocess.run(command, cwd=path, check=True, capture_output=True)
        response = subprocess.run([str(exe)], input='\n'.join(inputs) + '\n',
                                  text=True, check=True, capture_output=True)

    lines = iter(response.stdout.splitlines())
    max_scaled_error = 0.0
    for index, (reference_b, reference_a) in enumerate(expected):
        status, length = map(int, next(lines).split())
        actual_b = np.fromstring(next(lines), sep=' ')
        actual_a = np.fromstring(next(lines), sep=' ')
        if status != 0 or length != len(reference_b):
            raise AssertionError(f'case {index}: status or coefficient length mismatch')
        if len(actual_b) != length or len(actual_a) != length:
            raise AssertionError(f'case {index}: malformed coefficient output')
        scale = max(1.0, float(np.max(np.abs(reference_b))),
                    float(np.max(np.abs(reference_a))))
        error = max(float(np.max(np.abs(actual_b - reference_b))),
                    float(np.max(np.abs(actual_a - reference_a)))) / scale
        max_scaled_error = max(max_scaled_error, error)
        if error > 5e-13:
            raise AssertionError(f'case {index}: relative-scale difference {error:.3g}')
    print(f'sos2tf: {len(expected)} SciPy 1.17 comparison cases passed; '
          f'worst scaled coefficient error = {max_scaled_error:.3g}')


if __name__ == '__main__':
    main()
