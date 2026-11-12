#!/usr/bin/env python3
"""Compare SciFort parameterized polyphase resampling against SciPy 1.17."""
from __future__ import annotations

import subprocess
import sys
import numpy as np
from scipy import signal
from scipy.signal import windows
from math import gcd


WINDOWS = {
    1: ('gaussian', 2.75),
    2: ('kaiser', 4.5),
    3: ('tukey', 0.65),
    4: ('general_gaussian', 1.3, 4.4),
    5: ('general_cosine', [0.42, 0.50, 0.08]),
    6: ('taylor', 4, 32),
    7: 'hann',
    8: ('kaiser', 5.0),
}
PAD = ['constant', 'mean', 'line', 'reflect', 'median']


def compare(exe: str) -> None:
    rng = np.random.default_rng(20261028)
    cases = []
    for kind in (0, 1):
        for family in WINDOWS:
            for num in range(16):
                up, down = [(2, 3), (3, 2), (5, 4), (7, 6)][num % 4]
                n = [9, 15, 27, 42][(num // 2) % 4]
                pad = (num + family) % (4 if kind else 5)
                x = rng.normal(size=n) + np.cos(np.arange(n) * 0.23) / 5
                if kind:
                    x = x + 1j * (rng.normal(size=n) + np.sin(np.arange(n) * .47))
                cases.append((kind, family, n, up, down, pad, x))
    payload = [str(len(cases))]
    for kind, family, n, up, down, pad, x in cases:
        payload.extend((f'{kind} {family} {n} {up} {down} {pad+1}',
                        ' '.join(f'{v:.17e}' for v in x.real),
                        ' '.join(f'{v:.17e}' for v in x.imag)))
    proc = subprocess.run([exe], input='\n'.join(payload) + '\n', capture_output=True,
                          text=True, check=True, timeout=70)
    lines = iter(proc.stdout.splitlines())
    values = 0
    maximum = 0.0
    for i, (kind, family, n, up, down, pad, x) in enumerate(cases):
        window = WINDOWS[family]
        if family == 5:
            # SciPy 1.17's firwin -> get_window forbids general_cosine
            # through its array-API xp/device options. Reconstruct its
            # documented symmetric taper using SciPy windows.general_cosine.
            rate = max(up // gcd(up, down), down // gcd(up, down))
            nw = 20 * rate + 1
            win = windows.general_cosine(nw, window[1], sym=True)
            taps = signal.firwin(nw, 1.0 / rate, window='boxcar', scale=False) * win
            taps /= taps.sum()
            window = taps
        expected = signal.resample_poly(x, up, down, window=window, padtype=PAD[pad])
        status, nout = map(int, next(lines).split())
        if status != 0 or nout != len(expected):
            raise AssertionError(f'case {i}: {kind=} {family=} status/length {(status,nout)} '
                                 f'expected {len(expected)}')
        actual = np.array([complex(*map(float, next(lines).split())) for _ in range(nout)])
        diffs = abs(actual - expected)
        delta = float(np.max(diffs))
        if delta > 1e-10:
            j = int(np.argmax(diffs))
            raise AssertionError(f'case {i}: {kind=} {family=} {up=} {down=} {PAD[pad]=} '
                                 f'{delta=}, sample {j}: {actual[j]} != {expected[j]}')
        maximum = max(maximum, delta)
        values += len(actual)
    print(f'PASS: {len(cases)} SciPy polyphase resampling cases / {values} samples; '
          f'maximum absolute error {maximum:.3e}', flush=True)


if __name__ == '__main__':
    compare(sys.argv[1])
