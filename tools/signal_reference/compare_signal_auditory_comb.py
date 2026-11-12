#!/usr/bin/env python3
"""Compare SciFort IIR comb and gammatone coefficients to SciPy 1.17.

Usage: python tools/signal_reference/compare_signal_auditory_comb.py
Requires Python numpy/scipy and GNU Fortran; uses an isolated temporary build.
"""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.signal import iircomb, gammatone

ROOT = Path(__file__).resolve().parents[2]
MODULES = [
    'kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral',
    'signal_analog_design', 'signal_elliptic', 'signal_design', 'signal_iir_extra',
]


def main() -> None:
    rng = np.random.default_rng(20261021)
    commands = []
    reference = []
    for _ in range(450):
        fs = float(rng.choice([2, 200, 1000, 16000, 48000]))
        n = int(rng.choice([3, 4, 5, 7, 8, 10, 20, 33, 80, 125, 250]))
        f = fs / n
        q = float(rng.uniform(2, 250))
        family = str(rng.choice(['notch', 'peak']))
        dc = bool(rng.integers(0, 2))
        commands.append(f'c {family} {f:.17g} {q:.17g} {fs:.17g} 4 16 {int(dc)}')
        reference.append(iircomb(f, q, ftype=family, fs=fs, pass_zero=dc))
    for _ in range(450):
        fs = float(rng.choice([800, 1000, 1600, 8000, 16000, 48000]))
        f = float(rng.uniform(0.002, 0.49) * fs)
        order = int(rng.integers(1, 25))
        length = int(rng.integers(1, 181))
        commands.append(f'g fir {f:.17g} 0 {fs:.17g} {order} {length} 0')
        reference.append(gammatone(f, 'fir', order=order, numtaps=length, fs=fs))
    for _ in range(450):
        fs = float(rng.choice([800, 1000, 1600, 8000, 16000, 48000]))
        f = float(rng.uniform(0.002, 0.49) * fs)
        commands.append(f'g iir {f:.17g} 0 {fs:.17g} 4 16 0')
        reference.append(gammatone(f, 'iir', fs=fs))

    with tempfile.TemporaryDirectory(prefix='scifort-auditory-') as directory:
        work = Path(directory)
        objects = []
        for name in MODULES:
            obj = work / (name + '.o')
            subprocess.run([
                'gfortran', '-std=f2018', '-O2', '-Wimplicit-interface',
                '-Werror=implicit-interface', '-I', str(work), '-J', str(work),
                '-c', str(ROOT/'src'/f'{name}.f90'), '-o', str(obj),
            ], cwd=work, check=True)
            objects.append(str(obj))
        executable = work / 'reference'
        subprocess.run([
            'gfortran', '-std=f2018', '-O2', '-I', str(work),
            str(ROOT/'tools/signal_reference/signal_auditory_comb_driver.f90'),
            *objects, '-o', str(executable),
        ], cwd=work, check=True)
        run = subprocess.run([str(executable)], input='\n'.join(commands)+'\n',
                             text=True, stdout=subprocess.PIPE, check=True, cwd=work)
    lines = run.stdout.splitlines()
    if len(lines) != 3*len(commands):
        raise AssertionError(f'Expected {3*len(commands)} lines; got {len(lines)}')
    worst_absolute = 0.0
    worst_relative = 0.0
    for i, (expected_b, expected_a) in enumerate(reference):
        status, nb, na = map(int, lines[3*i].split())
        if status != 0:
            raise AssertionError(f'case {i}: unexpected status {status} for {commands[i]}')
        b = np.fromstring(lines[3*i+1], dtype=float, sep=' ')
        a = np.fromstring(lines[3*i+2], dtype=float, sep=' ')
        if len(b) != nb or len(a) != na:
            raise AssertionError(f'case {i}: size mismatch for {commands[i]}')
        if b.shape != expected_b.shape or a.shape != expected_a.shape:
            raise AssertionError(f'case {i}: shape mismatch for {commands[i]}')
        for found, expected, label in [(b, expected_b, 'b'), (a, expected_a, 'a')]:
            if not np.allclose(found, expected, rtol=5e-10, atol=1e-12):
                delta = np.max(np.abs(found - expected))
                raise AssertionError(f'case {i} {label} mismatch {delta} for {commands[i]}')
            absolute = float(np.max(np.abs(found - expected)))
            worst_absolute = max(worst_absolute, absolute)
            norm = float(np.max(np.abs(expected)))
            if norm > 1e-18:
                worst_relative = max(worst_relative, absolute/norm)
    print(f'PASS: {len(commands)}/{len(commands)} SciPy filter comparisons; '
          f'maximum absolute error {worst_absolute:.3g}; '
          f'maximum relative coefficient error {worst_relative:.3g}')


if __name__ == '__main__':
    main()
