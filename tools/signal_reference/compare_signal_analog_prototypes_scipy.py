#!/usr/bin/env python3
"""Compare SciFort analog lowpass prototype ZPK functions with SciPy 1.17+."""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import signal
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]


def run_driver(exe: Path, family: str, n: int, arg: str | None):
    command = [str(exe), family, str(n)]
    if arg is not None:
        command.append(arg)
    lines = subprocess.check_output(command, text=True).splitlines()
    status, nz, npoles, gain = lines[0].split()
    if int(status) != 0:
        raise AssertionError(f'{family}({n}, {arg}) returned status {status}')
    nz, npoles = int(nz), int(npoles)

    def read_roots(rows):
        return np.asarray([complex(*map(float, row.split())) for row in rows], dtype=complex)

    z = read_roots(lines[1:1+nz])
    p = read_roots(lines[1+nz:1+nz+npoles])
    if len(lines) != 1+nz+npoles:
        raise AssertionError('unexpected data from Fortran driver')
    return z, p, float(gain)


def compare_roots(actual, expected, desc):
    if len(actual) != len(expected):
        raise AssertionError(f'{desc}: root counts differ')
    if not len(actual):
        return 0.0, 0
    costs = abs(actual[:, None] - expected[None, :])
    row, col = linear_sum_assignment(costs)
    differences = costs[row, col]
    # Higher-order Bessel polynomial poles are somewhat ill conditioned.
    bound = (3e-8 if desc.startswith('besselap') else 2e-12)
    if not np.all(differences <= bound*np.maximum(1.0, abs(expected[col]))):
        raise AssertionError(f'{desc}: root mismatch {max(differences)}')
    return float(max(differences)), len(actual)


def cases():
    for n in range(0, 25):
        yield 'buttap', n, None
        for rp in (0.01, 0.1, 1.0, 5.0, 30.0, 80.0):
            yield 'cheb1ap', n, str(rp)
        for rs in (0.01, 10.0, 20.0, 50.0, 90.0):
            yield 'cheb2ap', n, str(rs)
    for n in range(0, 21):
        for mode in ('phase', 'delay', 'mag'):
            yield 'besselap', n, mode


def main():
    checked = 0
    n_values = 0
    max_root_error = 0.0
    max_gain_relative_error = 0.0
    max_frequency_error = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-prototypes-') as temporary:
        work = Path(temporary)
        flags = ['-std=f2018', '-Werror=implicit-interface', '-O0',
                 '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                 '-I', str(work), '-J', str(work)]
        objects = []
        for name in ('kinds', 'constants', 'math', 'special/bessel_i', 'special/bessel_k',
                     'special/bessel_jy', 'special/bessel_order', 'special/bessel_derivatives',
                     'special/hankel', 'special/bessel_complex', 'signal_analog_prototypes'):
            object_path = work/(name.replace('/','_')+'.o')
            subprocess.run(['gfortran', *flags, '-c', str(ROOT/'src'/(name+'.f90')),
                            '-o',str(object_path)], cwd=work, check=True)
            objects.append(str(object_path))
        executable = work/'analog_prototype_driver'
        subprocess.run(['gfortran', *flags, str(ROOT/'tools/signal_reference'/
                        'analog_prototype_driver.f90'), *objects, '-o',str(executable)],
                       cwd=work, check=True)
        frequencies = np.asarray([0.0, 0.1, 0.5, 1.0, 2.0, 5.0], dtype=float)
        for family, n, arg in cases():
            actual = run_driver(executable, family, n, arg)
            scipy_function = getattr(signal, family)
            expected = scipy_function(n) if arg is None else scipy_function(n, arg if family == 'besselap' else float(arg))
            description = f'{family}({n}, {arg})'
            a_z, a_p, a_k = actual
            e_z, e_p, e_k = expected
            e_z = np.asarray(e_z,dtype=complex)
            e_p = np.asarray(e_p,dtype=complex)
            root_error, count = compare_roots(a_z,e_z,description)
            max_root_error = max(max_root_error,root_error)
            n_values += count
            root_error, count = compare_roots(a_p,e_p,description)
            max_root_error = max(max_root_error,root_error)
            n_values += count
            gain_relative_error = abs(a_k-e_k)/max(1.0,abs(e_k))
            max_gain_relative_error = max(max_gain_relative_error,gain_relative_error)
            if gain_relative_error > (1e-8 if family == 'besselap' else 5e-11):
                raise AssertionError(f'{description}: gain mismatch {gain_relative_error}')
            n_values += 1
            response = signal.freqs_zpk(*actual, worN=frequencies)[1]
            expected_response = signal.freqs_zpk(*expected, worN=frequencies)[1]
            difference = abs(response-expected_response)
            relative = difference / np.maximum(1.0,abs(expected_response))
            max_frequency_error = max(max_frequency_error,float(relative.max()))
            if np.any(relative > (2e-8 if family == 'besselap' else 3e-10)):
                raise AssertionError(f'{description}: response mismatch {relative.max()}')
            n_values += len(frequencies)
            checked += 1
    print(f'PASS: {checked} analog prototype cases, {n_values} root/gain/response values')
    print(f'maximum root error: {max_root_error:.3e}')
    print(f'maximum normalized gain error: {max_gain_relative_error:.3e}')
    print(f'maximum normalized frequency-response error: {max_frequency_error:.3e}')


if __name__ == '__main__':
    main()
