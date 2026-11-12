#!/usr/bin/env python3
"""Compare independently composed Fortran analog frequency transforms with SciPy."""
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from numpy.testing import assert_allclose
from scipy import signal
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
METHODS = ('lp2lp', 'lp2hp', 'lp2bp', 'lp2bs')


def test_cases():
    rng = np.random.default_rng(10429)
    cases = []
    for n in range(1, 7):
        for degree_b in range(n + 1):
            for typ in (0, 1):
                b = rng.uniform(-1, 1, degree_b + 1) + 1
                a = rng.uniform(-1, 1, n + 1) + 1
                if typ:
                    b = b + 1j * rng.uniform(-0.2, 0.2, len(b))
                    a = a + 1j * rng.uniform(-0.2, 0.2, len(a))
                for method in range(1, 5):
                    wo = 0.7 if n % 2 else 2.3
                    bw = 0.35 if degree_b % 2 else 1.4
                    cases.append((1, method, b, a, wo, bw, None))
    for n in range(1, 7):
        for degree_b in range(n + 1):
            # Real analog filters have real roots and/or conjugate pairs.
            # SciPy deliberately discards the imaginary part of the ZPK gain
            # ratio, so non-conjugate root sets are compared only against the
            # exact mathematical transfer function in the native tests.
            def paired_roots(count):
                roots = []
                if count % 2:
                    roots.append(complex(-0.4-rng.random(), 0.0))
                for _ in range(count//2):
                    root = complex(-0.4-rng.random(), 0.2+rng.random())
                    roots.extend((root, root.conjugate()))
                return np.array(roots, dtype=complex)
            z = paired_roots(degree_b)
            p = paired_roots(n)
            k = 1.7+0.31j if n % 2 else 1.7+0j
            for method in range(1, 5):
                cases.append((2, method, z, p, 0.75 + 0.1*n,
                              0.4 + 0.03*degree_b, k))
    # Also check equal-degree transforms with a zero lowpass numerator.
    for method in range(1, 5):
        cases.append((1, method, np.array([0.0]), np.array([1.0, 3.0, 2.0]), 1.3, 0.7, None))
    return cases


def encode(cases):
    lines = [str(len(cases))]
    for form, method, x, y, wo, bw, k in cases:
        lines.append(f'{form} {method} {len(x)} {len(y)} {wo:.17g} {bw:.17g}')
        for array in (x, y):
            for v in array:
                lines.append(f'{complex(v).real:.17g} {complex(v).imag:.17g}')
        if form == 2:
            lines.append(f'{k.real:.17g} {k.imag:.17g}')
    return '\n'.join(lines)+'\n'


def nearest_root_comparison(actual, expected):
    if len(actual) != len(expected):
        raise AssertionError(f'root counts differ: {len(actual)}, {len(expected)}')
    if len(actual) == 0:
        return 0.0
    costs = abs(actual[:, None] - expected[None, :])
    row, col = linear_sum_assignment(costs)
    assert_allclose(actual[row], expected[col], rtol=2e-11, atol=5e-12)
    return float(np.max(costs[row, col]))


def run():
    cases = test_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-transform-') as tmp:
        tmp = Path(tmp)
        compiler = ['gfortran', '-std=f2018', '-Werror=implicit-interface',
                    '-O0', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                    '-I', str(tmp), '-J', str(tmp)]
        objects = []
        for src in ('kinds', 'signal_frequency_transforms'):
            obj = tmp / (src + '.o')
            subprocess.run([*compiler, '-c', str(ROOT/'src'/(src+'.f90')),
                            '-o', str(obj)], cwd=tmp, check=True)
            objects.append(obj)
        exe = tmp/'analog_transform_driver'
        subprocess.run([*compiler, str(ROOT/'tools/signal_reference' /
                             'analog_frequency_transform_driver.f90'),
                        *map(str, objects), '-o', str(exe)], cwd=tmp, check=True)
        result = subprocess.run([str(exe)], input=encode(cases), text=True,
                                capture_output=True, check=True, cwd=tmp)
    lines = iter(result.stdout.splitlines())
    nvalue, max_ba, max_root = 0, 0.0, 0.0
    for idx, (form, method, x, y, wo, bw, k) in enumerate(cases, 1):
        header = next(lines).split()
        if header[:2] != ['CASE', str(idx)] or int(header[2]) != 0:
            raise AssertionError(f'invalid driver header: {header}')
        first, second = int(header[3]), int(header[4])
        def read_complex(n):
            return np.array([complex(*map(float, next(lines).split())) for _ in range(n)], dtype=complex)
        r1 = read_complex(first)
        r2 = read_complex(second)
        name = METHODS[method - 1] + ('_zpk' if form == 2 else '')
        call = getattr(signal, name)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', signal.BadCoefficients)
            if method < 3:
                expected = call(x, y, k, wo=wo) if form == 2 else call(x, y, wo=wo)
            else:
                expected = call(x, y, k, wo=wo, bw=bw) if form == 2 else call(x, y, wo=wo, bw=bw)
        if form == 1:
            e1, e2 = expected
            assert_allclose(r1, e1, rtol=5e-10, atol=2e-10,
                            err_msg=f'{name} numerator case {idx}')
            assert_allclose(r2, e2, rtol=5e-10, atol=2e-10,
                            err_msg=f'{name} denominator case {idx}')
            max_ba = max(max_ba, max(abs(r1 - e1), default=0.0),
                         max(abs(r2 - e2), default=0.0))
            nvalue += len(e1) + len(e2)
        else:
            gain = read_complex(1)[0]
            e1, e2, egain = expected
            max_root = max(max_root, nearest_root_comparison(r1,e1),
                           nearest_root_comparison(r2,e2))
            assert_allclose(gain, egain, rtol=1e-11, atol=1e-11,
                            err_msg=f'{name} gain case {idx}')
            nvalue += len(e1) + len(e2) + 1
    try:
        next(lines)
        raise AssertionError('unexpected driver output')
    except StopIteration:
        pass
    print(f'PASS: {len(cases)} SciPy analog transform cases, {nvalue} complex values; '
          f'max BA diff {max_ba:.3g}, max root diff {max_root:.3g}')


if __name__ == '__main__':
    run()
