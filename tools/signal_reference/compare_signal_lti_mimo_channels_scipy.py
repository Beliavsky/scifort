#!/usr/bin/env python3
"""Cross-check selected prepared MIMO frequency responses with SciPy."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.linalg import solve
from run_signal_suite import MODULES, ROOT


def generate_cases():
    rng = np.random.default_rng(20261018)
    for mode in (1, 2):
        for complex_system in (False, True):
            for i in range(70):
                ns = int(rng.choice([0, 1, 2, 5, 12, 30]))
                ni = int(rng.integers(1, 7))
                no = int(rng.integers(1, 6))
                w = rng.uniform(-2.9, 2.9, 7) if mode == 2 else rng.uniform(-5, 5, 7)
                scale = 0.03
                if ns:
                    a = np.diag((-rng.uniform(.2, 4, ns) if mode == 1
                                 else rng.uniform(-.8, .8, ns))).astype(complex)
                    a += scale*rng.standard_normal((ns, ns))
                    b = rng.standard_normal((ns, ni)).astype(complex)
                    c = rng.standard_normal((no, ns)).astype(complex)
                    if complex_system:
                        a += .02j*rng.standard_normal((ns, ns))
                        b += 1j*rng.standard_normal((ns, ni))
                        c += 1j*rng.standard_normal((no, ns))
                else:
                    a = np.empty((0, 0), complex)
                    b = np.empty((0, ni), complex)
                    c = np.empty((no, 0), complex)
                d = rng.standard_normal((no, ni)).astype(complex)
                if complex_system:
                    d += 1j*rng.standard_normal((no, ni))
                inputs = rng.integers(0, ni, size=int(rng.integers(1, ni + 3)))
                outputs = rng.integers(0, no, size=int(rng.integers(1, no + 3)))
                yield (mode, a, b, c, d, w, inputs, outputs)


def input_text(cases):
    lines = []
    for mode, a, b, c, d, w, inputs, outputs in cases:
        ns, ni, no = len(a), b.shape[1], c.shape[0]
        lines.append(f'{mode} {ns} {ni} {no} {len(w)} {len(inputs)} {len(outputs)}')
        for matrix in (a, b, c, d):
            lines.extend(f'{v.real:.17g} {v.imag:.17g}'
                         for v in matrix.ravel(order='F'))
        lines.extend(f'{v:.17g}' for v in w)
        lines.extend(str(int(v) + 1) for v in inputs)
        lines.extend(str(int(v) + 1) for v in outputs)
    return '\n'.join(lines) + '\n'


def check(driver):
    cases = list(generate_cases())
    p = subprocess.run([str(driver)], input=input_text(cases), text=True,
                       capture_output=True, check=True, timeout=240)
    rows = iter(p.stdout.splitlines())
    worst = 0.0
    values = 0
    for case_number, (mode, a, b, c, d, w, inputs, outputs) in enumerate(cases):
        statuses = tuple(map(int, next(rows).split()))
        if statuses != (0, 0, 0, len(w)):
            raise AssertionError((case_number, statuses))
        nr = len(outputs)
        nc = len(inputs)
        data = np.array([list(map(float, next(rows).split()))
                         for _ in range(len(w)*nr*nc)])
        if data.shape != (len(w)*nr*nc, 4):
            raise AssertionError((case_number, 'bad response shape'))
        data = data.reshape((len(w), nc, nr, 4)).transpose(0, 2, 1, 3)
        actual = data[..., 0] + 1j*data[..., 1]
        direct = data[..., 2] + 1j*data[..., 3]
        expect = []
        for omega in w:
            z = np.exp(1j*omega) if mode == 2 else 1j*omega
            if len(a):
                x = solve(z*np.eye(len(a))-a,b,check_finite=False)
                h = c @ x + d
            else:
                h = d
            expect.append(h[np.ix_(outputs, inputs)])
        expect = np.array(expect)
        err = float(np.max(np.abs(actual-expect)/(1+np.abs(expect))))
        d_err = float(np.max(np.abs(direct-expect)/(1+np.abs(expect))))
        worst = max(worst, err, d_err)
        if max(err,d_err) > 4e-11:
            raise AssertionError((case_number,mode,len(a),err,d_err))
        values += actual.size
    if list(rows):
        raise AssertionError('unexpected trailing driver output')
    print(f'PASS: {len(cases)} selected MIMO systems, {values} channel-frequency '
          f'values; maximum scaled SciPy difference {worst:.3e}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--driver', type=Path)
    args = parser.parse_args()
    if args.driver:
        check(args.driver.resolve())
        return
    with tempfile.TemporaryDirectory(prefix='scifort-mimo-channels-') as tmp:
        work = Path(tmp)
        objects = []
        for module in MODULES:
            source = ROOT/'src'/f'{module}.f90'
            obj = work/f'{module}.o'
            subprocess.run(['gfortran', '-O2', '-std=f2018',
                            '-Wimplicit-interface', '-Werror=implicit-interface',
                            '-J', str(work), '-I', str(work), '-c', str(source),
                            '-o', str(obj)], check=True, cwd=work)
            objects.append(str(obj))
        exe = work/'lti_mimo_channels_driver'
        subprocess.run(['gfortran', '-O2', '-J', str(work), '-I', str(work),
                        str(ROOT/'tools/signal_reference/lti_mimo_channels_driver.f90'),
                        *objects, '-o', str(exe)], check=True, cwd=work)
        check(exe)


if __name__ == '__main__':
    main()
