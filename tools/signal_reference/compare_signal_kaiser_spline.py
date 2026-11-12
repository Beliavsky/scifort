#!/usr/bin/env python3
"""Independent comparison of SciFort B-spline/Kaiser helpers with SciPy 1.17.

Requires gfortran, scipy and numpy. Builds only the relevant Fortran modules.
"""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.signal import gauss_spline, kaiser_beta, kaiser_atten, kaiserord, firwin

ROOT = Path(__file__).resolve().parents[2]
MODULES = ('kinds', 'signal_dpss', 'signal_windows_extra', 'signal_spectral',
           'signal_analog_design', 'signal_design')


def main():
    rng = np.random.default_rng(20261021)
    commands = []
    expected = []
    for n in [0, 1, 2, 3, 5, 13, 100, 500, 10000]:
        for x in rng.uniform(-35, 35, 80):
            commands.append(f'g {x:.17g} 0 {n}')
            expected.append(('scalar', float(gauss_spline(x, n))))
    for a in np.r_[-10.0, 0., 8., 21., 21.001, 40., 50., 50.001,
                     rng.uniform(0, 190, 450)]:
        commands.append(f'b {a:.17g} 0 0')
        expected.append(('scalar', float(kaiser_beta(float(a)))))
    for n, width in zip(rng.integers(1,500,450), rng.uniform(0.0002,0.8,450)):
        commands.append(f'a {width:.17g} 0 {n}')
        expected.append(('scalar', float(kaiser_atten(int(n), float(width)))))
    for attenuation, width in zip(rng.uniform(8, 180, 500), rng.uniform(0.005, 0.9, 500)):
        if rng.uniform() < 0.2:
            attenuation = -attenuation
        commands.append(f'o {attenuation:.17g} {width:.17g} 0')
        expected.append(('order', kaiserord(float(attenuation), float(width))))
    for _ in range(260):
        n = int(rng.integers(3, 90))
        cutoff = float(rng.uniform(0.02, 0.94))
        beta = float(rng.uniform(0, 16))
        commands.append(f'f {cutoff:.17g} {beta:.17g} {n}')
        expected.append(('fir', firwin(n, cutoff, window=('kaiser',beta))))

    with tempfile.TemporaryDirectory(prefix='scifort-kaiser-spline-') as d:
        wd = Path(d)
        objs = []
        for name in MODULES:
            obj = wd/f'{name}.o'
            subprocess.run(['gfortran', '-std=f2018', '-O2', '-Wimplicit-interface',
                            '-Werror=implicit-interface', '-J', str(wd), '-I', str(wd),
                            '-c', str(ROOT/'src'/f'{name}.f90'), '-o', str(obj)],
                           cwd=wd, check=True, capture_output=True)
            objs.append(str(obj))
        executable = wd/'kaiser_reference'
        subprocess.run(['gfortran', '-std=f2018', '-O2', '-J', str(wd), '-I', str(wd),
                        str(ROOT/'tools/signal_reference/kaiser_spline_driver.f90'),
                        *objs, '-o', str(executable)], cwd=wd, check=True, capture_output=True)
        result = subprocess.run([str(executable)], input='\n'.join(commands)+'\n',
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=wd, check=True)
    lines = iter(result.stdout.splitlines())
    worst = 0.0
    for i, ((typ, ref), command) in enumerate(zip(expected, commands)):
        if typ == 'scalar':
            got = float(next(lines))
            worst = max(worst, abs(got-ref))
            if not np.isclose(got, ref, rtol=2e-13, atol=2e-13):
                raise AssertionError(f'case {i} {command}: got {got} versus {ref}')
        elif typ == 'order':
            st, n, beta = next(lines).split()
            if int(st) != 0 or int(n) != ref[0] or abs(float(beta)-ref[1]) > 1e-12:
                raise AssertionError(f'case {i} {command}: got {st} {n} {beta} versus {ref}')
            worst = max(worst, abs(float(beta)-ref[1]))
        else:
            st, n = map(int,next(lines).split())
            vals = np.fromstring(next(lines), sep=' ')
            if st != 0 or n != len(ref) or vals.size != n:
                raise AssertionError(f'case {i} bad FIR shape {st=} {n=} {vals.size=}')
            error = float(np.max(np.abs(vals-ref)))
            worst = max(worst,error)
            if not np.allclose(vals,ref, rtol=2e-12,atol=2e-12):
                raise AssertionError(f'case {i} {command} FIR error {error}')
    try:
        extra = next(lines)
    except StopIteration:
        pass
    else:
        raise AssertionError(f'unexpected output line: {extra}')
    print(f'PASS: {len(commands)}/{len(commands)} SciPy 1.17 comparisons, '
          f'maximum absolute difference {worst:.3g}')

if __name__ == '__main__':
    main()
