# SPDX-License-Identifier: MIT
"""Check modified Fresnel phases against independent high-precision erfc values.

Requires mpmath and an existing FPM build. Pass that build's module directory
and libscifort.a using --modules and --library. Only the changed Faddeeva module,
a small input driver, and four related native tests are compiled here.
"""

import argparse
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
DRIVER = """program phase_driver
    use scifort_kinds, only : dp
    use scifort_faddeeva, only : modfresnelp, modfresnelm
    implicit none
    real(dp) :: x
    complex(dp) :: fp, kp, fm, km
    integer :: ios
    do
        read(*, *, iostat=ios) x
        if (ios < 0) exit
        if (ios /= 0) error stop 1
        call modfresnelp(x, fp, kp)
        call modfresnelm(x, fm, km)
        write(*, '(8es26.17e3)') real(fp, dp), aimag(fp), real(kp, dp), aimag(kp), &
            real(fm, dp), aimag(fm), real(km, dp), aimag(km)
    end do
end program phase_driver
"""


def reference(x):
    """Rotate the defining Fresnel integral to erfc; no SciFort kernels used."""
    x = mp.mpf(x)  # exact binary64 input, not its short decimal representation
    phase = mp.exp(mp.j * mp.pi / 4)
    fp = mp.sqrt(mp.pi) / 2 * phase * mp.erfc(x / phase)
    kp = mp.exp(-mp.j * x*x) / phase * fp / mp.sqrt(mp.pi)
    return fp, kp, mp.conj(fp), mp.conj(kp)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modules', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    args = parser.parse_args()
    modules = args.modules.resolve()
    library = args.library.resolve()
    subprocess.run(['gfortran', '--version'], check=True)
    mp.mp.dps = 400  # resolves exact squares even near binary64 overflow
    rng = random.Random(20261119)
    positive = [0.0, 0.125, 1.0, 3.0, 6.0, 10.0, 67108864.125,
                100000000.25, 1000000000000.125, 1e100, 1e150,
                math.sqrt(float.fromhex('0x1.fffffffffffffp1023')) / 2]
    positive += [math.nextafter(x, math.inf) for x in positive[6:10]]
    positive += [10**rng.uniform(-2, 153) for _ in range(64)]
    xs = positive + [-x for x in positive if x != 0]
    refs = [reference(x) for x in xs]
    payload = ''.join(f'{x:.17e}\n' for x in xs)
    native_tests = ['modified_fresnel', 'faddeeva', 'complex_error', 'complex_fresnel']
    for profile, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                      '-ffpe-trap=invalid,zero,overflow', '-fbacktrace']),
                           ('release', ['-O3', '-funroll-loops'])]:
        directory = ROOT / 'build' / 'modified-fresnel-phase' / profile
        directory.mkdir(parents=True, exist_ok=True)
        obj = directory / 'faddeeva.o'
        common = ['gfortran', *flags, '-I', str(directory), '-I', str(modules)]
        subprocess.run([*common, '-J', str(directory), '-c',
                        str(ROOT / 'src/special/faddeeva.f90'), '-o', str(obj)], check=True)
        for test in native_tests:
            exe = directory / (test + '.exe')
            subprocess.run([*common, str(ROOT / f'test/test_special_{test}.f90'),
                            str(obj), str(library), '-o', str(exe)], check=True)
            subprocess.run([str(exe)], check=True)
        source = directory / 'driver.f90'
        source.write_text(DRIVER, encoding='ascii')
        exe = directory / 'driver.exe'
        subprocess.run([*common, str(source), str(obj), str(library), '-o', str(exe)],
                       check=True)
        result = subprocess.run([str(exe)], input=payload, text=True,
                                capture_output=True, check=True)
        lines = result.stdout.splitlines()
        assert len(lines) == len(xs), (len(lines), len(xs))
        worst = mp.mpf(0)
        for x, expected, line in zip(xs, refs, lines):
            parts = [mp.mpf(s) for s in line.split()]
            assert len(parts) == 8, line
            for j, ref in enumerate(expected):
                actual = mp.mpc(parts[2*j], parts[2*j+1])
                error = abs(actual - ref) / abs(ref)
                assert error < mp.mpf('5e-14'), (profile, x, j, error)
                worst = max(worst, error)
        print(f'{profile}: {len(xs)} inputs, {4*len(xs)} complex comparisons passed; '
              f'maximum relative error {float(worst):.6e}', flush=True)


if __name__ == '__main__':
    main()
