#!/usr/bin/env python3
"""Build/test shaped Fourier-resampling windows without FPM.

Requires GNU Fortran, NumPy, and SciPy. The optimized run recompiles the two
changed Fortran modules with -O3 and reuses checked unchanged dependencies.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'kinds', 'signal_windows_extra', 'signal_dpss', 'signal_spectral',
    'signal_analog_design', 'signal_design', 'signal_lfilter', 'signal_sos',
    'signal_filtfilt_gust', 'signal_filtfilt', 'signal_resample',
    'signal_fourier_resample', 'signal_multirate_nd_into',
)
MODIFIED = {'signal_fourier_resample', 'signal_multirate_nd_into'}


def build_and_check(fc: str, build: Path, optimized: bool,
                    checked: Path | None = None) -> None:
    build.mkdir()
    base_flags = ['-std=f2018', '-Werror=implicit-interface', '-I.', '-J.']
    checked_flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    optimized_flags = ['-O3', '-fcheck=all', '-fbacktrace']
    if checked is not None:
        for old in checked.iterdir():
            if old.suffix in ('.mod', '.o'):
                shutil.copy2(old, build/old.name)
    for name in SOURCES:
        if checked is not None and name not in MODIFIED:
            continue
        flags = optimized_flags if optimized and name in MODIFIED else checked_flags
        subprocess.run([fc, *base_flags, *flags, '-c',
                        str(ROOT/'src'/(name+'.f90')), '-o', name+'.o'],
                       cwd=build, check=True)
    objects = [name+'.o' for name in SOURCES]
    flags = optimized_flags if optimized else checked_flags
    for driver, comparator in (
        ('test/test_signal_resample_param_windows.f90', None),
        ('tools/reference_signal_resample_param_windows.f90',
         'tools/compare_signal_resample_param_windows.py'),
    ):
        exe = build / Path(driver).stem
        subprocess.run([fc, *base_flags, *flags, str(ROOT/driver),
                        *objects, '-o', str(exe)], cwd=build, check=True)
        if comparator:
            subprocess.run([sys.executable, str(ROOT/comparator), str(exe)],
                           cwd=build, check=True)
        else:
            subprocess.run([str(exe)], cwd=build, check=True)
    for driver, comparator in (
        ('multirate_nd_driver', 'check_multirate_nd'),
        ('multirate_complex_driver', 'check_multirate_complex'),
    ):
        exe=build/driver
        subprocess.run([fc, *base_flags, *flags,
                        str(ROOT/'tools'/'signal_reference'/(driver+'.f90')),
                        *objects, '-o', str(exe)], cwd=build, check=True)
        subprocess.run([sys.executable,
                        str(ROOT/'tools'/'signal_reference'/(comparator+'.py')),
                        str(exe)], cwd=build, check=True)
    print('PASS:', 'targeted -O3' if optimized else 'checked',
          'Fourier window and multirate regressions', flush=True)


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--checked-only',action='store_true')
    args=p.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_fourier_windows_') as temp:
        checked=Path(temp)/'checked'
        build_and_check(args.compiler,checked,optimized=False)
        if not args.checked_only:
            build_and_check(args.compiler,Path(temp)/'optimized',optimized=True,checked=checked)


if __name__=='__main__':
    main()
