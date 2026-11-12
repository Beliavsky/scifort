#!/usr/bin/env python3
"""Build and validate complex SOS decimation against native and SciPy reference tests.

The optimized mode recompiles the two modified production modules with -O3;
unchanged dependencies retain runtime-checking flags. Requires gfortran, SciPy,
and NumPy, with no SciPy dependency in the Fortran library itself.
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
MODIFIED = {'signal_resample', 'signal_multirate_nd_into'}
TESTS = (
    ('test/test_signal_decimate_complex_sos.f90', None),
    ('tools/reference_signal_decimate_complex_sos.f90',
     'tools/compare_signal_decimate_complex_sos.py'),
    ('test/test_signal_decimate_transfer_function.f90', None),
    ('tools/reference_signal_decimate_transfer_function.f90',
     'tools/compare_signal_decimate_transfer_function.py'),
    ('test/test_signal_decimate_custom_filters.f90', None),
    ('tools/reference_signal_decimate_custom_filters.f90',
     'tools/compare_signal_decimate_custom_filters.py'),
)


def validate(compiler: str, build: Path, optimized: bool,
             checked_build: Path | None = None) -> None:
    build.mkdir()
    common = ['-std=f2018', '-Werror=implicit-interface', '-I.', '-J.']
    checked = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    optimized_flags = ['-O3', '-fcheck=all', '-fbacktrace']
    objects = [name+'.o' for name in SOURCES]
    if checked_build is not None:
        # Retain unchanged checked dependencies without recompiling their
        # largest modules. Replace both modified objects and module files.
        for source in checked_build.iterdir():
            if source.suffix.lower() in ('.o', '.mod'):
                shutil.copy2(source, build/source.name)
    for name in SOURCES:
        if checked_build is not None and name not in MODIFIED:
            continue
        flags = optimized_flags if optimized and name in MODIFIED else checked
        subprocess.run([compiler, *common, *flags, '-c',
                        str(ROOT/'src'/f'{name}.f90'), '-o', name+'.o'],
                       cwd=build, check=True)
    test_flags = optimized_flags if optimized else checked
    for source, comparator in TESTS:
        executable = build/Path(source).stem
        subprocess.run([compiler, *common, *test_flags, str(ROOT/source),
                        *objects, '-o', str(executable)], cwd=build, check=True)
        if comparator:
            subprocess.run([sys.executable, str(ROOT/comparator), str(executable)],
                           cwd=build, check=True)
        else:
            subprocess.run([str(executable)], cwd=build, check=True)
    for driver, comparator in (
        ('multirate_nd_driver', 'check_multirate_nd'),
        ('multirate_complex_driver', 'check_multirate_complex'),
    ):
        executable = build/driver
        subprocess.run([compiler, *common, *test_flags,
                        str(ROOT/'tools'/'signal_reference'/f'{driver}.f90'),
                        *objects, '-o', str(executable)], cwd=build, check=True)
        subprocess.run([sys.executable,
                        str(ROOT/'tools'/'signal_reference'/f'{comparator}.py'),
                        str(executable)], cwd=build, check=True)
    print('PASS:', 'targeted optimized -O3' if optimized else 'checked',
          'complex-SOS decimation and prior multirate regressions',flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--checked-only', action='store_true')
    parser.add_argument('--optimized-only', action='store_true')
    args = parser.parse_args()
    if args.checked_only and args.optimized_only:
        parser.error('select at most one build mode')
    with tempfile.TemporaryDirectory(prefix='scifort_complex_sos_') as temp:
        checked = Path(temp)/'checked'
        if not args.optimized_only:
            validate(args.compiler, checked, False)
        if not args.checked_only:
            validate(args.compiler, Path(temp)/'optimized', True,
                     checked_build=checked if not args.optimized_only else None)


if __name__ == '__main__':
    main()
