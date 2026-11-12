#!/usr/bin/env python3
"""Checked and targeted optimized GNU Fortran / SciPy transfer-function regression."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('kinds', 'signal_windows_extra', 'signal_dpss', 'signal_spectral',
           'signal_analog_design', 'signal_design', 'signal_lfilter',
           'signal_sos', 'signal_filtfilt_gust', 'signal_filtfilt',
           'signal_resample', 'signal_fourier_resample', 'signal_multirate_nd_into')


def run(compiler: str, build: Path, optimized: bool) -> None:
    build.mkdir()
    shared = ['-std=f2018', '-Werror=implicit-interface', '-I.', '-J.']
    checked = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    # The unchanged dependencies are checked; both changed modules use -O3.
    for source in SOURCES:
        flags = ['-O3', '-fcheck=all'] if optimized and source in {
            'signal_resample', 'signal_multirate_nd_into'} else checked
        subprocess.run([compiler, *shared, *flags, '-c',
                        str(ROOT/'src'/f'{source}.f90'), '-o', f'{source}.o'],
                       cwd=build, check=True)
    objects = [f'{name}.o' for name in SOURCES]
    flags = ['-O3', '-fcheck=all'] if optimized else checked
    for source, compare in (
        ('test/test_signal_decimate_transfer_function.f90', None),
        ('tools/reference_signal_decimate_transfer_function.f90',
         'tools/compare_signal_decimate_transfer_function.py'),
        ('test/test_signal_decimate_custom_filters.f90', None),
        ('tools/reference_signal_decimate_custom_filters.f90',
         'tools/compare_signal_decimate_custom_filters.py'),
    ):
        executable = build/Path(source).stem
        subprocess.run([compiler, *shared, *flags, str(ROOT/source), *objects,
                        '-o', str(executable)], cwd=build, check=True)
        if compare is None:
            subprocess.run([str(executable)], cwd=build, check=True)
        else:
            subprocess.run([sys.executable, str(ROOT/compare), str(executable)],
                           cwd=build, check=True)
    for driver, checker in (
        ('multirate_nd_driver', 'check_multirate_nd'),
        ('multirate_complex_driver', 'check_multirate_complex'),
    ):
        exe = build/driver
        subprocess.run([compiler, *shared, *flags,
                        str(ROOT/'tools'/'signal_reference'/f'{driver}.f90'),
                        *objects, '-o', str(exe)], cwd=build, check=True)
        subprocess.run([sys.executable,
                        str(ROOT/'tools'/'signal_reference'/f'{checker}.py'),
                        str(exe)], cwd=build, check=True)
    print('PASS:', 'targeted -O3' if optimized else 'checked',
          'transfer-function decimation and previous multirate regressions', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--checked-only', action='store_true')
    parser.add_argument('--optimized-only', action='store_true')
    args = parser.parse_args()
    if args.checked_only and args.optimized_only:
        parser.error('select at most one build mode')
    with tempfile.TemporaryDirectory(prefix='scifort_decimate_transfer_') as temp:
        if not args.optimized_only:
            run(args.compiler, Path(temp)/'checked', optimized=False)
        if not args.checked_only:
            run(args.compiler, Path(temp)/'optimized', optimized=True)


if __name__ == '__main__':
    main()
