#!/usr/bin/env python3
"""Build and compare SciFort decimation with caller-supplied FIR/SOS filters."""
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


def check(compiler: str, build: Path, optimize: bool) -> None:
    build.mkdir()
    shared = ['-std=f2018', '-Werror=implicit-interface', '-I.', '-J.']
    debug = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    changed = {'signal_resample', 'signal_multirate_nd_into'}
    objs = []
    for name in SOURCES:
        flags = ['-O3', '-fcheck=all'] if optimize and name in changed else debug
        subprocess.run([compiler, *shared, *flags, '-c', str(ROOT/'src'/f'{name}.f90'),
                        '-o', name+'.o'], cwd=build, check=True)
        objs.append(name+'.o')
    flags = ['-O3', '-fcheck=all'] if optimize else debug
    for source in ('test/test_signal_decimate_custom_filters.f90',
                   'tools/reference_signal_decimate_custom_filters.f90'):
        exe = build/Path(source).stem
        subprocess.run([compiler, *shared, *flags, str(ROOT/source),
                        *objs, '-o', str(exe)], cwd=build, check=True)
        if source.startswith('test/'):
            subprocess.run([str(exe)], cwd=build, check=True)
        else:
            subprocess.run([sys.executable, str(ROOT/'tools'/'compare_signal_decimate_custom_filters.py'),
                            str(exe)], cwd=build, check=True)
    for driver, checker in (
        ('multirate_nd_driver', 'check_multirate_nd'),
        ('multirate_complex_driver', 'check_multirate_complex'),
    ):
        exe = build/driver
        subprocess.run([compiler, *shared, *flags, str(ROOT/'tools'/'signal_reference'/f'{driver}.f90'),
                        *objs, '-o', str(exe)], cwd=build, check=True)
        subprocess.run([sys.executable, str(ROOT/'tools'/'signal_reference'/f'{checker}.py'),
                        str(exe)], cwd=build, check=True)
    print('PASS:', 'targeted -O3' if optimize else 'checked',
          'custom FIR/SOS decimation', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--checked-only', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_custom_decimate_') as temp:
        check(args.compiler, Path(temp)/'checked', False)
        if not args.checked_only:
            check(args.compiler, Path(temp)/'optimized', True)


if __name__ == '__main__':
    main()
