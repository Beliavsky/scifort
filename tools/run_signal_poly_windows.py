#!/usr/bin/env python3
"""Build and validate SciFort resample_poly named-window parameters (without FPM)."""
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
    flags = ['-std=f2018', '-Werror=implicit-interface', '-I.', '-J.']
    debug = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    # Compiling the very large, unchanged FFT/spectral dependencies at -O3
    # takes substantially longer than these tests. Optimize the two changed
    # modules, retaining checked compilation of their unmodified dependencies.
    changed = {'signal_resample', 'signal_multirate_nd_into'}
    objects = []
    for name in SOURCES:
        obj = f'{name}.o'
        source_flags = ['-O3', '-fcheck=all'] if optimize and name in changed else debug
        subprocess.run([compiler, *flags, *source_flags, '-c', str(ROOT/'src'/f'{name}.f90'),
                        '-o', obj], cwd=build, check=True)
        objects.append(obj)
    if optimize:
        flags += ['-O3', '-fcheck=all']
    else:
        flags += debug
    for name, directory in [('test_signal_resample_poly_param_windows', 'test'),
                            ('reference_signal_poly_windows', 'tools')]:
        exe = build/name
        subprocess.run([compiler, *flags, str(ROOT/directory/f'{name}.f90'),
                        *objects, '-o', str(exe)], cwd=build, check=True)
        if directory == 'test':
            subprocess.run([str(exe)], cwd=build, check=True)
        else:
            subprocess.run([sys.executable, str(ROOT/'tools'/'compare_signal_poly_windows.py'),
                            str(exe)], cwd=build, check=True)
    print('PASS:', 'targeted -O3' if optimize else 'checked', 'polyphase window tests', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--checked-only', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_poly_windows_') as temp:
        check(args.compiler, Path(temp)/'checked', False)
        if not args.checked_only:
            check(args.compiler, Path(temp)/'optimized', True)


if __name__ == '__main__':
    main()
