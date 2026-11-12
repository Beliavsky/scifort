#!/usr/bin/env python3
"""Compile and check the standalone windowed SciFort spectral slice.

Usage:
    python tools/run_signal_windowed_estimators.py
    python tools/run_signal_windowed_estimators.py --compiler gfortran --checked-only

Needs a Fortran 2018 compiler and (for the numerical reference) SciPy 1.17+
and NumPy. Does not require FPM or download any Fortran dependencies.
Builds in a temporary directory and never modifies the source tree.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'kinds',
    'signal_windows_extra',
    'signal_dpss',
    'signal_spectral',
    'signal_spectral_complex',
    'signal_spectral_nd',
    'signal_timefreq_nd',
)


def run(command: list[str], *, cwd: Path) -> None:
    print(' '.join(map(str, command)), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def compile_and_test(compiler: str, tmp: Path, *, optimized: bool) -> None:
    profile = 'optimized' if optimized else 'checked'
    build = tmp / profile
    build.mkdir()
    flags = ['-std=f2018', '-ffree-line-length-none', '-fimplicit-none',
             '-Werror=implicit-interface', '-J.', '-I.']
    if optimized:
        flags += ['-O3']
    else:
        flags += ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    objects = []
    for name in SOURCES:
        obj = name + '.o'
        run([compiler, *flags, '-c', str(ROOT / 'src' / (name + '.f90')),
             '-o', obj], cwd=build)
        objects.append(obj)
    extension = '.exe' if os.name == 'nt' else ''
    for name in ('test_signal_windowed_estimators', 'test_signal_stft_boundaries'):
        executable = build / (name + extension)
        run([compiler, *flags, str(ROOT / 'test' / (name + '.f90')),
             *objects, '-o', str(executable)], cwd=build)
        run([str(executable)], cwd=build)
    driver = build / ('reference_signal_windowed_estimators' + extension)
    run([compiler, *flags, str(ROOT / 'tools' / 'reference_signal_windowed_estimators.f90'),
         *objects, '-o', str(driver)], cwd=build)
    run([sys.executable, str(ROOT / 'tools' / 'compare_signal_windowed_estimators.py'),
         str(driver)], cwd=build)
    print(f'PASS: {profile} windowed spectral estimators', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--checked-only', action='store_true',
                        help='skip optimized build (still run checked SciPy comparison)')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_windows_') as td:
        compile_and_test(args.compiler, Path(td), optimized=False)
        if not args.checked_only:
            compile_and_test(args.compiler, Path(td), optimized=True)


if __name__ == '__main__':
    main()
