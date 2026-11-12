#!/usr/bin/env python3
"""Reproducible checked/O3 two-sided signal builds and SciPy 1.17 comparisons.

python tools/run_signal_twosided.py
python tools/run_signal_twosided.py --checked-only

This standalone slice needs GNU Fortran, Python, NumPy, and SciPy. No FPM,
network, or third-party Fortran libraries are required.
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
    'kinds', 'signal_windows_extra', 'signal_dpss',
    'signal_spectral', 'signal_spectral_complex', 'signal_spectral_nd',
    'signal_timefreq_nd', 'signal_istft', 'signal_istft_nd',
)
TESTS = (
    'test_signal_twosided', 'test_signal_windowed_estimators',
    'test_signal_stft_boundaries', 'test_signal_istft_parameter_windows',
    'test_signal_istft',
)


def run(argv: list[str], cwd: Path) -> None:
    subprocess.run(argv, cwd=cwd, check=True)


def verify(compiler: str, build_dir: Path, optimized: bool) -> None:
    build_dir.mkdir(parents=True)
    options = [
        '-std=f2018', '-ffree-line-length-none', '-fimplicit-none',
        '-Werror=implicit-interface', '-J.', '-I.',
    ]
    options += ['-O3'] if optimized else ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    objects = []
    for source in SOURCES:
        target = source + '.o'
        run([compiler, *options, '-c', str(ROOT / 'src' / (source + '.f90')),
             '-o', target], build_dir)
        objects.append(target)
    exe_suffix = '.exe' if os.name == 'nt' else ''
    for name in TESTS:
        executable = build_dir / (name + exe_suffix)
        run([compiler, *options, str(ROOT / 'test' / (name + '.f90')),
             *objects, '-o', str(executable)], build_dir)
        run([str(executable)], build_dir)
    executable = build_dir / ('reference_signal_twosided' + exe_suffix)
    run([compiler, *options, str(ROOT / 'tools' / 'reference_signal_twosided.f90'),
         *objects, '-o', str(executable)], build_dir)
    run([sys.executable, str(ROOT / 'tools' / 'compare_signal_twosided.py'),
         str(executable)], build_dir)
    print('PASS:', 'optimized -O3' if optimized else 'checked -fcheck=all',
          'two-sided spectral suite', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--checked-only', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_two_sided_') as temp:
        verify(args.compiler, Path(temp)/'checked', optimized=False)
        if not args.checked_only:
            verify(args.compiler, Path(temp)/'optimized', optimized=True)


if __name__ == '__main__':
    main()
