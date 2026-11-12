#!/usr/bin/env python3
"""Build and compare ShortTimeFFT callable detrending with SciPy.

Requires GNU Fortran, Python, NumPy and SciPy. No FPM is necessary.
Run checked and targeted optimized builds in separate clean directories.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'kinds', 'signal_windows_extra', 'signal_dpss',
    'signal_spectral', 'signal_spectral_complex', 'signal_short_time_fft',
)
NATIVE = (
    'test_signal_stft_callable_detrend',
    'test_signal_short_time_fft',
    'test_signal_short_time_fft_helpers',
    'test_signal_short_time_fft_nd',
    'test_signal_short_time_fft_rank45',
    'test_signal_short_time_fft_high_rank',
    'test_signal_mixed_cross',
)
REFERENCE = 'signal_stft_callable_detrend_driver'


def run_variant(compiler: str, folder: Path, optimized: bool) -> None:
    """Compile an independent test suite and verify SciPy reference values."""
    folder.mkdir(parents=True, exist_ok=True)
    flags = ['-std=f2018', '-ffree-line-length-none', '-fimplicit-none',
             '-Werror=implicit-interface', '-J.', '-I.']
    checked = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    objects: list[str] = []
    for source in SOURCES:
        # Optimize the changed ShortTimeFFT implementation, but leave
        # unchanged large dependency modules at -O0 in targeted -O3 builds.
        mode = (['-O3'] if optimized and source == 'signal_short_time_fft' else
                ['-O0'] if optimized else checked)
        obj = f'{source}.o'
        subprocess.run([compiler, *flags, *mode, '-c',
                        str(ROOT / 'src' / f'{source}.f90'), '-o', obj],
                       cwd=folder, check=True)
        objects.append(obj)
    test_flags = ['-O3'] if optimized else checked
    for test in NATIVE:
        subprocess.run([compiler, *flags, *test_flags,
                        str(ROOT / 'test' / f'{test}.f90'),
                        *objects, '-o', test], cwd=folder, check=True)
        subprocess.run([str(folder / test)], cwd=folder, check=True)
    subprocess.run([compiler, *flags, *test_flags,
                    str(ROOT / 'reference' / f'{REFERENCE}.f90'),
                    *objects, '-o', REFERENCE], cwd=folder, check=True)
    data = folder / 'stft_callable_results.txt'
    with data.open('w', encoding='utf-8') as out:
        subprocess.run([str(folder / REFERENCE)], cwd=folder,
                       stdout=out, check=True)
    subprocess.run([sys.executable,
                    str(ROOT / 'tools' / 'compare_signal_stft_callable_detrend.py'),
                    str(data)], cwd=folder, check=True)
    print(f"PASS: {'targeted -O3' if optimized else 'checked -fcheck=all'} "
          'ShortTimeFFT callable tests and prior regressions', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--checked-only', action='store_true')
    modes.add_argument('--optimized-only', action='store_true')
    options = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_stft_callable_') as scratch:
        build = Path(scratch).resolve()
        if not options.optimized_only:
            run_variant(options.compiler, build / 'checked', False)
        if not options.checked_only:
            run_variant(options.compiler, build / 'optimized', True)


if __name__ == '__main__':
    main()
