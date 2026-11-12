#!/usr/bin/env python3
"""Compile/test spectral detrend callbacks and compare their results with SciPy."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'kinds', 'signal_windows_extra', 'signal_dpss', 'signal_spectral',
    'signal_spectral_complex', 'signal_spectral_nd', 'signal_timefreq_nd',
)
OBJECTS = tuple(f'{name}.o' for name in SOURCES)


def run_variant(compiler: str, folder: Path, optimized: bool) -> None:
    """Compile dependencies, run native tests, and check SciPy output."""
    folder.mkdir(parents=True, exist_ok=True)
    flags = [
        '-std=f2018', '-fimplicit-none', '-Werror=implicit-interface',
        '-ffree-line-length-none', '-I.', '-J.',
    ]
    if optimized:
        flags += ['-O3']
    else:
        flags += ['-O0', '-g', '-fcheck=all', '-fbacktrace']

    for name, output in zip(SOURCES, OBJECTS):
        subprocess.run(
            [compiler, *flags, '-c', str(ROOT / 'src' / f'{name}.f90'),
             '-o', output],
            cwd=folder, check=True,
        )

    native = 'test_signal_spectral_callback_nd'
    subprocess.run(
        [compiler, *flags, str(ROOT / 'test' / f'{native}.f90'),
         *OBJECTS, '-o', native],
        cwd=folder, check=True,
    )
    subprocess.run([str(folder / native)], cwd=folder, check=True)

    reference = 'signal_detrend_callback_reference_driver'
    subprocess.run(
        [compiler, *flags, str(ROOT / 'reference' / f'{reference}.f90'),
         *OBJECTS, '-o', reference],
        cwd=folder, check=True,
    )
    output = folder / 'reference_results.txt'
    with output.open('w') as stream:
        subprocess.run([str(folder / reference)], cwd=folder, stdout=stream, check=True)
    subprocess.run(
        [sys.executable, str(ROOT / 'tools' / 'compare_signal_detrend_callbacks.py'),
         str(output)], cwd=folder, check=True,
    )
    print(f"PASS: {'optimized -O3' if optimized else 'checked -fcheck=all'} callback tests", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--checked-only', action='store_true')
    parser.add_argument('--optimized-only', action='store_true')
    args = parser.parse_args()
    if args.checked_only and args.optimized_only:
        parser.error('Select at most one of --checked-only and --optimized-only')
    with tempfile.TemporaryDirectory(prefix='scifort_spectral_callbacks_') as tmp:
        tmp = Path(tmp).resolve()
        if not args.optimized_only:
            run_variant(args.compiler, tmp / 'checked', False)
        if not args.checked_only:
            run_variant(args.compiler, tmp / 'optimized', True)


if __name__ == '__main__':
    main()
