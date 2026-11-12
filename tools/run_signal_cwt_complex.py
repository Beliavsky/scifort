#!/usr/bin/env python3
"""Validate complex Morlet and callback CWT with SciPy, no FPM needed.

Builds the two local Fortran modules and the native and reference drivers
from clean temporary directories with both runtime checks and optimization.
Use the general signal standalone runner for the large public module tests.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def execute(command: list[str], folder: Path, output=None) -> None:
    subprocess.run(command, cwd=folder, stdout=output, check=True)


def check(compiler: str, folder: Path, optimized: bool) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    base = ['-std=f2018', '-fimplicit-none', '-Werror=implicit-interface', '-J.', '-I.']
    flags = ['-O3'] if optimized else ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    objects = ['kinds.o', 'signal_cwt_peaks.o']
    for name in ('kinds', 'signal_cwt_peaks'):
        execute([compiler, *base, *flags, '-c', str(ROOT / 'src' / f'{name}.f90'),
                 '-o', f'{name}.o'], folder)
    native = 'test_signal_cwt_complex'
    execute([compiler, *base, *flags, str(ROOT / 'test' / f'{native}.f90'),
             *objects, '-o', native], folder)
    execute([str(folder / native)], folder)
    driver = 'signal_cwt_complex_driver'
    execute([compiler, *base, *flags, str(ROOT / 'reference' / f'{driver}.f90'),
             *objects, '-o', driver], folder)
    path = folder / 'cwt_complex_reference.txt'
    with path.open('w', encoding='utf-8') as file:
        execute([str(folder / driver)], folder, output=file)
    execute([sys.executable, str(ROOT / 'tools' / 'compare_signal_cwt_complex.py'),
             str(path)], folder)
    print('PASS complex CWT ' + ('optimized -O3' if optimized else 'checked -fcheck=all'), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--checked-only', action='store_true')
    modes.add_argument('--optimized-only', action='store_true')
    opts = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_cwt_complex_') as folder:
        root = Path(folder).resolve()
        if not opts.optimized_only:
            check(opts.compiler, root / 'checked', False)
        if not opts.checked_only:
            check(opts.compiler, root / 'optimized', True)


if __name__ == '__main__':
    main()
