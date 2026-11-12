#!/usr/bin/env python3
"""Build native Wiener regressions and compare complex/legacy filters with SciPy."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BASE = ('kinds', 'signal_denoise')


def command(args: list[str], folder: Path) -> None:
    subprocess.run(args, cwd=folder, check=True)


def run_variant(compiler: str, folder: Path, optimized: bool) -> None:
    """Build in an empty folder and run both new and existing tests."""
    folder.mkdir(parents=True, exist_ok=True)
    flags = ['-std=f2018', '-fimplicit-none', '-Werror=implicit-interface',
             '-ffree-line-length-none', '-J.', '-I.']
    flags.extend(['-O3'] if optimized else ['-O0', '-g', '-fcheck=all', '-fbacktrace'])
    for name in BASE:
        command([compiler, *flags, '-c', str(ROOT / 'src' / f'{name}.f90'),
                 '-o', f'{name}.o'], folder)
    objects = [f'{name}.o' for name in BASE]
    for test in ('test_signal_denoise', 'test_signal_wiener_complex'):
        command([compiler, *flags, str(ROOT / 'test' / f'{test}.f90'),
                 *objects, '-o', test], folder)
        command([str(folder / test)], folder)
    ref = 'signal_wiener_complex_reference_driver'
    command([compiler, *flags, str(ROOT / 'reference' / f'{ref}.f90'),
             *objects, '-o', ref], folder)
    with (folder / 'complex_results.txt').open('w') as out:
        subprocess.run([str(folder / ref)], cwd=folder, stdout=out, check=True)
    command([sys.executable, str(ROOT / 'tools' / 'compare_signal_wiener_complex.py'),
             str(folder / 'complex_results.txt')], folder)
    legacy = 'denoise_reference_driver'
    command([compiler, *flags, str(ROOT / 'tools' / 'signal_reference' / f'{legacy}.f90'),
             *objects, '-o', legacy], folder)
    command([sys.executable,
             str(ROOT / 'tools' / 'signal_reference' / 'compare_signal_denoise_scipy.py'),
             str(folder / legacy)], folder)
    print(f"PASS: {'optimized -O3' if optimized else 'checked -fcheck=all'} Wiener suites", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--checked-only', action='store_true')
    parser.add_argument('--optimized-only', action='store_true')
    opts = parser.parse_args()
    if opts.checked_only and opts.optimized_only:
        parser.error('Choose only one build mode')
    with tempfile.TemporaryDirectory(prefix='scifort_wiener_complex_') as scratch:
        folder = Path(scratch).resolve()
        if not opts.optimized_only:
            run_variant(opts.compiler, folder / 'checked', False)
        if not opts.checked_only:
            run_variant(opts.compiler, folder / 'optimized', True)


if __name__ == '__main__':
    main()
