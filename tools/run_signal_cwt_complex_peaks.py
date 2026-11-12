#!/usr/bin/env python3
"""Build and validate complex magnitude ridge peaks and legacy wavelet behavior.

Requires gfortran, NumPy, and SciPy (tested with SciPy 1.17); no FPM needed.
Includes both a checked build and a fully optimized two-module -O3 build.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / 'tools' / 'signal_reference'


def execute(command: list[str], work: Path, output=None) -> None:
    subprocess.run(command, cwd=work, stdout=output, check=True)


def check(compiler: str, work: Path, optimized: bool) -> None:
    work.mkdir(parents=True, exist_ok=True)
    base = ['-std=f2018', '-fimplicit-none', '-Werror=implicit-interface', '-J.', '-I.']
    flags = ['-O3'] if optimized else ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    objects = ['kinds.o', 'signal_cwt_peaks.o']
    for name in ('kinds', 'signal_cwt_peaks'):
        execute([compiler, *base, *flags, '-c', str(ROOT/'src'/f'{name}.f90'),
                 '-o', f'{name}.o'], work)

    # Native tests cover malformed callbacks and status behavior, while the
    # independent numerical drivers cover CWT values and ridge locations.
    for name in ('test_signal_cwt_complex_peaks', 'test_signal_cwt_complex'):
        execute([compiler, *base, *flags, str(ROOT/'test'/f'{name}.f90'),
                 *objects, '-o', name], work)
        execute([str(work/name)], work)

    driver = 'cwt_complex_peaks_driver'
    execute([compiler, *base, *flags, str(REFERENCE/f'{driver}.f90'),
             *objects, '-o', driver], work)
    execute([sys.executable, str(REFERENCE/'compare_cwt_complex_peaks_scipy.py'),
             str(work/driver)], work)

    # The original real Ricker ridge detector must retain exact SciPy
    # behavior across all 300 existing seeded cases after ridge-code reuse.
    driver = 'cwt_peaks_driver'
    execute([compiler, *base, *flags, str(REFERENCE/f'{driver}.f90'),
             *objects, '-o', driver], work)
    sys.path.insert(0, str(REFERENCE))
    from verify_cwt_place_poles_scipy import check_wavelets
    check_wavelets(work/driver)

    # Also recheck the previous 240 complex CWT convolution comparisons.
    driver = 'signal_cwt_complex_driver'
    execute([compiler, *base, *flags, str(ROOT/'reference'/f'{driver}.f90'),
             *objects, '-o', driver], work)
    reference_file = work / 'complex_cwt_reference.txt'
    with reference_file.open('w', encoding='utf-8') as handle:
        execute([str(work/driver)], work, output=handle)
    execute([sys.executable, str(ROOT/'tools'/'compare_signal_cwt_complex.py'),
             str(reference_file)], work)
    print('PASS complex CWT ridge peaks: ' +
          ('fully optimized -O3' if optimized else 'checked -fcheck=all'), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--checked-only', action='store_true')
    modes.add_argument('--optimized-only', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_complex_ridges_') as tmp:
        root = Path(tmp).resolve()
        if not args.optimized_only:
            check(args.compiler, root/'checked', False)
        if not args.checked_only:
            check(args.compiler, root/'optimized', True)


if __name__ == '__main__':
    main()
