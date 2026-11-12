#!/usr/bin/env python3
"""Check native SOS step/impulse responses in checked and optimized builds."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MODULES = ('kinds', 'signal_lfilter', 'signal_sos', 'signal_windows_extra',
           'signal_dpss', 'signal_spectral', 'signal_analog_design', 'signal_design',
           'signal_analog', 'signal_zpk_conversion', 'signal_state_space',
           'signal_discrete_simulation')
TESTS = ('test_signal_sos_responses', 'test_signal_group_delay_zpk',
         'test_signal_group_delay_sos', 'test_signal_group_delay')


def run(args: list[str], cwd: Path) -> None:
    subprocess.run(args, cwd=cwd, check=True)


def check(compiler: str, build: Path, optimized: bool) -> None:
    build.mkdir()
    flags = ['-std=f2018', '-fimplicit-none', '-Werror=implicit-interface', '-I.', '-J.']
    flags += ['-O3'] if optimized else ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    objects = []
    for source in MODULES:
        obj = f'{source}.o'
        run([compiler, *flags, '-c', str(ROOT/'src'/f'{source}.f90'), '-o', obj], build)
        objects.append(obj)
    for program in TESTS:
        run([compiler, *flags, str(ROOT/'test'/f'{program}.f90'),
             *objects, '-o', program], build)
        run([str(build/program)], build)
    reference = 'sos_response_driver'
    run([compiler, *flags, str(ROOT/'tools'/'signal_reference'/f'{reference}.f90'),
         *objects, '-o', reference], build)
    run([sys.executable, str(ROOT/'tools'/'signal_reference'/
                             'compare_sos_responses_scipy.py'),str(build/reference)], build)
    print('PASS SOS responses', '-O3' if optimized else '-fcheck=all',flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--checked-only', action='store_true')
    mode.add_argument('--optimized-only', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_sos_responses_') as tmp:
        work=Path(tmp).resolve()
        if not args.optimized_only:
            check(args.compiler,work/'checked',False)
        if not args.checked_only:
            check(args.compiler,work/'optimized',True)


if __name__ == '__main__':
    main()
