#!/usr/bin/env python3
"""Build and compare dlsim_sos against SciPy (checked and optimized modes)."""
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
TESTS = ('test_signal_dlsim_sos', 'test_signal_sos_responses')


def run(argv: list[str], cwd: Path) -> None:
    subprocess.run(argv, cwd=cwd, check=True)


def build_and_test(compiler: str, build: Path, optimized: bool) -> None:
    build.mkdir()
    flags = ['-std=f2018', '-fimplicit-none', '-Werror=implicit-interface', '-I.', '-J.']
    checked_flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    targeted_flags = ['-O3'] if optimized else checked_flags
    objs = []
    for name in MODULES:
        obj = name + '.o'
        # Keep the large unchanged dependencies checked in optimized mode.
        mod_flags = targeted_flags if name == 'signal_discrete_simulation' else checked_flags
        run([compiler, *flags, *mod_flags, '-c',
             str(ROOT/'src'/f'{name}.f90'), '-o', obj], build)
        objs.append(obj)
    for name in TESTS:
        run([compiler, *flags, *targeted_flags,
             str(ROOT/'test'/f'{name}.f90'), *objs, '-o', name], build)
        run([str(build/name)], build)
    ref_dir = ROOT/'tools'/'signal_reference'
    for driver, compare in (
            ('sos_simulation_driver', 'compare_sos_simulation_scipy.py'),
            ('sos_response_driver', 'compare_sos_responses_scipy.py')):
        run([compiler, *flags, *targeted_flags,
             str(ref_dir/f'{driver}.f90'), *objs, '-o', driver], build)
        run([sys.executable, str(ref_dir/compare), str(build/driver)], build)
    print('PASS SOS arbitrary-input simulation',
          'targeted -O3' if optimized else '-fcheck=all', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--checked-only', action='store_true')
    mode.add_argument('--optimized-only', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_dlsim_sos_') as tmp:
        root = Path(tmp).resolve()
        if not args.optimized_only:
            build_and_test(args.compiler, root/'checked', False)
        if not args.checked_only:
            build_and_test(args.compiler, root/'optimized', True)


if __name__ == '__main__':
    main()
