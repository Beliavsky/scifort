#!/usr/bin/env python3
"""Validate multichannel SOS simulation against SciPy in checked/-O3 modes."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MODULES = (
    'kinds', 'signal_lfilter', 'signal_sos', 'signal_windows_extra',
    'signal_dpss', 'signal_spectral', 'signal_analog_design', 'signal_design',
    'signal_analog', 'signal_zpk_conversion', 'signal_state_space',
    'signal_discrete_simulation',
)
TESTS = ('test_signal_dlsim_sos_bank', 'test_signal_dlsim_sos',
         'test_signal_sos_responses')
COMPARISONS = (
    ('sos_bank_driver', 'compare_sos_bank_scipy.py'),
    ('sos_simulation_driver', 'compare_sos_simulation_scipy.py'),
    ('sos_response_driver', 'compare_sos_responses_scipy.py'),
)


def run(command: list[str], cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def build_and_test(compiler: str, path: Path, optimized: bool) -> None:
    path.mkdir()
    flags = ['-std=f2018', '-fimplicit-none', '-Werror=implicit-interface',
             '-I.', '-J.']
    base_flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    changed_flags = ['-O3'] if optimized else base_flags
    objects = []
    for source in MODULES:
        obj = source + '.o'
        module_flags = changed_flags if source == 'signal_discrete_simulation' else base_flags
        run([compiler, *flags, *module_flags, '-c',
             str(ROOT / 'src' / f'{source}.f90'), '-o', obj], path)
        objects.append(obj)
    for test in TESTS:
        run([compiler, *flags, *changed_flags,
             str(ROOT / 'test' / f'{test}.f90'), *objects, '-o', test], path)
        run([str(path / test)], path)
    for driver, comparer in COMPARISONS:
        run([compiler, *flags, *changed_flags,
             str(ROOT / 'tools' / 'signal_reference' / f'{driver}.f90'),
             *objects, '-o', driver], path)
        run([sys.executable, str(ROOT / 'tools' / 'signal_reference' / comparer),
             str(path / driver)], path)
    print('PASS multichannel SOS simulation:',
          'targeted -O3' if optimized else 'checked', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--checked-only', action='store_true')
    group.add_argument('--optimized-only', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_dlsim_bank_') as tmp:
        root = Path(tmp).resolve()
        if not args.optimized_only:
            build_and_test(args.compiler, root / 'checked', False)
        if not args.checked_only:
            build_and_test(args.compiler, root / 'optimized', True)


if __name__ == '__main__':
    main()
