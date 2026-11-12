#!/usr/bin/env python3
"""Compile and validate SOS group delay against SciPy without requiring FPM."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MODULES = ('kinds', 'signal_windows_extra', 'signal_dpss', 'signal_spectral',
           'signal_analog_design', 'signal_design')


def run(command: list[str], work: Path) -> None:
    subprocess.run(command,cwd=work,check=True)


def check(compiler: str, work: Path, optimized: bool) -> None:
    work.mkdir()
    flags = ['-std=f2018','-fimplicit-none','-Werror=implicit-interface',
             '-I.','-J.']
    flags += ['-O3'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    objects = []
    for module in MODULES:
        obj = f'{module}.o'
        run([compiler,*flags,'-c',str(ROOT/'src'/f'{module}.f90'),'-o',obj],work)
        objects.append(obj)
    for test in ('test_signal_group_delay_sos','test_signal_group_delay'):
        run([compiler,*flags,str(ROOT/'test'/f'{test}.f90'),*objects,
             '-o',test],work)
        run([str(work/test)],work)
    driver = 'group_delay_sos_driver'
    run([compiler,*flags,str(ROOT/'tools'/'signal_reference'/f'{driver}.f90'),
         *objects,'-o',driver],work)
    run([sys.executable,str(ROOT/'tools'/'signal_reference'/
                            'compare_group_delay_sos_scipy.py'),str(work/driver)],work)
    print('PASS SOS group delay', 'optimized -O3' if optimized else 'checked -fcheck=all',
          flush=True)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--checked-only',action='store_true')
    mode.add_argument('--optimized-only',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_sos_group_delay_') as tmp:
        root=Path(tmp).resolve()
        if not args.optimized_only:
            check(args.compiler,root/'checked',False)
        if not args.checked_only:
            check(args.compiler,root/'optimized',True)


if __name__ == '__main__':
    main()
