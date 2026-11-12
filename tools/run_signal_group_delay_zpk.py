#!/usr/bin/env python3
"""Build checked/optimized GNU Fortran and validate ZPK group delay versus SciPy."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MODULES = ('kinds','signal_windows_extra','signal_dpss','signal_spectral',
           'signal_analog_design','signal_design')
TESTS = ('test_signal_group_delay_zpk','test_signal_group_delay_sos',
         'test_signal_group_delay')


def run(command: list[str], work: Path) -> None:
    subprocess.run(command,cwd=work,check=True)


def check(compiler: str, work: Path, optimized: bool) -> None:
    work.mkdir()
    flags = ['-std=f2018','-fimplicit-none','-Werror=implicit-interface','-I.','-J.']
    flags += ['-O3'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    objs = []
    for module in MODULES:
        obj = f'{module}.o'
        run([compiler,*flags,'-c',str(ROOT/'src'/f'{module}.f90'),'-o',obj],work)
        objs.append(obj)
    for test in TESTS:
        run([compiler,*flags,str(ROOT/'test'/f'{test}.f90'),*objs,'-o',test],work)
        run([str(work/test)],work)
    name = 'group_delay_zpk_driver'
    run([compiler,*flags,str(ROOT/'tools'/'signal_reference'/f'{name}.f90'),
         *objs,'-o',name],work)
    run([sys.executable,str(ROOT/'tools'/'signal_reference'/
                            'compare_group_delay_zpk_scipy.py'),str(work/name)],work)
    print('PASS ZPK group delay', 'optimized -O3' if optimized else 'checked -fcheck=all',
          flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--checked-only',action='store_true')
    mode.add_argument('--optimized-only',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_group_delay_zpk_') as tmp:
        work=Path(tmp).resolve()
        if not args.optimized_only:
            check(args.compiler,work/'checked',False)
        if not args.checked_only:
            check(args.compiler,work/'optimized',True)


if __name__ == '__main__':
    main()
