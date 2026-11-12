#!/usr/bin/env python3
"""Check SciFort scalar/logical/multiband firwin against SciPy.

Usage: python tools/run_signal_firwin_multiband.py [--checked-only]

Builds the independent SciFort FIR kernel, native assertions, and a seeded
reference driver in both checked and optimized GNU Fortran configurations.
No FPM or external Fortran library is required.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('kinds','signal_windows_extra','signal_dpss','signal_spectral',
           'signal_analog_design','signal_design')


def check(compiler: str, build: Path, optimized: bool) -> None:
    build.mkdir()
    flags = ['-std=f2018','-fimplicit-none','-Werror=implicit-interface',
             '-ffree-line-length-none','-I.','-J.']
    flags += ['-O3'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    objects = []
    for source in SOURCES:
        name = source+'.o'
        subprocess.run([compiler,*flags,'-c',str(ROOT/'src'/f'{source}.f90'),
                        '-o',name],cwd=build,check=True)
        objects.append(name)
    for name,srcdir in [('test_signal_firwin_multiband','test'),
                        ('reference_signal_firwin_multiband','tools')]:
        executable = build/name
        subprocess.run([compiler,*flags,str(ROOT/srcdir/f'{name}.f90'),
                        *objects,'-o',str(executable)],cwd=build,check=True)
        if srcdir=='test':
            subprocess.run([str(executable)],cwd=build,check=True)
        else:
            subprocess.run([sys.executable,str(ROOT/'tools'/'compare_signal_firwin_multiband.py'),
                            str(executable)],cwd=build,check=True)
    print('PASS: FIR multiband', 'optimized -O3' if optimized else 'checked -fcheck=all',flush=True)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--checked-only',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_firwin_multi_') as d:
        check(args.compiler,Path(d)/'checked',False)
        if not args.checked_only:
            check(args.compiler,Path(d)/'optimized',True)


if __name__=='__main__':
    main()
