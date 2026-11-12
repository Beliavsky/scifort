#!/usr/bin/env python3
"""Build and verify mixed 1D convolutions against SciPy 1.17 without FPM.

Builds GNU Fortran dependencies in clean temporary directories. The optimized
variant builds both edited modules and tests with -O3, leaving dependencies
at -O0. Native and independently computed SciPy reference tests are run.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = ('kinds', 'signal_windows_extra', 'signal_dpss', 'signal_spectral')
EDITED = ('signal_fftconvolve', 'signal_convolve')


def execute(args: list[str], cwd: Path, output=None) -> None:
    subprocess.run(args,cwd=cwd,stdout=output,check=True)


def build(compiler: str, destination: Path, optimized: bool) -> None:
    destination.mkdir(parents=True,exist_ok=True)
    flags=['-std=f2018','-ffree-line-length-none','-fimplicit-none',
           '-Werror=implicit-interface','-J.','-I.']
    checked=['-O0','-g','-fcheck=all','-fbacktrace']
    for module in DEPENDENCIES + EDITED:
        mode= ['-O3'] if optimized and module in EDITED else (
            ['-O0'] if optimized else checked)
        execute([compiler,*flags,*mode,'-c',str(ROOT/'src'/f'{module}.f90'),
                 '-o',f'{module}.o'],destination)
    objects=[f'{n}.o' for n in DEPENDENCIES + EDITED]
    mode=['-O3'] if optimized else checked
    for test in ('test_signal_mixed_convolution',):
        execute([compiler,*flags,*mode,str(ROOT/'test'/f'{test}.f90'),
                 *objects,'-o',test],destination)
        execute([str(destination/test)],destination)
    reference='signal_mixed_convolution_driver'
    execute([compiler,*flags,*mode,
             str(ROOT/'reference'/f'{reference}.f90'),*objects,'-o',reference],destination)
    data=destination/'mixed_reference.txt'
    with data.open('w',encoding='utf-8') as f:
        execute([str(destination/reference)],destination,output=f)
    execute([sys.executable,str(ROOT/'tools'/'compare_signal_mixed_convolution.py'),
             str(data)],destination)
    print(f"PASS: {'targeted -O3' if optimized else 'checked -fcheck=all'} mixed convolution suite",flush=True)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    flags=parser.add_mutually_exclusive_group()
    flags.add_argument('--checked-only',action='store_true')
    flags.add_argument('--optimized-only',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_mixed_convolution_') as work:
        folder=Path(work).resolve()
        if not args.optimized_only:
            build(args.compiler,folder/'checked',False)
        if not args.checked_only:
            build(args.compiler,folder/'optimized',True)

if __name__=='__main__':
    main()
