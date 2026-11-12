#!/usr/bin/env python3
"""Validate reusable banded LU factors in checked and optimized GNU Fortran builds."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
from compare_banded_lu_scipy import build_cases, run

ROOT=Path(__file__).resolve().parents[2]
SOURCES=['src/kinds.f90','src/linalg.f90','src/linalg_banded.f90']


def compile_program(work: Path, source: str, name: str, flags: list[str]) -> Path:
    exe=work/name
    subprocess.run(['gfortran','-std=f2018','-ffree-line-length-none',*flags,
                   *[str(ROOT/s) for s in SOURCES],str(ROOT/source),'-o',str(exe)],
                   cwd=work,check=True)
    return exe


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    parser.add_argument('--cases',type=int,default=500)
    parser.add_argument('--seed',type=int,default=171)
    parser.add_argument('--max-order',type=int,default=70)
    parser.add_argument('--stress',action='store_true')
    parser.add_argument('--tolerance',type=float,default=3.e-10)
    args=parser.parse_args()
    cases=build_cases(args.cases,args.seed,args.max_order,args.stress)
    modes=['checked','optimized'] if args.mode=='both' else [args.mode]
    with tempfile.TemporaryDirectory(prefix='scifort-banded-lu-') as tmp:
        work=Path(tmp)
        for mode in modes:
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            native=compile_program(work,'test/test_linalg_banded_lu.f90',f'native_{mode}',flags)
            completed=subprocess.run([str(native)],check=True,capture_output=True,text=True)
            print(mode+': '+completed.stdout.strip())
            driver=compile_program(work,'tools/linalg_reference/driver_banded_lu.f90',f'driver_{mode}',flags)
            run(driver,cases,args.tolerance)

if __name__=='__main__':
    main()
