#!/usr/bin/env python3
"""Build and test real/complex packed-band SciFort solvers against SciPy."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

from compare_banded_scipy import build_cases, run

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ['src/kinds.f90','src/linalg.f90','src/linalg_banded.f90']


def compile_program(workdir: Path, sources: list[str], name: str, flags: list[str]) -> Path:
    exe = workdir / name
    subprocess.run(['gfortran','-ffree-line-length-none',*flags,
                    *[str(ROOT / f) for f in sources],'-o',str(exe)],
                   cwd=workdir,check=True)
    return exe


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    parser.add_argument('--cases',type=int,default=400)
    parser.add_argument('--max-order',type=int,default=70)
    parser.add_argument('--seed',type=int,default=109)
    parser.add_argument('--stress',action='store_true')
    parser.add_argument('--tolerance',type=float,default=3e-10)
    args=parser.parse_args()
    modes=['checked','optimized'] if args.mode=='both' else [args.mode]
    cases=build_cases(args.cases,args.seed,args.max_order,args.stress)
    with tempfile.TemporaryDirectory(prefix='scifort-linalg-banded-') as d:
        work=Path(d)
        for mode in modes:
            flags=['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
            native=compile_program(work,SOURCES+['test/test_linalg_banded.f90'],f'native_{mode}',flags)
            result=subprocess.run([str(native)],capture_output=True,text=True,check=True)
            print(mode+': '+result.stdout.strip())
            driver=compile_program(work,SOURCES+['tools/linalg_reference/driver_banded.f90'],f'driver_{mode}',flags)
            run(driver,cases,args.tolerance)


if __name__=='__main__':
    main()
