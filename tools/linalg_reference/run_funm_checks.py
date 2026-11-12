#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build and validate scipy.linalg.funm compatibility in checked/O3 modes."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
from compare_funm_scipy import compare
from compare_funm_stress import compare as stress
from compare_funm_reordered_scipy import compare as reordered
from compare_funm_diagonalizable_scipy import compare as diagonalizable
from compare_funm_taylor_scipy import compare as taylor

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ['src/kinds.f90', 'src/linalg.f90',
           'src/linalg_matrix_functions.f90', 'src/linalg_funm.f90']


def compile_driver(compiler, flags, folder, filename, output):
    destination = folder/output
    cmd = [compiler,'-std=f2018',*flags,'-J',str(folder),'-I',str(folder),
           *map(lambda x:str(ROOT/x), SOURCES),str(ROOT/filename),
           '-o',str(destination)]
    subprocess.run(cmd,cwd=folder,check=True)
    return str(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=400)
    parser.add_argument('--seed',type=int,default=20261107)
    parser.add_argument('--max-order',type=int,default=8)
    parser.add_argument('--skip-stress',action='store_true')
    parser.add_argument('--reorder-cases',type=int,default=120)
    parser.add_argument('--diagonalizable-cases',type=int,default=100)
    parser.add_argument('--taylor-cases',type=int,default=180)
    args=parser.parse_args()
    if (args.cases < 1 or args.reorder_cases < 1 or args.diagonalizable_cases < 1 or
            args.taylor_cases < 1 or not (1<=args.max_order<=40)):
        parser.error('cases and reorder-cases must be positive; max-order must be 1..40')
    modes = ('checked','optimized') if args.mode=='both' else (args.mode,)
    with tempfile.TemporaryDirectory(prefix='scifort-funm-') as temp:
        for mode in modes:
            folder = Path(temp)/mode
            folder.mkdir()
            flags=(['-O0','-g','-fcheck=all','-fbacktrace']
                   if mode=='checked' else ['-O3'])
            native=compile_driver(args.compiler,flags,folder,
                                  'test/test_linalg_funm.f90','native')
            subprocess.run([native],check=True)
            example=compile_driver(args.compiler,flags,folder,
                                   'example/linalg_funm_example.f90','example')
            subprocess.run([example],check=True)
            driver=compile_driver(args.compiler,flags,folder,
                                  'tools/linalg_reference/driver_funm.f90','driver')
            print(mode+':',flush=True)
            compare(driver, args.cases, args.seed, args.max_order)
            if not args.skip_stress:
                stress(driver)
                reordered(driver, args.reorder_cases, args.seed)
                diagonalizable(driver, args.diagonalizable_cases, args.seed)
            taylor(driver, args.taylor_cases, args.seed)


if __name__=='__main__':
    main()
