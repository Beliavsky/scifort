#!/usr/bin/env python3
"""Compile and compare complex pivoted sparse LU against SciPy without FPM.

Examples:
    python tools/sparse_reference/run_complex_direct_checks.py --mode both
    python tools/sparse_reference/run_complex_direct_checks.py --mode optimized --count 500
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
MODULES = (
    'src/kinds.f90',
    'src/random.f90',
    'src/sparse.f90',
    'src/sparse_coo.f90',
    'src/sparse_coo_nd.f90',
    'src/sparse_coo_nd_complex.f90',
    'src/sparse_complex.f90',
    'src/sparse_linalg_complex_direct.f90',
)


def check(mode: str, count: int, seeds: list[int]) -> None:
    with tempfile.TemporaryDirectory(prefix=f'scifort-sparse-direct-{mode}-') as directory:
        build = Path(directory)
        flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
        fc = ['gfortran','-std=f2018',*flags,'-J',str(build),'-I',str(build)]
        objs = []
        for file in MODULES:
            dest = build / (Path(file).stem+'.o')
            subprocess.run([*fc,'-c',str(ROOT / file),'-o',str(dest)],
                           cwd=build,check=True)
            objs.append(str(dest))
        exe = build/'test_sparse_complex_direct'
        subprocess.run([*fc,str(ROOT/'test/test_sparse_complex_direct.f90'),
                        *objs,'-o',str(exe)],cwd=build,check=True)
        output = subprocess.run([str(exe)],capture_output=True,text=True,check=True)
        print(f'{mode}: {output.stdout.strip()}',flush=True)
        driver = build/'complex_sparse_direct_driver'
        subprocess.run([*fc,str(ROOT/'tools/sparse_reference/complex_sparse_direct_driver.f90'),
                        *objs,'-o',str(driver)],cwd=build,check=True)
        for seed in seeds:
            comparison = subprocess.run([
                sys.executable,
                str(ROOT/'tools/sparse_reference/compare_sparse_complex_direct_scipy.py'),
                str(driver),'--count',str(count),'--seed',str(seed)],
                capture_output=True,text=True,check=True)
            print(f'{mode}, seed={seed}: {comparison.stdout.strip()}',flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--count',type=int,default=350)
    parser.add_argument('--seeds',type=int,nargs='+',default=[2667,4113])
    args = parser.parse_args()
    if args.count < 1:
        parser.error('count must be positive')
    for mode in (['checked','optimized'] if args.mode == 'both' else [args.mode]):
        check(mode,args.count,args.seeds)


if __name__ == '__main__':
    main()
