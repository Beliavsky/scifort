#!/usr/bin/env python3
"""Compile/test SciFort's dense matrix functions, comparing with SciPy.

Usage:
    python tools/linalg_reference/run_matrix_functions_checks.py --mode both

Requires gfortran and Python with numpy/scipy. Uses a temporary build directory.
"""
import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from compare_matrix_functions_scipy import compare


ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    'src/kinds.f90',
    'src/linalg.f90',
    'src/linalg_matrix_functions.f90',
]


def compile_exe(compiler, flags, suffix, directory, driver):
    output = directory / suffix
    cmd = [compiler, '-std=f2018', *flags, '-J', str(directory),
           '-I', str(directory), *[str(ROOT / s) for s in SOURCES],
           str(ROOT / driver), '-o', str(output)]
    subprocess.run(cmd, check=True, cwd=directory)
    return str(output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    parser.add_argument('--cases', type=int, default=500,
                        help='randomized matrices per build')
    parser.add_argument('--seed', type=int, default=20261106)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error('cases must be positive')
    modes = ('checked', 'optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix=f'scifort_matrix_functions_{mode}_') as d:
            directory = Path(d)
            flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
            native = compile_exe(args.compiler, flags, 'native', directory,
                                 'test/test_linalg_matrix_functions.f90')
            subprocess.run([native], check=True)
            reference = compile_exe(args.compiler, flags, 'reference', directory,
                                    'tools/linalg_reference/matrix_functions_driver.f90')
            print(f'{mode}:', flush=True)
            compare(reference, args.cases, args.seed, False)


if __name__ == '__main__':
    main()
