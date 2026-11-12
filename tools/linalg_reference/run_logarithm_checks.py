#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compile checked/optimized SciFort logm and fractional-power tests."""
import argparse
import subprocess
import tempfile
from pathlib import Path

from compare_logarithm_scipy import compare

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    'src/kinds.f90',
    'src/linalg.f90',
    'src/linalg_matrix_functions.f90',
    'src/linalg_expm_frechet.f90',
    'src/linalg_logarithm.f90',
]


def build(compiler, flags, directory, output, driver):
    cmd = [compiler, '-std=f2018', *flags, '-J', str(directory), '-I', str(directory),
           *[str(ROOT / name) for name in SOURCES], str(ROOT / driver),
           '-o', str(directory / output)]
    subprocess.run(cmd, check=True, cwd=directory)
    return str(directory / output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    parser.add_argument('--cases', type=int, default=240)
    parser.add_argument('--seed', type=int, default=20261107)
    parser.add_argument('--max-order', type=int, default=5)
    args = parser.parse_args()
    for mode in (('checked', 'optimized') if args.mode == 'both' else (args.mode,)):
        with tempfile.TemporaryDirectory(prefix=f'scifort_logm_{mode}_') as d:
            directory = Path(d)
            flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
            native = build(args.compiler, flags, directory, 'native', 'test/test_linalg_logarithm.f90')
            subprocess.run([native], check=True)
            driver = build(args.compiler, flags, directory, 'reference',
                           'tools/linalg_reference/driver_logarithm.f90')
            print(f'{mode}:', flush=True)
            compare(driver, args.cases, args.seed, max_order=args.max_order)


if __name__ == '__main__':
    main()
