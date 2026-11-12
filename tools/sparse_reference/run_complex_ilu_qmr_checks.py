#!/usr/bin/env python3
"""Compile and validate the complex ILU(0)/BiCG/QMR subset without FPM.

All required Fortran modules are compiled with the requested mode; the
full sparse regression suite remains a separate command.

Examples:
    python tools/sparse_reference/run_complex_ilu_qmr_checks.py --mode checked
    python tools/sparse_reference/run_complex_ilu_qmr_checks.py --mode optimized --per-solver 180
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
    'src/sparse_complex_ilu.f90',
    'src/sparse_linalg_complex.f90',
)


def check_mode(mode: str, per_solver: int, seeds: list[int]) -> None:
    with tempfile.TemporaryDirectory(prefix=f'scifort-complex-{mode}-') as directory:
        build = Path(directory)
        flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked'
                 else ['-O3'])
        cc = ['gfortran', '-std=f2018', *flags, '-J', str(build), '-I', str(build)]
        objects = []
        for relative in MODULES:
            obj = build / (Path(relative).stem + '.o')
            subprocess.run([*cc, '-c', str(ROOT / relative), '-o', str(obj)],
                           cwd=build, check=True)
            objects.append(str(obj))
        test = build / 'test_complex_ilu_qmr'
        subprocess.run([*cc, str(ROOT / 'test/test_sparse_complex_ilu_bicg_qmr.f90'),
                        *objects, '-o', str(test)], cwd=build, check=True)
        result = subprocess.run([str(test)], cwd=build, text=True,
                                capture_output=True, check=True)
        print(f'{mode}: {result.stdout.strip()}', flush=True)
        driver = build / 'complex_sparse_ilu_qmr_driver'
        subprocess.run([*cc, str(ROOT / 'tools/sparse_reference/complex_sparse_ilu_qmr_driver.f90'),
                        *objects, '-o', str(driver)], cwd=build, check=True)
        for seed in seeds:
            result = subprocess.run([sys.executable,
                                     str(ROOT / 'tools/sparse_reference/compare_sparse_complex_ilu_qmr_scipy.py'),
                                     str(driver), '--per-solver', str(per_solver),
                                     '--seed', str(seed)], text=True,
                                    capture_output=True, check=True)
            print(f'{mode}: {result.stdout.strip()}', flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    ap.add_argument('--per-solver', type=int, default=80)
    ap.add_argument('--seeds', type=int, nargs='+', default=[1728, 9894])
    args = ap.parse_args()
    if args.per_solver < 1:
        ap.error('per-solver must be positive')
    modes = ['checked', 'optimized'] if args.mode == 'both' else [args.mode]
    for mode in modes:
        check_mode(mode, args.per_solver, args.seeds)


if __name__ == '__main__':
    main()
