#!/usr/bin/env python3
"""Compile/test dense Cholesky and compare with scipy.linalg in both modes."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
from compare_dense_cholesky_scipy import build_cases, run

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ['src/kinds.f90', 'src/linalg.f90', 'src/linalg_dense.f90']


def compile_fortran(work: Path, source: str, name: str, flags: list[str]) -> Path:
    exe = work / name
    cmd = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
           *[str(ROOT/s) for s in SOURCES], str(ROOT/source), '-o', str(exe)]
    subprocess.run(cmd, cwd=work, check=True, capture_output=True, text=True)
    return exe


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['checked','optimized','both'], default='both')
    parser.add_argument('--cases', type=int, default=360)
    parser.add_argument('--seed', type=int, default=389)
    parser.add_argument('--max-order', type=int, default=70)
    parser.add_argument('--tolerance', type=float, default=5e-9)
    parser.add_argument('--stress', action='store_true')
    args = parser.parse_args()
    if args.cases < 1 or args.max_order < 1 or args.tolerance <= 0:
        parser.error('cases/max-order/tolerance must be positive')
    cases = build_cases(args.cases, args.seed, args.max_order, args.stress)
    modes = ['checked', 'optimized'] if args.mode == 'both' else [args.mode]
    with tempfile.TemporaryDirectory(prefix='scifort-dense-cho-') as name:
        work = Path(name)
        for mode in modes:
            flags = (['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked'
                     else ['-O3'])
            native = compile_fortran(work, 'test/test_linalg_dense_cholesky.f90',
                                     'native_'+mode, flags)
            r = subprocess.run([str(native)], cwd=work, check=True,
                               capture_output=True, text=True)
            print(mode+': '+r.stdout.strip())
            driver = compile_fortran(work, 'tools/linalg_reference/driver_dense_cholesky.f90',
                                     'driver_'+mode, flags)
            print(mode+': ', end='', flush=True)
            run(driver,cases,args.tolerance)


if __name__=='__main__':
    main()
