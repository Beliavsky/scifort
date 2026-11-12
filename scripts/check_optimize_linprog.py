#!/usr/bin/env python3
"""Compile and test the standalone two-phase SciFort linprog implementation.

Requirements: gfortran, Python with NumPy and SciPy. No FPM needed. This
checks the new module only, not the complete SciFort dependency graph.
"""
from pathlib import Path
import os
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent
GFORTRAN = os.environ.get('FC', 'gfortran')
SOURCES = [ROOT / 'src/kinds.f90', ROOT / 'src/optimize/linprog.f90']


def run(command, cwd=None):
    subprocess.run([str(x) for x in command], cwd=cwd, check=True)


def build(mode, flags):
    with tempfile.TemporaryDirectory(prefix=f'scifort_linprog_{mode}_') as t:
        work = Path(t)
        run([GFORTRAN, '-std=f2008','-Wall','-Wextra', *flags,
             '-c', *SOURCES], cwd=work)
        objects = [work / 'kinds.o', work / 'linprog.o']
        suffix = '.exe' if os.name == 'nt' else ''
        native = work / ('test_linprog' + suffix)
        driver = work / ('linprog_driver' + suffix)
        example = work / ('linprog_example' + suffix)
        for src, dest in [
            (ROOT / 'test/test_optimize_linprog.f90', native),
            (ROOT / 'tools/optimize_reference/linprog_reference_driver.f90', driver),
            (ROOT / 'example/optimize_linprog_example.f90', example),
        ]:
            run([GFORTRAN, '-std=f2008', *flags, '-I',work, src, *objects, '-o',dest],cwd=work)
        print(f'=== {mode}: native tests ===',flush=True)
        run([native])
        print(f'=== {mode}: example ===',flush=True)
        run([example])
        print(f'=== {mode}: SciPy reference comparisons ===',flush=True)
        run([sys.executable, ROOT / 'test/compare_linprog_scipy.py', driver])


if __name__ == '__main__':
    build('checked',['-O0','-g','-fcheck=all','-fbacktrace'])
    build('optimized',['-O3'])
