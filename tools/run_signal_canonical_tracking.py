#!/usr/bin/env python3
"""Checked and optimized native tests and SciPy reference for canonical tracking.

Requires gfortran, Python, NumPy and SciPy. No external Fortran BLAS/LAPACK.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_signal_canonical_coherence import SOURCES as CANONICAL_SOURCES

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (*CANONICAL_SOURCES, 'signal_canonical_tracking.f90')
TESTS = (
    'test_signal_canonical_tracking.f90',
    'test_signal_canonical_vectors.f90',
    'test_signal_canonical_coherence.f90',
    'test_signal_group_coherence.f90',
)


def run(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-canonical-tracking-') as tmp:
        build = Path(tmp)
        flags = ['gfortran', '-std=f2018', '-fimplicit-none',
                 '-ffree-line-length-none', '-fbacktrace']
        flags.extend(['-O3', '-march=native'] if optimized else ['-O0', '-g', '-fcheck=all'])
        objects = []
        for source in SOURCES:
            path = ROOT / 'src' / source
            obj = build / (path.stem + '.o')
            subprocess.run([*flags, '-c', str(path), '-o', str(obj)], cwd=build, check=True)
            objects.append(str(obj))
        for name in TESTS:
            executable = build / Path(name).stem
            subprocess.run([*flags, *objects, str(ROOT / 'test' / name),
                            '-o', str(executable)], cwd=build, check=True)
            subprocess.run([str(executable)], cwd=build, check=True)
        driver = build / 'canonical_tracking_driver'
        subprocess.run([*flags, *objects, str(ROOT / 'tools' / 'signal_reference' /
                                              'canonical_tracking_driver.f90'),
                        '-o', str(driver)], cwd=build, check=True)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, str(ROOT / 'tools' / 'signal_reference' /
                                            'verify_canonical_tracking.py'), str(driver)],
                       cwd=build, env=env, check=True)
        print(f'PASS {"optimized" if optimized else "checked"} canonical tracking suite',
              flush=True)


if __name__ == '__main__':
    run(False)
    run(True)
