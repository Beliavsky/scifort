#!/usr/bin/env python3
"""Run spectral canonical-coefficient regressions in checked and optimized builds.

Requirements: GNU Fortran, Python, NumPy. No external Fortran BLAS or LAPACK.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_signal_canonical_coherence import SOURCES

ROOT = Path(__file__).resolve().parents[1]
TESTS = (
    'test_signal_canonical_vectors.f90',
    'test_signal_canonical_coherence.f90',
    'test_signal_group_coherence.f90',
    'test_signal_partial_shrinkage_cv.f90',
)


def run(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-canonical-vectors-') as tmp:
        build = Path(tmp)
        flags = ['gfortran', '-std=f2018', '-fimplicit-none',
                 '-ffree-line-length-none', '-fbacktrace']
        flags.extend(['-O3', '-march=native'] if optimized else ['-O0', '-g', '-fcheck=all'])
        objects = []
        for source in SOURCES:
            path = ROOT / 'src' / source
            obj = build / (path.stem+'.o')
            subprocess.run([*flags, '-c', str(path), '-o', str(obj)], cwd=build, check=True)
            objects.append(str(obj))
        for name in TESTS:
            executable = build / Path(name).stem
            subprocess.run([*flags, *objects, str(ROOT / 'test' / name),
                            '-o', str(executable)], cwd=build, check=True)
            subprocess.run([str(executable)], cwd=build, check=True)
        driver = build / 'canonical_vectors_driver'
        subprocess.run([*flags, *objects, str(ROOT / 'tools' / 'signal_reference' /
                                              'canonical_vectors_driver.f90'),
                        '-o', str(driver)], cwd=build, check=True)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, str(ROOT / 'tools' / 'signal_reference' /
                                            'verify_canonical_vectors.py'), str(driver)],
                       cwd=build, env=env, check=True)
        print(f'PASS {"optimized" if optimized else "checked"} canonical vectors suite',
              flush=True)


if __name__ == '__main__':
    run(False)
    run(True)
