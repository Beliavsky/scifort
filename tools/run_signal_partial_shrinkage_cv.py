#!/usr/bin/env python3
"""Build checked and optimized signal CV shrinkage tests and SciPy/NumPy references.

Requires GNU Fortran, NumPy and SciPy; no third-party Fortran libraries,
BLAS/LAPACK or FPM needed. Use check_signal_partial_shrinkage_auto_public.py
for full top-level scifort_signal API compilation.
"""
from __future__ import annotations
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('kinds.f90', 'signal_dpss.f90', 'signal_windows_extra.f90',
           'signal_spectral.f90', 'signal_multitaper.f90', 'signal_partial_coherence.f90')
NATIVE = ('test_signal_partial_shrinkage_cv.f90',
          'test_signal_eigen_auto_shrinkage.f90',
          'test_signal_partial_shrinkage_auto.f90',
          'test_signal_partial_coherency.f90',
          'test_signal_tf_partial_coherence.f90',
          'test_signal_multitaper_tf_cross.f90',
          'test_signal_multitaper_multichannel.f90')


def run(cmd: list[str], build: Path, env: dict[str,str] | None = None) -> None:
    subprocess.run(cmd,cwd=build,env=env,check=True)


def validate(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-cv-shrink-') as temp:
        build = Path(temp)
        flags = ['gfortran', '-std=f2018', '-ffree-line-length-none',
                 '-fimplicit-none', '-fbacktrace']
        flags += ['-O3','-march=native'] if optimized else ['-O0','-g','-fcheck=all']
        objects = []
        for name in SOURCES:
            source = ROOT/'src'/name
            run([*flags,'-c',str(source)],build)
            objects.append(str(build/(source.stem+'.o')))
        for name in NATIVE:
            path = ROOT/'test'/name
            exe = build/path.stem
            run([*flags,*objects,str(path),'-o',str(exe)],build)
            run([str(exe)],build)
        source = ROOT/'tools'/'reference_signal_partial_shrinkage_cv.f90'
        exe = build/source.stem
        run([*flags,*objects,str(source),'-o',str(exe)],build)
        env = dict(os.environ,OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
        run([sys.executable,str(ROOT/'tools'/'compare_signal_partial_shrinkage_cv.py'),str(exe)],build,env)
        print(f'PASS {"optimized" if optimized else "checked"} CV shrinkage suite',flush=True)


if __name__ == '__main__':
    validate(False)
    validate(True)
