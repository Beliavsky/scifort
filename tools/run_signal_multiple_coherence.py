#!/usr/bin/env python3
"""Checked/optimized GNU Fortran regressions plus independent NumPy validation.

No FPM or third-party Fortran dependencies required. To also verify the
scifort_signal public API, run check_signal_multiple_coherence_public.py.
"""
from __future__ import annotations
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('kinds.f90','signal_dpss.f90','signal_windows_extra.f90',
           'signal_spectral.f90','signal_multitaper.f90',
           'signal_partial_coherence.f90','signal_multiple_coherence.f90')
TESTS = ('test_signal_multiple_coherence.f90','test_signal_partial_shrinkage_cv.f90',
         'test_signal_eigen_auto_shrinkage.f90',
         'test_signal_partial_coherence.f90','test_signal_tf_partial_coherence.f90',
         'test_signal_partial_shrinkage_auto.f90')

def run(cmd: list[str], folder: Path) -> None:
    subprocess.run(cmd,cwd=folder,check=True)


def validate(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-multiple-coherence-') as tmp:
        path = Path(tmp)
        flags = ['gfortran','-std=f2018','-fimplicit-none','-ffree-line-length-none',
                 '-fbacktrace']
        flags += ['-O3','-march=native'] if optimized else ['-O0','-g','-fcheck=all']
        objects=[]
        for source in SOURCES:
            run([*flags,'-c',str(ROOT/'src'/source)],path)
            objects.append(str(path/(Path(source).stem+'.o')))
        for name in TESTS:
            exe = path/Path(name).stem
            run([*flags,*objects,str(ROOT/'test'/name),'-o',str(exe)],path)
            run([str(exe)],path)
        driver=path/'reference_multiple_coherence'
        run([*flags,*objects,str(ROOT/'tools'/'reference_signal_multiple_coherence.f90'),
             '-o',str(driver)],path)
        env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        subprocess.run([sys.executable,str(ROOT/'tools'/'compare_signal_multiple_coherence.py'),
                        str(driver)],cwd=path,env=env,check=True)
        print(f'PASS {"optimized" if optimized else "checked"} multiple coherence suite',flush=True)


if __name__ == '__main__':
    validate(False)
    validate(True)
