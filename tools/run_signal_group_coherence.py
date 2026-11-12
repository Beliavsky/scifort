#!/usr/bin/env python3
"""Rebuild and independently validate conditional group coherence."""
from __future__ import annotations
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'kinds.f90','signal_dpss.f90','signal_windows_extra.f90',
    'signal_spectral.f90','signal_multitaper.f90',
    'signal_partial_coherence.f90','signal_multiple_coherence.f90',
    'signal_band_coherence.f90','signal_group_coherence.f90',
)
TESTS = ('test_signal_group_coherence.f90', 'test_signal_band_coherence.f90',
         'test_signal_multiple_coherence.f90', 'test_signal_partial_coherence.f90',
         'test_signal_partial_shrinkage_auto.f90', 'test_signal_partial_shrinkage_cv.f90')

def run(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-group-coherence-') as temporary:
        folder = Path(temporary)
        flags = ['gfortran','-std=f2018','-fimplicit-none','-ffree-line-length-none','-fbacktrace']
        flags += ['-O3','-march=native'] if optimized else ['-O0','-g','-fcheck=all']
        objects=[]
        for name in SOURCES:
            source=ROOT/'src'/name
            subprocess.run([*flags,'-c',str(source)],cwd=folder,check=True)
            objects.append(str(folder/(source.stem+'.o')))
        for name in TESTS:
            program=folder/Path(name).stem
            subprocess.run([*flags,*objects,str(ROOT/'test'/name),'-o',str(program)],
                           cwd=folder,check=True)
            subprocess.run([str(program)],cwd=folder,check=True)
        driver=folder/'group_coherence_driver'
        subprocess.run([*flags,*objects,str(ROOT/'tools'/'signal_reference'/
                          'group_coherence_driver.f90'),'-o',str(driver)],cwd=folder,check=True)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        subprocess.run([sys.executable,str(ROOT/'tools'/'signal_reference'/
                        'verify_group_coherence.py'),str(driver)],
                       cwd=folder,env=env,check=True)
        print('PASS '+('optimized' if optimized else 'checked')+' group-coherence suite',flush=True)

if __name__ == '__main__':
    run(False)
    run(True)
