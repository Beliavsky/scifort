#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify block-aware BCa intervals, old interval APIs, and NumPy references."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_signal_canonical_interval_methods import SOURCES

ROOT = Path(__file__).resolve().parents[1]
TESTS = ('test_signal_canonical_bca.f90',
         'test_signal_canonical_studentized.f90',
         'test_signal_canonical_interval_methods.f90',
         'test_signal_canonical_moving_bootstrap.f90',
         'test_signal_canonical_block_jackknife.f90',
         'test_signal_canonical_jackknife.f90',
         'test_signal_canonical_coherence.f90')


def run(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-bca-') as directory:
        build = Path(directory)
        flags = ['gfortran','-std=f2018','-fimplicit-none',
                 '-ffree-line-length-none','-fbacktrace','-J.','-I.']
        flags += ['-O2'] if optimized else ['-O0','-g','-fcheck=all']
        objects = []
        for name in SOURCES:
            target = build/(Path(name).stem+'.o')
            subprocess.run([*flags,'-c',str(ROOT/'src'/name),'-o',str(target)],
                           cwd=build,check=True)
            objects.append(str(target))
        for name in TESTS:
            exe = build/Path(name).stem
            subprocess.run([*flags,*objects,str(ROOT/'test'/name),'-o',str(exe)],
                           cwd=build,check=True)
            subprocess.run([str(exe)],cwd=build,check=True)
        reference_dir = ROOT/'tools'/'signal_reference'
        env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
                   MKL_NUM_THREADS='1')
        for driver_file,script,limit in (
            ('canonical_bootstrap_bca_driver.f90',
             'verify_canonical_bca_bootstrap.py','16' if optimized else '30'),
            ('canonical_bootstrap_intervals_driver.f90',
             'verify_canonical_bootstrap_intervals.py','4' if optimized else '8')):
            driver = build/Path(driver_file).stem
            subprocess.run([*flags,*objects,str(reference_dir/driver_file),'-o',str(driver)],
                           cwd=build,check=True)
            subprocess.run([sys.executable,str(reference_dir/script),str(driver),
                            '--max-cases',limit],cwd=build,env=env,check=True)
        if not optimized:
            coverage_driver = build/'canonical_bca_coverage_driver'
            subprocess.run([*flags,*objects,str(reference_dir/
                'canonical_bca_coverage_driver.f90'),'-o',str(coverage_driver)],
                cwd=build,check=True)
            subprocess.run([sys.executable,str(reference_dir/'study_canonical_bca_var.py'),
                str(coverage_driver),'--trials','32'],cwd=build,env=env,check=True)
        print('PASS '+('optimized' if optimized else 'checked')+
              ' BCa/native/legacy regression tests',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--checked-only',action='store_true')
    mode.add_argument('--optimized-only',action='store_true')
    args = parser.parse_args()
    if not args.optimized_only:
        run(False)
    if not args.checked_only:
        run(True)
