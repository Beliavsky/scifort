#!/usr/bin/env python3
"""Test block-jackknife studentized canonical moving-bootstrap intervals."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_signal_canonical_interval_methods import SOURCES, TESTS as OLD_TESTS

ROOT = Path(__file__).resolve().parents[1]
TESTS = ('test_signal_canonical_studentized.f90',*OLD_TESTS)


def run(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-studentized-') as tmp:
        build = Path(tmp)
        flags = ['gfortran','-std=f2018','-fimplicit-none','-ffree-line-length-none',
                 '-fbacktrace','-J.','-I.']
        flags += ['-O2'] if optimized else ['-O0','-g','-fcheck=all']
        objects = []
        for source in SOURCES:
            path = ROOT/'src'/source
            obj = build/(path.stem+'.o')
            subprocess.run([*flags,'-c',str(path),'-o',str(obj)],cwd=build,check=True)
            objects.append(str(obj))
        for name in TESTS:
            exe = build/Path(name).stem
            subprocess.run([*flags,*objects,str(ROOT/'test'/name),'-o',str(exe)],
                           cwd=build,check=True)
            subprocess.run([str(exe)],cwd=build,check=True)
        reference_dir = ROOT/'tools'/'signal_reference'
        driver = build/'canonical_bootstrap_studentized_driver'
        subprocess.run([*flags,*objects,str(reference_dir/
            'canonical_bootstrap_studentized_driver.f90'),'-o',str(driver)],
            cwd=build,check=True)
        env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
                   MKL_NUM_THREADS='1')
        subprocess.run([sys.executable,str(reference_dir/'verify_canonical_studentized_bootstrap.py'),
                        str(driver),'--max-cases','16' if optimized else '30'],
                       cwd=build,env=env,check=True)
        # The previous NumPy reference is retained to verify exactly the
        # original Park-Miller draws and percentile/basic/BC calculations.
        old_driver = build/'canonical_bootstrap_intervals_driver'
        subprocess.run([*flags,*objects,str(reference_dir/'canonical_bootstrap_intervals_driver.f90'),
                        '-o',str(old_driver)],cwd=build,check=True)
        subprocess.run([sys.executable,str(reference_dir/'verify_canonical_bootstrap_intervals.py'),
                        str(old_driver),'--max-cases','2' if optimized else '8'],
                       cwd=build,env=env,check=True)
        if not optimized:
            study_driver = build/'canonical_studentized_coverage_driver'
            subprocess.run([*flags,*objects,str(reference_dir/
                'canonical_studentized_coverage_driver.f90'),'-o',str(study_driver)],
                cwd=build,check=True)
            subprocess.run([sys.executable,str(reference_dir/'study_canonical_studentized_var.py'),
                            str(study_driver),'--trials','32'],
                           cwd=build,env=env,check=True)
        print('PASS '+('optimized' if optimized else 'checked')+
              ' studentized and prior interval methods',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--checked-only',action='store_true')
    modes.add_argument('--optimized-only',action='store_true')
    args = parser.parse_args()
    if not args.optimized_only:
        run(False)
    if not args.checked_only:
        run(True)
