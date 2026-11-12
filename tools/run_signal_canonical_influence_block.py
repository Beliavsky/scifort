#!/usr/bin/env python3
"""Checked and optimized statistic-specific canonical block pilot and old native suites."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_signal_canonical_moving_bootstrap import SOURCES

ROOT = Path(__file__).resolve().parents[1]
TESTS = (
    'test_signal_canonical_influence_block.f90',
    'test_signal_canonical_auto_block.f90',
    'test_signal_canonical_moving_bootstrap.f90',
    'test_signal_canonical_block_jackknife.f90',
    'test_signal_canonical_coherence.f90',
)


def run(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-canonical-influence-') as directory:
        build = Path(directory)
        flags = ['gfortran','-std=f2018','-fimplicit-none',
                 '-ffree-line-length-none','-fbacktrace']
        flags += ['-O2'] if optimized else ['-O0','-g','-fcheck=all']
        objects = []
        for name in SOURCES:
            path = ROOT / 'src' / name
            obj = build / (path.stem + '.o')
            subprocess.run([*flags,'-c',str(path),'-o',str(obj)],cwd=build,check=True)
            objects.append(str(obj))
        for name in TESTS:
            exe = build / Path(name).stem
            subprocess.run([*flags,*objects,str(ROOT/'test'/name),'-o',str(exe)],
                           cwd=build,check=True)
            subprocess.run([str(exe)],cwd=build,check=True)
        driver = build/'canonical_influence_block_driver'
        subprocess.run([*flags,*objects,str(ROOT/'tools'/'signal_reference'/
                                            'canonical_influence_block_driver.f90'),
                        '-o',str(driver)],cwd=build,check=True)
        env = dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
        subprocess.run([sys.executable,str(ROOT/'tools'/'signal_reference'/
                                            'verify_canonical_influence_block.py'),str(driver)],
                       cwd=build,env=env,check=True)
        if not optimized:
            var_driver = build/'canonical_influence_coverage_driver'
            subprocess.run([*flags,*objects,str(ROOT/'tools'/'signal_reference'/
                                                'canonical_influence_coverage_driver.f90'),
                            '-o',str(var_driver)],cwd=build,check=True)
            subprocess.run([sys.executable,str(ROOT/'tools'/'signal_reference'/
                                                'study_canonical_influence_var.py'),
                            str(var_driver),'--trials','48'],
                           cwd=build,env=env,check=True)
        print('PASS '+('optimized' if optimized else 'checked')+
              ' statistic-specific canonical block pilot and preceding native suites',flush=True)


if __name__ == '__main__':
    run(False)
    run(True)
