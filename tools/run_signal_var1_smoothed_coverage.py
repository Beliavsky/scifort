#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build and validate exact VAR(1) tapered-population targets; optionally study coverage."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_signal_canonical_interval_methods import SOURCES

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT/'tools'/'signal_reference'


def run(optimized: bool, study: bool, trials: int) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-var1-taper-') as tmp:
        build = Path(tmp)
        flags = ['gfortran','-std=f2018','-fimplicit-none','-fbacktrace',
                 '-ffree-line-length-none','-J.','-I.']
        flags += ['-O2'] if optimized else ['-O0','-g','-fcheck=all']
        objects = []
        for name in ('kinds.f90','signal_var1_taper_spectrum.f90'):
            obj = build/(Path(name).stem+'.o')
            subprocess.run([*flags,'-c',str(ROOT/'src'/name),'-o',str(obj)],
                           cwd=build,check=True)
            objects.append(str(obj))
        exe = build/'test_signal_var1_taper_spectrum'
        subprocess.run([*flags,*objects,str(ROOT/'test'/f'{exe.name}.f90'),'-o',str(exe)],
                       cwd=build,check=True)
        subprocess.run([str(exe)],cwd=build,check=True)
        driver = build/'var1_taper_spectrum_driver'
        subprocess.run([*flags,*objects,str(REF/f'{driver.name}.f90'),'-o',str(driver)],
                       cwd=build,check=True)
        env = dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
        subprocess.run([sys.executable,str(REF/'verify_var1_taper_spectrum.py'),
                        str(driver),'--cases','80'],cwd=build,env=env,check=True)
        subprocess.run([sys.executable,str(REF/'verify_var1_taper_monte_carlo.py')],
                       cwd=build,env=env,check=True)
        if study:
            # The earlier BCa coverage driver consumes the actual windowed
            # spectral estimates. These modules are unchanged in this update.
            all_objects = []
            for name in SOURCES:
                obj = build/(Path(name).stem+'.o')
                if name != 'kinds.f90':
                    subprocess.run([*flags,'-c',str(ROOT/'src'/name),'-o',str(obj)],
                                   cwd=build,check=True)
                all_objects.append(str(obj))
            bc_driver = build/'canonical_bca_coverage_driver'
            subprocess.run([*flags,*all_objects,str(REF/f'{bc_driver.name}.f90'),
                            '-o',str(bc_driver)],cwd=build,check=True)
            subprocess.run([sys.executable,str(REF/'study_canonical_smoothed_var.py'),
                            str(bc_driver),'--trials',str(trials)],cwd=build,env=env,check=True)
        print('PASS '+('optimized' if optimized else 'checked')+
              ' VAR(1) tapered-spectrum numerical tests',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--checked-only',action='store_true')
    group.add_argument('--optimized-only',action='store_true')
    parser.add_argument('--study',action='store_true',
                        help='also run four-scenario canonical bootstrap coverage experiment')
    parser.add_argument('--trials',type=int,default=64)
    args = parser.parse_args()
    if not args.optimized_only:
        run(False,args.study,args.trials)
    if not args.checked_only:
        run(True,False,args.trials)
