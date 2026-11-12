#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compile and run conditional multimode VAR(p) bootstrap coverage diagnostics.

Default is a small deterministic smoke test; --study runs a larger exploratory
simulation. These tests require NumPy/SciPy and GNU Fortran, not FPM.
"""
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


def run(optimized: bool, trials: int, nw: int, length: int,
        resamples: int, block_length: int, seed: int, show: bool,
        csv_path: Path | None = None) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-varp-coverage-') as tmp:
        build = Path(tmp)
        flags = ['gfortran','-std=f2018','-fimplicit-none','-fbacktrace',
                 '-ffree-line-length-none','-J.','-I.']
        flags += ['-O2'] if optimized else ['-O0','-g','-fcheck=all']
        objects=[]
        for name in SOURCES:
            target=build/(Path(name).stem+'.o')
            subprocess.run([*flags,'-c',str(ROOT/'src'/name),'-o',str(target)],
                           cwd=build,check=True)
            objects.append(str(target))
        driver=build/'canonical_varp_coverage_driver'
        subprocess.run([*flags,*objects,
                        str(REF/'canonical_varp_coverage_driver.f90'),
                        '-o',str(driver)],cwd=build,check=True)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
                 MKL_NUM_THREADS='1')
        subprocess.run([sys.executable,str(REF/'test_study_canonical_varp_conditional.py')],
                       cwd=build,env=env,check=True)
        command=[sys.executable,str(REF/'study_canonical_varp_conditional.py'),
                 str(driver),'--trials',str(trials),'--windows',str(nw),
                 '--samples-per-window',str(length),'--resamples',str(resamples),
                 '--block-length',str(block_length),'--seed',str(seed)]
        if csv_path is not None:
            command.extend(['--csv',str(csv_path.resolve())])
        if show:
            subprocess.run(command,cwd=build,env=env,check=True)
        else:
            completed=subprocess.run(command,cwd=build,env=env,check=True,
                                     capture_output=True,text=True)
            lines=completed.stdout.splitlines()
            print('\n'.join(lines[:4]))
        print('PASS '+('optimized' if optimized else 'checked')+
              ' conditional VAR(p) two-mode bootstrap study',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--checked-only',action='store_true')
    group.add_argument('--optimized-only',action='store_true')
    parser.add_argument('--study',action='store_true',
                        help='run more than the two-dataset-per-model smoke test')
    parser.add_argument('--trials',type=int,default=16,
                        help='datasets per VAR model when --study is used')
    parser.add_argument('--windows',type=int,default=40)
    parser.add_argument('--samples-per-window',type=int,default=16)
    parser.add_argument('--resamples',type=int,default=99)
    parser.add_argument('--block-length',type=int,default=4)
    parser.add_argument('--seed',type=int,default=20261030)
    parser.add_argument('--csv',type=Path,help='save study coverage, validity and point bias')
    args=parser.parse_args()
    # Run the heavier experiment only once; re-test the checked build in
    # optimized mode using the compact but complete method/configuration set.
    if not args.optimized_only:
        run(False,args.trials if args.study else 2,args.windows,
            args.samples_per_window,args.resamples if args.study else 39,
            args.block_length,args.seed,args.study,args.csv if args.study else None)
    if not args.checked_only:
        run(True,2,args.windows,args.samples_per_window,39,
            args.block_length,args.seed,False)
