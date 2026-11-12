#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build/test exact fixed-taper multichannel spectra and VAR(p) DPSS coverage.

Requirements: GNU Fortran, NumPy, and SciPy. Does not require FPM or BLAS.
Default runs checked and optimized native and independent Python tests; use
--study with --csv to generate exploratory coverage tables.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'tools'/'signal_reference'
MODULE=re.compile(r'^\s*module\s+(scifort_\w+)',re.I|re.M)
USE=re.compile(r'^\s*use\s+(scifort_\w+)',re.I|re.M)


def ordered_dependencies():
    modules={}
    for file in (ROOT/'src').rglob('*.f90'):
        for mod in MODULE.findall(file.read_text()):
            modules[mod.lower()]=file
    visited=set()
    result=[]
    def include(name):
        if name in visited:
            return
        visited.add(name)
        source=modules[name]
        for dependency in USE.findall(source.read_text()):
            include(dependency.lower())
        if source not in result:
            result.append(source)
    for name in ('scifort_signal_canonical_bootstrap',
                 'scifort_signal_fixed_taper_spectrogram'):
        include(name)
    return result


def run(optimized, *, cases, trials, nw, length, resamples, blocks, settings,seed,
        show, csv_path):
    with tempfile.TemporaryDirectory(prefix='scifort-fixed-dpss-') as temp:
        build=Path(temp)
        flags=['gfortran','-std=f2018','-fimplicit-none',
               '-ffree-line-length-none','-fbacktrace','-J.','-I.']
        flags+=['-O2'] if optimized else ['-O0','-g','-fcheck=all']
        objects=[]
        for source in ordered_dependencies():
            output=build/(source.stem+'.o')
            subprocess.run([*flags,'-c',str(source),'-o',str(output)],cwd=build,check=True)
            objects.append(str(output))
        def executable(name, location):
            out=build/name
            subprocess.run([*flags,str(location/(name+'.f90')),*objects,'-o',str(out)],
                           cwd=build,check=True)
            return out
        test=executable('test_signal_fixed_taper_spectrogram',ROOT/'test')
        subprocess.run([str(test)],cwd=build,check=True)
        test_driver=executable('fixed_taper_spectrogram_driver',REF)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
        subprocess.run([sys.executable,str(REF/'verify_fixed_taper_spectrogram.py'),
                        str(test_driver),'--cases',str(cases)],cwd=build,env=env,check=True)
        subprocess.run([sys.executable,str(REF/'test_study_canonical_varp_multitaper.py')],
                       cwd=build,env=env,check=True)
        study_driver=executable('canonical_varp_multitaper_coverage_driver',REF)
        command=[sys.executable,str(REF/'study_canonical_varp_multitaper.py'),
                 str(study_driver),'--trials',str(trials),'--windows',str(nw),
                 '--samples-per-window',str(length),'--resamples',str(resamples),
                 '--block-lengths',blocks,'--settings',settings,'--seed',str(seed)]
        if csv_path is not None:
            command+=['--csv',str(csv_path.resolve())]
        if show:
            subprocess.run(command,cwd=build,env=env,check=True)
        else:
            completed=subprocess.run(command,cwd=build,env=env,check=True,
                                     capture_output=True,text=True)
            print('\n'.join(completed.stdout.splitlines()[:4]),flush=True)
        print('PASS '+('optimized' if optimized else 'checked')+
              ' native exact fixed-taper/VAR(p) DPSS bootstrap',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--checked-only',action='store_true')
    mode.add_argument('--optimized-only',action='store_true')
    parser.add_argument('--study',action='store_true',help='exploratory grid instead of smoke run')
    parser.add_argument('--trials',type=int,default=8)
    parser.add_argument('--windows',type=int,default=40)
    parser.add_argument('--samples-per-window',type=int,default=16)
    parser.add_argument('--resamples',type=int,default=79)
    parser.add_argument('--block-lengths',default='2,4,8')
    parser.add_argument('--settings',default='2.5:1:equal,2.5:3:eigen,3.5:5:equal,3.5:5:eigen')
    parser.add_argument('--seed',type=int,default=20261101)
    parser.add_argument('--cases',type=int,default=50)
    parser.add_argument('--csv',type=Path,help='write per-scalar exploratory coverage statistics')
    args=parser.parse_args()
    smoke_settings='2.5:1:equal,3.5:5:eigen'
    if not args.optimized_only:
        run(False,cases=args.cases,trials=args.trials if args.study else 1,
            nw=args.windows,length=args.samples_per_window,
            resamples=args.resamples if args.study else 27,
            blocks=args.block_lengths if args.study else '2,4',
            settings=args.settings if args.study else smoke_settings,
            seed=args.seed,show=args.study,csv_path=args.csv if args.study else None)
    if not args.checked_only:
        run(True,cases=args.cases,trials=1,nw=args.windows,
            length=args.samples_per_window,resamples=27,blocks='2,4',
            settings=smoke_settings,seed=args.seed,show=False,csv_path=None)
