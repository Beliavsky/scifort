#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Checked and optimized Fortran and independent SciPy/NumPy VAR(p) tests."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'tools'/'signal_reference'


def run(optimized: bool, cases: int):
    with tempfile.TemporaryDirectory(prefix='scifort-varp-') as tmp:
        b=Path(tmp)
        flags=['gfortran','-std=f2018','-fimplicit-none','-ffree-line-length-none',
               '-fbacktrace','-J.','-I.']
        flags+=['-O2'] if optimized else ['-O0','-g','-fcheck=all']
        objs=[]
        for source in ('kinds.f90','signal_var1_taper_spectrum.f90',
                       'signal_varp_taper_spectrum.f90'):
            obj=b/(Path(source).stem+'.o')
            subprocess.run([*flags,'-c',str(ROOT/'src'/source),'-o',str(obj)],
                           cwd=b,check=True)
            objs.append(str(obj))
        for name in ('test_signal_varp_taper_spectrum','varp_taper_spectrum_driver'):
            path=(ROOT/'test' if name.startswith('test_') else REF)/(name+'.f90')
            exe=b/name
            subprocess.run([*flags,*objs,str(path),'-o',str(exe)],cwd=b,check=True)
            if name.startswith('test_'):
                subprocess.run([str(exe)],cwd=b,check=True)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
        subprocess.run([sys.executable,str(REF/'verify_varp_taper_spectrum.py'),
                        str(b/'varp_taper_spectrum_driver'),'--cases',str(cases)],
                       cwd=b,env=env,check=True)
        subprocess.run([sys.executable,str(REF/'verify_varp_taper_monte_carlo.py'),
                        '--trials','20000'],cwd=b,env=env,check=True)
        print('PASS '+('optimized' if optimized else 'checked')+' VAR(p) spectral references',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--checked-only',action='store_true')
    group.add_argument('--optimized-only',action='store_true')
    parser.add_argument('--cases',type=int,default=80)
    args=parser.parse_args()
    if not args.optimized_only:
        run(False,args.cases)
    if not args.checked_only:
        run(True,args.cases)
