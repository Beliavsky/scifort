#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Clean public signal closure; test VAR(p) spectral and canonical APIs."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from check_signal_canonical_bca_public import TESTS as OLD_TESTS
from check_signal_mixed_public import closure

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'tools'/'signal_reference'


def check(build: Path, cases: int):
    order=closure()
    build.mkdir(parents=True,exist_ok=True)
    flags=['gfortran','-std=f2018','-fimplicit-none','-ffree-line-length-none',
           '-O0','-J.','-I.']
    objects=[]
    print(f'VAR(p) public signal dependency closure: {len(order)} modules',flush=True)
    for k,source in enumerate(order,1):
        obj=build/(source.stem+'.o')
        subprocess.run([*flags,'-c',str(source),'-o',str(obj)],cwd=build,check=True)
        objects.append(str(obj))
        if k%20==0 or k==len(order):
            print(f'Compiled {k}/{len(order)} public modules',flush=True)
    for name in ('test_signal_var1_taper_spectrum_public.f90',
                 'test_signal_varp_population_public.f90',
                 'test_signal_varp_multitaper_public.f90',
                 'test_signal_fixed_taper_spectrogram_public.f90',*OLD_TESTS):
        exe=build/Path(name).stem
        subprocess.run([*flags,str(ROOT/'test'/name),*objects,'-o',str(exe)],
                       cwd=build,check=True)
        subprocess.run([str(exe)],cwd=build,check=True)
    print(f'PASS {len(OLD_TESTS)+4} public signal tests',flush=True)

    driver=build/'varp_canonical_target_driver'
    subprocess.run([*flags,*objects,str(REF/(driver.name+'.f90')),'-o',str(driver)],
                   cwd=build,check=True)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
    subprocess.run([sys.executable,str(REF/'verify_varp_canonical_target.py'),
                    str(driver),'--cases',str(cases)],cwd=build,env=env,check=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir',type=Path)
    parser.add_argument('--cases',type=int,default=50)
    args=parser.parse_args()
    if args.build_dir:
        check(args.build_dir.resolve(),args.cases)
    else:
        with tempfile.TemporaryDirectory(prefix='scifort-varp-public-') as tmp:
            check(Path(tmp),args.cases)
