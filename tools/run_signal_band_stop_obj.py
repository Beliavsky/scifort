#!/usr/bin/env python3
"""Build and test SciFort scipy.signal.band_stop_obj without FPM/BLAS.

Usage: python tools/run_signal_band_stop_obj.py [--checked-only]
Requires gfortran, NumPy, SciPy.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('kinds','signal_elliptic','signal_band_stop_obj')

def run(build: Path, optimized: bool) -> None:
    build.mkdir()
    flags=['-std=f2018','-fimplicit-none','-Werror=implicit-interface',
           '-ffree-line-length-none','-I.','-J.']
    flags+=['-O3'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    objects=[]
    for stem in SOURCES:
        obj=stem+'.o'
        subprocess.run(['gfortran',*flags,'-c',str(ROOT/'src'/(stem+'.f90')),
                        '-o',obj],cwd=build,check=True)
        objects.append(obj)
    for stem, folder in [('test_signal_band_stop_obj','test'),
                         ('band_stop_obj_driver','tools/signal_reference')]:
        exe=build/stem
        subprocess.run(['gfortran',*flags,str(ROOT/folder/(stem+'.f90')),*objects,
                        '-o',str(exe)],cwd=build,check=True)
        if folder=='test':
            subprocess.run([str(exe)],cwd=build,check=True)
        else:
            subprocess.run([sys.executable,str(ROOT/'tools/signal_reference/verify_band_stop_obj.py'),
                            str(exe)],cwd=build,check=True)
    print('PASS: band-stop objective', 'optimized -O3' if optimized else 'checked -fcheck=all')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checked-only',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_bandstop_') as temp:
        run(Path(temp)/'checked',False)
        if not args.checked_only:
            run(Path(temp)/'optimized',True)
