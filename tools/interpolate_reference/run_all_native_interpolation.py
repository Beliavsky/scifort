#!/usr/bin/env python3
"""Build the public interpolation facade and run all native interpolation tests.

Uses GNU Fortran directly: no internet connection, FPM, or third-party root
backend is required. Build objects once and link each test against them.
"""
from __future__ import annotations
import argparse
import subprocess
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SOURCES=('kinds','interpolate_piecewise','interpolate_bspline','interpolate_fitpack_bivariate','interpolate_fitpack_compat','interpolate_smoothing',
         'interpolate_regular_grid','interpolate_barycentric',
         'interpolate_global','interpolate_kdtree','interpolate_rbf','interpolate_legacy',
         'interpolate_scattered','interpolate_clough_tocher','interpolate')


def run(cmd, cwd):
    result=subprocess.run(cmd,cwd=cwd,text=True,capture_output=True)
    if result.returncode:
        raise RuntimeError('Failed command: '+ ' '.join(map(str,cmd)) +
                           '\nstdout:\n'+result.stdout[-2500:]+
                           '\nstderr:\n'+result.stderr[-3000:])
    return result.stdout.strip()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--only',choices=('checked','optimized','both'),default='both')
    args=ap.parse_args()
    tests=sorted((ROOT/'test').glob('test_interpolate*.f90'))
    modes=[]
    if args.only!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.only!='checked':modes.append(('optimized',['-O3']))
    for title,flags in modes:
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            common=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
            for name in SOURCES:
                run([*common,'-I'+str(d),'-J'+str(d),'-c',str(ROOT/'src'/f'{name}.f90'),
                     '-o',str(d/f'{name}.o')],d)
            objs=[str(d/f'{name}.o') for name in SOURCES]
            for test in tests:
                exe=d/'native_test'
                run([*common,'-I'+str(d),str(test),*objs,'-o',str(exe)],d)
                result=run([str(exe)],d)
                print(f'{title}: {test.name}: {result}',flush=True)
        print(f'{title}: {len(tests)}/{len(tests)} native interpolation programs PASS',flush=True)

if __name__=='__main__': main()
