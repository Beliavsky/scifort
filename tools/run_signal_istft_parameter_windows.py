#!/usr/bin/env python3
"""Build and compare parameterized ISTFT windows in both checked and O3 configurations.

    python tools/run_signal_istft_parameter_windows.py
    python tools/run_signal_istft_parameter_windows.py --checked-only

Requires Python, NumPy, SciPy, and a Fortran 2018 compiler; no fpm needed.
Builds out of tree and does not download Fortran dependencies.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'kinds', 'signal_windows_extra', 'signal_dpss',
    'signal_spectral', 'signal_spectral_complex',
    'signal_spectral_nd', 'signal_timefreq_nd',
    'signal_istft', 'signal_istft_nd',
)


def run(command, directory):
    subprocess.run(command, cwd=directory, check=True)


def check(compiler, build, optimize):
    build.mkdir(parents=True)
    flags = ['-std=f2018', '-ffree-line-length-none', '-fimplicit-none',
             '-Werror=implicit-interface', '-J.', '-I.']
    flags.extend(['-O3'] if optimize else ['-O0', '-g', '-fcheck=all', '-fbacktrace'])
    objects = []
    for source in SOURCES:
        object_name = source + '.o'
        run([compiler, *flags, '-c', str(ROOT/'src'/ (source+'.f90')),
             '-o', object_name],build)
        objects.append(object_name)
    exeext = '.exe' if os.name=='nt' else ''
    for name in ('test_signal_istft_parameter_windows', 'test_signal_istft',
                 'test_signal_windowed_estimators', 'test_signal_stft_boundaries'):
        exe = build/(name+exeext)
        run([compiler,*flags,str(ROOT/'test'/(name+'.f90')), *objects,'-o',str(exe)],build)
        run([str(exe)],build)
    reference=build/('reference_signal_istft_parameter_windows'+exeext)
    run([compiler,*flags,str(ROOT/'tools'/'reference_signal_istft_parameter_windows.f90'),
         *objects,'-o',str(reference)],build)
    run([sys.executable,str(ROOT/'tools'/'compare_signal_istft_parameter_windows.py'),
         str(reference)],build)
    print('PASS:', 'O3' if optimize else 'checked', 'ISTFT window test slice',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--compiler',default='gfortran')
    ap.add_argument('--checked-only',action='store_true')
    options=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_istft_windows_') as temp:
        check(options.compiler,Path(temp)/'checked',False)
        if not options.checked_only:
            check(options.compiler,Path(temp)/'optimized',True)


if __name__=='__main__':
    main()
