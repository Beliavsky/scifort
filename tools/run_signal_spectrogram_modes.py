#!/usr/bin/env python3
"""Validate SciFort legacy spectrogram output modes against SciPy (checked and -O3)."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ['kinds', 'signal_windows_extra', 'signal_dpss', 'signal_spectral',
           'signal_spectral_complex', 'signal_spectral_nd', 'signal_timefreq_nd']


def check(compiler: str, path: Path, optimized: bool):
    path.mkdir()
    flags = ['-std=f2018','-fimplicit-none','-Werror=implicit-interface',
             '-ffree-line-length-none','-I.','-J.']
    flags += ['-O3'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    objects = []
    for source in SOURCES:
        object_name = source+'.o'
        subprocess.run([compiler,*flags,'-c',str(ROOT/'src'/f'{source}.f90'),'-o',object_name],
                       cwd=path,check=True)
        objects.append(object_name)
    for source,folder in [('test_signal_spectrogram_modes','test'),
                          ('test_signal_spectrogram_modes_nd','test'),
                          ('reference_signal_spectrogram_modes','tools')]:
        exe = path/source
        subprocess.run([compiler,*flags,str(ROOT/folder/f'{source}.f90'),*objects,'-o',str(exe)],
                       cwd=path,check=True)
        if folder == 'test':
            subprocess.run([str(exe)],cwd=path,check=True)
        else:
            subprocess.run([sys.executable,str(ROOT/'tools'/'compare_signal_spectrogram_modes.py'),
                            str(exe)],cwd=path,check=True)
    print('PASS:', '-O3' if optimized else '-fcheck=all', 'spectrogram modes', flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--checked-only',action='store_true')
    args=p.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_spectrogram_modes_') as path:
        check(args.compiler,Path(path)/'checked',False)
        if not args.checked_only:
            check(args.compiler,Path(path)/'optimized',True)


if __name__ == '__main__':
    main()
