#!/usr/bin/env python3
"""Build and independently verify SciFort FIR width/window additions.

    python tools/run_signal_fir_width_windows.py
    python tools/run_signal_fir_width_windows.py --checked-only

Requires GNU Fortran, Python, NumPy and SciPy; no FPM or network needed.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('kinds','signal_windows_extra','signal_dpss','signal_spectral',
           'signal_analog_design','signal_design','signal_fir_advanced')


def verify(compiler: str, dest: Path, optimized: bool) -> None:
    dest.mkdir()
    options = ['-std=f2018','-ffree-line-length-none','-fimplicit-none',
               '-Werror=implicit-interface','-I.','-J.']
    options += ['-O3'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    objects = []
    for name in SOURCES:
        obj = name+'.o'
        subprocess.run([compiler, *options, '-c', str(ROOT/'src'/f'{name}.f90'),
                        '-o',obj], cwd=dest, check=True)
        objects.append(obj)
    for name in ('test_signal_fir_width_windows', 'reference_signal_fir_width_windows'):
        file = (ROOT/'test' if name.startswith('test_') else ROOT/'tools')/f'{name}.f90'
        binary = dest/name
        subprocess.run([compiler,*options,str(file),*objects,'-o',str(binary)],
                       cwd=dest,check=True)
        if name.startswith('test_'):
            subprocess.run([str(binary)],cwd=dest,check=True)
        else:
            subprocess.run([sys.executable,str(ROOT/'tools'/'compare_signal_fir_width_windows.py'),
                            str(binary)],cwd=dest,check=True)
    print('PASS:', 'optimized -O3' if optimized else 'checked -fcheck=all',
          'FIR design suite',flush=True)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--checked-only',action='store_true')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='scifort_fir_width_') as tmp:
        verify(args.compiler, Path(tmp)/'checked',False)
        if not args.checked_only:
            verify(args.compiler, Path(tmp)/'optimized',True)


if __name__ == '__main__':
    main()
