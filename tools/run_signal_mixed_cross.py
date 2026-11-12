#!/usr/bin/env python3
"""Build checked and optimized mixed-type cross-spectrum regression suites.

Requirements: gfortran, Python with NumPy and SciPy. FPM is not required.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('kinds', 'signal_windows_extra', 'signal_dpss',
           'signal_spectral', 'signal_spectral_complex', 'signal_short_time_fft')
NATIVE = ('test_signal_mixed_cross', 'test_signal_short_time_fft',
          'test_signal_short_time_fft_helpers', 'test_signal_spectrogram_modes')
REF = 'signal_mixed_cross_reference_driver'


def run_variant(compiler: str, folder: Path, optimized: bool) -> None:
    """Build in an independent empty directory and run native/SciPy tests."""
    folder.mkdir(parents=True, exist_ok=True)
    flags = ['-std=f2018', '-ffree-line-length-none', '-fimplicit-none',
             '-Werror=implicit-interface', '-J.', '-I.']
    # Only the two changed Fortran modules need optimized compilation here.
    # The large DPSS dependency is compiled once at -O0 to keep clean builds
    # practical while still exercising optimized spectral/cross-STFT kernels.
    checked_flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace']
    objs = []
    for name in SOURCES:
        obj = f'{name}.o'
        mode = (['-O3'] if optimized and name in
                {'signal_spectral_complex','signal_short_time_fft'} else
                (['-O0'] if optimized else checked_flags))
        subprocess.run([compiler, *flags, *mode, '-c',
                        str(ROOT/'src'/f'{name}.f90'), '-o',obj],
                       cwd=folder,check=True)
        objs.append(obj)
    test_flags = ['-O3'] if optimized else checked_flags
    for name in NATIVE:
        subprocess.run([compiler,*flags,*test_flags,str(ROOT/'test'/f'{name}.f90'),*objs,'-o',name],
                       cwd=folder,check=True)
        subprocess.run([str(folder/name)],cwd=folder,check=True)
    subprocess.run([compiler,*flags,*test_flags,str(ROOT/'reference'/f'{REF}.f90'),*objs,'-o',REF],
                   cwd=folder,check=True)
    output = folder/'mixed_results.txt'
    with output.open('w') as stream:
        subprocess.run([str(folder/REF)],cwd=folder,stdout=stream,check=True)
    subprocess.run([sys.executable,str(ROOT/'tools'/'compare_signal_mixed_cross.py'),
                    str(output)],cwd=folder,check=True)
    print(f"PASS: {'targeted optimized -O3' if optimized else 'checked -fcheck=all'} "
          'mixed spectral tests and regressions',flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--checked-only',action='store_true')
    parser.add_argument('--optimized-only',action='store_true')
    opts = parser.parse_args()
    if opts.checked_only and opts.optimized_only:
        parser.error('Choose only one build mode')
    with tempfile.TemporaryDirectory(prefix='scifort_mixed_cross_') as scratch:
        folder = Path(scratch).resolve()
        if not opts.optimized_only:
            run_variant(opts.compiler,folder/'checked',False)
        if not opts.checked_only:
            run_variant(opts.compiler,folder/'optimized',True)


if __name__ == '__main__':
    main()
