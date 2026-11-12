#!/usr/bin/env python3
"""Validate sliding multitaper partial coherence using gfortran, NumPy and SciPy.

The optimized pass recompiles the modified module and native tests at -O3,
using unchanged signal dependencies built with checks, to keep CI work small.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'kinds.f90',
    'signal_dpss.f90',
    'signal_windows_extra.f90',
    'signal_spectral.f90',
    'signal_multitaper.f90',
    'signal_partial_coherence.f90',
)
NATIVE = (
    'test_signal_tf_partial_coherence.f90',
    'test_signal_partial_coherence.f90',
    'test_signal_multitaper_tf_cross.f90',
    'test_signal_multitaper_time_frequency.f90',
    'test_signal_multitaper_multichannel.f90',
)


def run(args, cwd, env=None):
    subprocess.run(args, cwd=cwd, env=env, check=True)


def main():
    for name, flags in (
        ('checked', ['-O0', '-g', '-fcheck=all', '-fbacktrace']),
        ('targeted-optimized', ['-O3', '-march=native']),
    ):
        with tempfile.TemporaryDirectory(prefix='scifort-tf-partial-') as temp:
            work = Path(temp)
            base = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags]
            deps = ['gfortran', '-std=f2018', '-ffree-line-length-none',
                    '-O0', '-g', '-fcheck=all', '-fbacktrace']
            objects = []
            for source_name in SOURCES:
                source = ROOT / 'src' / source_name
                compile_flags = (base if name == 'checked' or
                                 source_name == 'signal_partial_coherence.f90'
                                 else deps)
                run([*compile_flags, '-c', str(source)], work)
                objects.append(str(work / (source.stem + '.o')))
            for native in NATIVE:
                binary = work / Path(native).stem
                run([*base, *objects, str(ROOT / 'test' / native),
                     '-o', str(binary)], work)
                run([str(binary)], work)
            driver = work / 'reference_signal_tf_partial_coherence'
            run([*base, *objects,
                 str(ROOT / 'tools/reference_signal_tf_partial_coherence.f90'),
                 '-o', str(driver)], work)
            env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                       MKL_NUM_THREADS='1')
            run([sys.executable,
                 str(ROOT / 'tools/compare_signal_tf_partial_coherence.py'),
                 str(driver)], work, env=env)
            print(f'PASS {name} build', flush=True)


if __name__ == '__main__':
    main()
