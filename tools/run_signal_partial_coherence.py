#!/usr/bin/env python3
"""Build and validate conditional multitaper coherence against NumPy/SciPy.

Requires GNU Fortran, NumPy and SciPy. No FPM or external Fortran libraries.
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
NATIVE_TESTS = (
    'test_signal_partial_coherence.f90',
    'test_signal_multitaper_multichannel.f90',
    'test_signal_dpss_multitaper.f90',
    'test_signal_multitaper_adaptive.f90',
)


def call(args, cwd, env=None):
    subprocess.run(args, cwd=cwd, env=env, check=True)


def main():
    for label, flags in (
        ('checked', ['-O0', '-g', '-fcheck=all', '-fbacktrace']),
        ('targeted-optimized', ['-O3', '-march=native']),
    ):
        with tempfile.TemporaryDirectory(prefix='scifort-partial-') as td:
            directory = Path(td)
            base = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags]
            dependencies = ['gfortran', '-std=f2018', '-ffree-line-length-none',
                            '-O0', '-fcheck=all']
            objects = []
            for filename in SOURCES:
                source = ROOT / 'src' / filename
                compiler = base if label == 'checked' or filename == 'signal_partial_coherence.f90' else dependencies
                call([*compiler, '-c', str(source)], cwd=directory)
                objects.append(str(directory / (source.stem + '.o')))
            for filename in NATIVE_TESTS:
                binary = directory / Path(filename).stem
                call([*base, *objects, str(ROOT / 'test' / filename), '-o', str(binary)], cwd=directory)
                call([str(binary)], cwd=directory)
            executable = directory / 'reference_partial_coherence'
            call([*base, *objects, str(ROOT / 'tools' / 'reference_signal_partial_coherence.f90'),
                  '-o', str(executable)], cwd=directory)
            env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
            call([sys.executable, str(ROOT / 'tools' / 'compare_signal_partial_coherence.py'),
                  str(executable)], cwd=directory, env=env)
            print(f'PASS {label} build', flush=True)


if __name__ == '__main__':
    main()
