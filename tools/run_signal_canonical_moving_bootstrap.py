#!/usr/bin/env python3
"""Run checked and optimized moving-block bootstrap and legacy numerical tests."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_signal_canonical_jackknife import SOURCES as JACKKNIFE_SOURCES

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (*JACKKNIFE_SOURCES, 'signal_canonical_block_selection.f90',
           'signal_canonical_bootstrap.f90')
TESTS = (
    'test_signal_canonical_moving_bootstrap.f90',
    'test_signal_canonical_block_jackknife.f90',
    'test_signal_canonical_jackknife.f90',
    'test_signal_canonical_coherence.f90',
    'test_signal_group_coherence.f90',
    'test_signal_band_coherence.f90',
)


def run(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-canonical-bootstrap-') as directory:
        build = Path(directory)
        flags = ['gfortran', '-std=f2018', '-fimplicit-none',
                 '-ffree-line-length-none', '-fbacktrace']
        flags += ['-O2'] if optimized else ['-O0', '-g', '-fcheck=all']
        objects = []
        for source in SOURCES:
            path = ROOT / 'src' / source
            obj = build / (path.stem + '.o')
            subprocess.run([*flags, '-c', str(path), '-o', str(obj)],
                           cwd=build, check=True)
            objects.append(str(obj))
        for name in TESTS:
            exe = build / Path(name).stem
            subprocess.run([*flags, *objects, str(ROOT / 'test' / name),
                            '-o', str(exe)], cwd=build, check=True)
            subprocess.run([str(exe)], cwd=build, check=True)
        driver = build / 'canonical_bootstrap_driver'
        subprocess.run([*flags, *objects, str(ROOT / 'tools' / 'signal_reference' /
                                              'canonical_bootstrap_driver.f90'),
                        '-o', str(driver)], cwd=build, check=True)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, str(ROOT / 'tools' / 'signal_reference' /
                                            'verify_canonical_moving_bootstrap.py'), str(driver)],
                       cwd=build, env=env, check=True)
        print('PASS ' + ('optimized' if optimized else 'checked') +
              ' canonical moving-block bootstrap and previous regressions', flush=True)


if __name__ == '__main__':
    run(False)
    run(True)
