#!/usr/bin/env python3
"""Test canonical coherence percentile/basic/bias-corrected bootstrap intervals."""
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
    'test_signal_canonical_interval_methods.f90',
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
        interval_driver = build / 'canonical_bootstrap_intervals_driver'
        subprocess.run([*flags, *objects, str(ROOT / 'tools' / 'signal_reference' /
                                               'canonical_bootstrap_intervals_driver.f90'),
                        '-o', str(interval_driver)], cwd=build, check=True)
        subprocess.run([sys.executable, str(ROOT / 'tools' / 'signal_reference' /
                                            'verify_canonical_bootstrap_intervals.py'),
                        str(interval_driver),'--max-cases',
                        '12' if optimized else '40'], cwd=build, env=env, check=True)
        if not optimized:
            coverage_driver = build / 'canonical_interval_coverage_driver'
            subprocess.run([*flags, *objects, str(ROOT / 'tools' / 'signal_reference' /
                                                'canonical_interval_coverage_driver.f90'),
                            '-o', str(coverage_driver)], cwd=build, check=True)
            subprocess.run([sys.executable, str(ROOT / 'tools' / 'signal_reference' /
                                               'study_canonical_interval_var.py'),
                            str(coverage_driver),'--trials','32'],
                           cwd=build, env=env, check=True)
        print('PASS ' + ('optimized' if optimized else 'checked') +
              ' canonical bootstrap interval methods and legacy regressions', flush=True)


if __name__ == '__main__':
    run(False)
    run(True)
