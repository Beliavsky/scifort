#!/usr/bin/env python3
"""Compile and validate metric-aware canonical subspace tracking."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_signal_canonical_subspace_tracking import SOURCES

ROOT = Path(__file__).resolve().parents[1]
TESTS = ('test_signal_canonical_covmetric_tracking.f90',
         'test_signal_canonical_subspace_tracking.f90',
         'test_signal_canonical_tracking.f90',
         'test_signal_canonical_vectors.f90',
         'test_signal_canonical_coherence.f90')


def run(optimized: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-covmetric-') as temp:
        build = Path(temp)
        flags = ['gfortran', '-std=f2018', '-fimplicit-none',
                 '-ffree-line-length-none', '-fbacktrace']
        flags += (['-O2'] if optimized else ['-O0', '-g', '-fcheck=all'])
        objects = []
        for source in SOURCES:
            path = ROOT / 'src' / source
            obj = build / (path.stem + '.o')
            subprocess.run([*flags, '-c', str(path), '-o', str(obj)], cwd=build, check=True)
            objects.append(str(obj))
        for name in TESTS:
            exe = build / Path(name).stem
            subprocess.run([*flags, *objects, str(ROOT / 'test' / name), '-o', str(exe)],
                           cwd=build, check=True)
            subprocess.run([str(exe)], cwd=build, check=True)
        driver = build / 'canonical_covmetric_tracking_driver'
        subprocess.run([*flags, *objects, str(ROOT / 'tools' / 'signal_reference' /
                                              'canonical_covmetric_tracking_driver.f90'),
                        '-o', str(driver)], cwd=build, check=True)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, str(ROOT / 'tools' / 'signal_reference' /
                                            'verify_canonical_covmetric_tracking.py'), str(driver)],
                       cwd=build, env=env, check=True)
        # The old independently verified reference is included in the checked
        # build, while its native regressions run in both configurations.
        if not optimized:
            legacy_driver = build / 'canonical_subspace_tracking_driver'
            subprocess.run([*flags, *objects, str(ROOT / 'tools' / 'signal_reference' /
                                                  'canonical_subspace_tracking_driver.f90'),
                            '-o', str(legacy_driver)], cwd=build, check=True)
            subprocess.run([sys.executable, str(ROOT / 'tools' / 'signal_reference' /
                                                'verify_canonical_subspace_tracking.py'),
                            str(legacy_driver)], cwd=build, env=env, check=True)
        print('PASS ' + ('optimized' if optimized else 'checked') +
              ' canonical covariance-metric and native regression suites', flush=True)


if __name__ == '__main__':
    run(False)
    run(True)
