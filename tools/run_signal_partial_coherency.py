#!/usr/bin/env python3
"""Run NumPy/SciPy checks for phase-aware, diagonal-shrunk partial coherency.

Builds checked and targeted -O3 versions of the modified Fortran module,
plus previous coherence and multitaper native regressions. The larger
legacy numerical comparisons can be included with --legacy-references. No FPM,
external Fortran BLAS, or installed project library is required.
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
    'kinds.f90', 'signal_dpss.f90', 'signal_windows_extra.f90',
    'signal_spectral.f90', 'signal_multitaper.f90',
    'signal_partial_coherence.f90',
)
NATIVE = (
    'test_signal_partial_coherency.f90',
    'test_signal_partial_coherence.f90',
    'test_signal_tf_partial_coherence.f90',
    'test_signal_multitaper_multichannel.f90',
    'test_signal_multitaper_tf_cross.f90',
)
REFERENCES = (
    ('reference_signal_partial_coherency.f90', 'compare_signal_partial_coherency.py'),
    ('reference_signal_partial_coherence.f90', 'compare_signal_partial_coherence.py'),
    ('reference_signal_tf_partial_coherence.f90', 'compare_signal_tf_partial_coherence.py'),
)


def run(args: list[str], cwd: Path, env: dict[str,str]|None = None) -> None:
    subprocess.run(args, cwd=cwd, env=env, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--legacy-references', action='store_true',
                        help='also repeat the two older NumPy/SciPy reference suites in the checked build')
    args = parser.parse_args()
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1')
    for name, optimized in [('checked', False), ('targeted-optimized', True)]:
        with tempfile.TemporaryDirectory(prefix='scifort-complex-partial-') as temp:
            build = Path(temp)
            base = ['gfortran', '-std=f2018', '-ffree-line-length-none']
            checked = [*base, '-O0', '-g', '-fcheck=all', '-fbacktrace']
            target = [*base, '-O3', '-march=native'] if optimized else checked
            objects = []
            for filename in SOURCES:
                flags = target if filename == 'signal_partial_coherence.f90' else checked
                source = ROOT / 'src' / filename
                run([*flags, '-c', str(source)], build)
                objects.append(str(build / (source.stem + '.o')))
            for filename in NATIVE:
                exe = build / Path(filename).stem
                run([*target, *objects, str(ROOT / 'test' / filename), '-o', str(exe)], build)
                run([str(exe)], build)
            comparisons = REFERENCES if args.legacy_references and not optimized else REFERENCES[:1]
            for fortran_file, compare_file in comparisons:
                exe = build / Path(fortran_file).stem
                run([*target, *objects, str(ROOT / 'tools' / fortran_file), '-o', str(exe)], build)
                run([sys.executable, str(ROOT / 'tools' / compare_file), str(exe)],
                    build, env=env)
            print(f'PASS {name} complete', flush=True)


if __name__ == '__main__':
    main()
