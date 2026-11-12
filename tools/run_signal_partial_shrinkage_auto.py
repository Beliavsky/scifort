#!/usr/bin/env python3
"""Build and test automatic Gaussian-Wishart diagonal spectral shrinkage.

Runs GNU Fortran runtime-checked and targeted -O3 builds, NumPy independent
precision/coherency comparisons, and preceding native regressions. Optional
--legacy-references also runs the previous independent numerical suites.
No installed SciFort, FPM, system BLAS or external Fortran package is needed.
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
    'test_signal_eigen_auto_shrinkage.f90',
    'test_signal_partial_shrinkage_auto.f90',
    'test_signal_partial_coherency.f90',
    'test_signal_partial_coherence.f90',
    'test_signal_tf_partial_coherence.f90',
    'test_signal_multitaper_multichannel.f90',
    'test_signal_multitaper_tf_cross.f90',
)
REFERENCES = (
    ('reference_signal_partial_shrinkage_auto.f90', 'compare_signal_partial_shrinkage_auto.py'),
    ('reference_signal_partial_coherency.f90', 'compare_signal_partial_coherency.py'),
    ('reference_signal_partial_coherence.f90', 'compare_signal_partial_coherence.py'),
    ('reference_signal_tf_partial_coherence.f90', 'compare_signal_tf_partial_coherence.py'),
)


def run(command: list[str], cwd: Path, env: dict[str,str] | None = None) -> None:
    subprocess.run(command, cwd=cwd, env=env, check=True)


def validate(optimized: bool, legacy: bool) -> None:
    with tempfile.TemporaryDirectory(prefix='scifort-auto-shrinkage-') as temp:
        build = Path(temp)
        base = ['gfortran', '-std=f2018', '-ffree-line-length-none']
        checked = base + ['-O0', '-g', '-fcheck=all', '-fbacktrace']
        target = base + ['-O3', '-march=native'] if optimized else checked
        objs = []
        for filename in SOURCES:
            source = ROOT / 'src' / filename
            flags = target if filename == 'signal_partial_coherence.f90' else checked
            run([*flags, '-c', str(source)], build)
            objs.append(str(build / (source.stem + '.o')))
        for filename in NATIVE:
            source = ROOT / 'test' / filename
            exe = build / source.stem
            run([*target, *objs, str(source), '-o', str(exe)], build)
            run([str(exe)], build)
        comparisons = REFERENCES if legacy and not optimized else REFERENCES[:1]
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   OMP_NUM_THREADS='1')
        for reference_file, python_file in comparisons:
            source = ROOT / 'tools' / reference_file
            exe = build / source.stem
            run([*target, *objs, str(source), '-o', str(exe)], build)
            run([sys.executable, str(ROOT / 'tools' / python_file), str(exe)], build, env)
        print(f'PASS {"targeted O3" if optimized else "checked"} build and regressions', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--legacy-references', action='store_true',
                        help='repeat previous 288, 252, and 120-case numerical suites in checked mode')
    args = parser.parse_args()
    validate(False, args.legacy_references)
    validate(True, False)


if __name__ == '__main__':
    main()
