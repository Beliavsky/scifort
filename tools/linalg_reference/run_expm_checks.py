#!/usr/bin/env python3
"""Build/run standalone SciFort expm tests with gfortran, no FPM required.

From the repository root:
    python tools/linalg_reference/run_expm_checks.py

Requires gfortran, NumPy, SciPy. Leaves the original tree unchanged.
Works on Windows, Linux, and macOS where gfortran is installed.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCES = [ROOT / 'src' / name for name in (
    'kinds.f90', 'linalg.f90', 'linalg_dense.f90', 'linalg_expm.f90')]
TOOL = ROOT / 'tools' / 'linalg_reference'


def run(args: list[str], cwd: Path) -> None:
    print('RUN:', ' '.join(args), flush=True)
    subprocess.run(args, cwd=cwd, check=True)


def main() -> None:
    compiler = shutil.which('gfortran')
    if compiler is None:
        raise SystemExit('gfortran was not found in PATH')
    with tempfile.TemporaryDirectory(prefix='scifort-expm-') as tmp:
        work = Path(tmp)
        exe = '.exe' if sys.platform == 'win32' else ''
        for config, flags in (
            ('checked', ['-O0', '-g', '-fcheck=all', '-fbacktrace']),
            ('optimized', ['-O3', '-march=native']),
        ):
            build = work / config
            build.mkdir()
            common = [compiler, '-std=f2018', *flags]
            # One source invocation builds modules in dependency order.
            run([*common, '-c', *map(str, SOURCES)], build)
            objects = [str(build / f'{s.stem}.o') for s in SOURCES]
            for name, source in (
                ('test_linalg_expm', ROOT/'test'/'test_linalg_expm.f90'),
                ('driver_expm', TOOL/'driver_expm.f90'),
                ('linalg_expm_example', ROOT/'example'/'linalg_expm_example.f90'),
                ('test_linalg_dense', ROOT/'test'/'test_linalg_dense.f90'),
            ):
                run([*common, '-o', str(build/(name+exe)), str(source), *objects], build)
            run([str(build/('test_linalg_expm'+exe))], build)
            run([str(build/('test_linalg_dense'+exe))], build)
            run([str(build/('linalg_expm_example'+exe))], build)
            run([sys.executable, str(TOOL/'compare_expm_scipy.py'),
                 str(build/('driver_expm'+exe))], build)
    print('PASS: standalone expm suite in checked and optimized builds')


if __name__ == '__main__':
    main()
