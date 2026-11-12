#!/usr/bin/env python3
"""Check automatic-shrinkage API through scifort_signal's full dependency closure.

Optional --build-dir permits reusing unchanged object files from a previous
validated build. Remove changed source objects before reusing any such build.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
from check_signal_mixed_public import closure

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_TESTS = (
    'test_signal_partial_shrinkage_cv_public.f90',
    'test_signal_partial_shrinkage_auto_public.f90',
    'test_signal_partial_coherency_public.f90',
    'test_signal_partial_coherence_public.f90',
    'test_signal_tf_partial_coherence_public.f90',
)


def check(build: Path, limit: int) -> None:
    ordered = closure()
    build.mkdir(parents=True, exist_ok=True)
    print(f'Public-module dependency closure: {len(ordered)} source files', flush=True)
    flags = ['gfortran', '-std=f2018', '-ffree-line-length-none',
             '-fimplicit-none', '-O0', '-J.', '-I.']
    for i, source in enumerate(ordered, 1):
        obj = build / (source.stem+'.o')
        if obj.exists():
            continue
        if i > limit:
            print(f'Partial build: stopped at {i-1}/{len(ordered)}', flush=True)
            return
        subprocess.run([*flags, '-c', str(source), '-o', str(obj)],
                       cwd=build, check=True)
        if i%15 == 0 or i == len(ordered):
            print(f'Compiled {i}/{len(ordered)}', flush=True)
    objs = [str(build / (source.stem+'.o')) for source in ordered]
    for name in PUBLIC_TESTS:
        exe = build / Path(name).stem
        subprocess.run([*flags, str(ROOT/'test'/name), *objs, '-o', str(exe)],
                       cwd=build, check=True)
        subprocess.run([str(exe)], cwd=build, check=True)
    print('PASS five public partial-coherence/shrinkage tests',flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path)
    parser.add_argument('--count', type=int, default=999)
    options = parser.parse_args()
    if options.build_dir is None:
        with tempfile.TemporaryDirectory(prefix='scifort-auto-shrinkage-public-') as temp:
            check(Path(temp), options.count)
    else:
        check(options.build_dir.resolve(), options.count)


if __name__ == '__main__':
    main()
