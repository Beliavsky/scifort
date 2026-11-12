#!/usr/bin/env python3
"""Compile the public scifort_signal module and test new/legacy partial coherence.

A --build-dir may be supplied to keep compiler objects between interrupted
runs. --count compiles only the first N ordered dependencies, if requested.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

from check_signal_mixed_public import closure

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_TESTS = (
    'test_signal_partial_coherency_public.f90',
    'test_signal_partial_coherence_public.f90',
    'test_signal_tf_partial_coherence_public.f90',
)


def check(build: Path, count: int) -> None:
    order = closure()
    print(f'Public scifort_signal dependency closure: {len(order)} modules', flush=True)
    build.mkdir(parents=True, exist_ok=True)
    flags = ['gfortran', '-std=f2018', '-ffree-line-length-none', '-O0',
             '-fimplicit-none', '-J.', '-I.']
    for i, src in enumerate(order, start=1):
        obj = build / (src.stem + '.o')
        if obj.exists():
            continue
        if i > count:
            print(f'Partial build: {i-1} of {len(order)} source modules available', flush=True)
            return
        subprocess.run([*flags, '-c', str(src), '-o', str(obj)], cwd=build, check=True)
        if i % 15 == 0 or i == len(order):
            print(f'Compiled {i}/{len(order)}', flush=True)
    objects = [str(build / (p.stem + '.o')) for p in order]
    for filename in PUBLIC_TESTS:
        exe = build / Path(filename).stem
        subprocess.run([*flags, str(ROOT / 'test' / filename), *objects, '-o', str(exe)],
                       cwd=build, check=True)
        subprocess.run([str(exe)], cwd=build, check=True)
    print('PASS all three public partial-coherence interfaces', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path)
    parser.add_argument('--count', type=int, default=999)
    args = parser.parse_args()
    if args.build_dir is None:
        with tempfile.TemporaryDirectory(prefix='scifort-public-coherency-') as tmp:
            check(Path(tmp), args.count)
    else:
        check(args.build_dir.resolve(), args.count)


if __name__ == '__main__':
    main()
