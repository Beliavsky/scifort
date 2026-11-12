#!/usr/bin/env python3
"""Build complete umbrella module and verify canonical-subspace public API."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

from check_signal_mixed_public import closure
from check_signal_canonical_tracking_public import TESTS as PREVIOUS_TESTS

ROOT = Path(__file__).resolve().parents[1]
TESTS = ('test_signal_canonical_subspace_tracking_public.f90', *PREVIOUS_TESTS)


def check(build: Path) -> None:
    build.mkdir(parents=True, exist_ok=True)
    order = closure()
    flags = ['gfortran', '-std=f2018', '-fimplicit-none',
             '-ffree-line-length-none', '-O0', '-J.', '-I.']
    print(f'Public-module dependency closure: {len(order)} source files', flush=True)
    objects = []
    for i, source in enumerate(order, 1):
        obj = build / (source.stem + '.o')
        subprocess.run([*flags, '-c', str(source), '-o', str(obj)], cwd=build, check=True)
        objects.append(str(obj))
        if i % 15 == 0 or i == len(order):
            print(f'Compiled {i}/{len(order)}', flush=True)
    for name in TESTS:
        exe = build / Path(name).stem
        subprocess.run([*flags, str(ROOT / 'test' / name), *objects, '-o', str(exe)],
                       cwd=build, check=True)
        subprocess.run([str(exe)], cwd=build, check=True)
    print(f'PASS {len(TESTS)} canonical subspace and preceding public tests', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path)
    args = parser.parse_args()
    if args.build_dir:
        check(args.build_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='scifort-subspace-public-') as temp:
            check(Path(temp))
