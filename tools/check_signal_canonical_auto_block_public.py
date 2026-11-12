#!/usr/bin/env python3
"""Compile full scifort_signal dependency closure and all bootstrap public tests."""
from __future__ import annotations

import argparse
from pathlib import Path
import tempfile

from check_signal_canonical_moving_bootstrap_public import TESTS as LEGACY_TESTS
from check_signal_mixed_public import closure
import subprocess

ROOT = Path(__file__).resolve().parents[1]
TESTS = ('test_signal_canonical_auto_block_public.f90', *LEGACY_TESTS)


def check(build: Path) -> None:
    build.mkdir(parents=True,exist_ok=True)
    order = closure()
    flags = ['gfortran','-std=f2018','-fimplicit-none',
             '-ffree-line-length-none','-O0','-J.','-I.']
    objects = []
    print(f'Canonical auto-block public dependency closure: {len(order)} modules',flush=True)
    for index,source in enumerate(order,1):
        obj = build / (source.stem + '.o')
        subprocess.run([*flags,'-c',str(source),'-o',str(obj)],cwd=build,check=True)
        objects.append(str(obj))
        if index % 20 == 0 or index == len(order):
            print(f'Compiled {index}/{len(order)} public modules',flush=True)
    for name in TESTS:
        exe = build / Path(name).stem
        subprocess.run([*flags,str(ROOT/'test'/name),*objects,'-o',str(exe)],
                       cwd=build,check=True)
        subprocess.run([str(exe)],cwd=build,check=True)
    print(f'PASS {len(TESTS)} canonical public-interface regression suites',flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir',type=Path)
    args = parser.parse_args()
    if args.build_dir is not None:
        check(args.build_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='scifort-auto-block-public-') as folder:
            check(Path(folder))
