#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build and validate exact VAR(p) fixed-weight multitaper population targets.

Requires gfortran; numerical references additionally need NumPy and SciPy.
Use --checked-only / --optimized-only or --cases N as required.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / 'tools' / 'signal_reference'
MODULE = re.compile(r'^\s*module\s+(?!procedure\b|function\b|subroutine\b)(\w+)', re.I | re.M)
USE = re.compile(r'^\s*use\s+(scifort_\w+)', re.I | re.M)


def ordered_dependencies():
    mapping = {}
    for file in (ROOT / 'src').rglob('*.f90'):
        for name in MODULE.findall(file.read_text()):
            mapping[name.lower()] = file
    done, order = set(), []

    def include(name):
        if name in done:
            return
        done.add(name)
        file = mapping[name]
        for dep in USE.findall(file.read_text()):
            include(dep.lower())
        if file not in order:
            order.append(file)

    include('scifort_signal_varp_multitaper_target')
    return order


def run(optimized: bool, cases: int):
    with tempfile.TemporaryDirectory(prefix='scifort-varp-multitaper-') as tmp:
        build = Path(tmp)
        flags = ['gfortran', '-std=f2018', '-fimplicit-none',
                 '-ffree-line-length-none', '-fbacktrace', '-J.', '-I.']
        flags += ['-O2'] if optimized else ['-O0', '-g', '-fcheck=all']
        objects = []
        for source in ordered_dependencies():
            obj = build / (source.stem + '.o')
            subprocess.run([*flags, '-c', str(source), '-o', str(obj)],
                           cwd=build, check=True)
            objects.append(str(obj))
        for name in ('test_signal_varp_multitaper_target',
                     'varp_multitaper_target_driver'):
            source = (ROOT / 'test' if name.startswith('test_') else REF) / (name + '.f90')
            exe = build / name
            subprocess.run([*flags, str(source), *objects, '-o', str(exe)],
                           cwd=build, check=True)
            if name.startswith('test_'):
                subprocess.run([str(exe)], cwd=build, check=True)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, str(REF / 'verify_varp_multitaper_target.py'),
                        str(build / 'varp_multitaper_target_driver'),
                        '--cases', str(cases)], cwd=build, env=env, check=True)
        print('PASS '+('optimized' if optimized else 'checked')+
              ' fixed-weight VAR(p) multitaper references', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--checked-only', action='store_true')
    mode.add_argument('--optimized-only', action='store_true')
    parser.add_argument('--cases', type=int, default=40)
    args = parser.parse_args()
    if not args.optimized_only:
        run(False, args.cases)
    if not args.checked_only:
        run(True, args.cases)
