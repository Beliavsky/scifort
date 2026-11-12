#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build all public signal modules and test the new VAR(1) spectral target API."""
from __future__ import annotations

import argparse
from pathlib import Path

from check_signal_canonical_bca_public import TESTS as OLD_TESTS
from check_signal_mixed_public import closure

ROOT = Path(__file__).resolve().parents[1]


def check(build: Path) -> None:
    # Keep the previous public regression suite and add the exact target.
    import subprocess
    order = closure()
    build.mkdir(parents=True,exist_ok=True)
    flags = ['gfortran','-std=f2018','-fimplicit-none','-ffree-line-length-none',
             '-O0','-J.','-I.']
    objects = []
    print(f'Finite-window VAR(1) public dependency closure: {len(order)} modules',flush=True)
    for i,source in enumerate(order,1):
        obj = build/(source.stem+'.o')
        subprocess.run([*flags,'-c',str(source),'-o',str(obj)],cwd=build,check=True)
        objects.append(str(obj))
        if i%20 == 0 or i == len(order):
            print(f'Compiled {i}/{len(order)} public modules',flush=True)
    for name in ('test_signal_var1_taper_spectrum_public.f90',*OLD_TESTS):
        exe = build/Path(name).stem
        subprocess.run([*flags,str(ROOT/'test'/name),*objects,'-o',str(exe)],
                       cwd=build,check=True)
        subprocess.run([str(exe)],cwd=build,check=True)
    print(f'PASS {len(OLD_TESTS)+1} public signal interface tests',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir',type=Path)
    args = parser.parse_args()
    if args.build_dir:
        check(args.build_dir.resolve())
    else:
        import tempfile
        with tempfile.TemporaryDirectory(prefix='scifort-var1-public-') as tmp:
            check(Path(tmp))
