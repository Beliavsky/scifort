#!/usr/bin/env python3
"""Validate full scifort_signal public dependency closure and interface tests."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
from check_signal_mixed_public import closure

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_TESTS = ('test_signal_multiple_coherence_public.f90',
                'test_signal_partial_shrinkage_cv_public.f90',
                'test_signal_partial_shrinkage_auto_public.f90',
                'test_signal_partial_coherency_public.f90',
                'test_signal_partial_coherence_public.f90',
                'test_signal_tf_partial_coherence_public.f90')


def check(folder: Path) -> None:
    folder.mkdir(parents=True,exist_ok=True)
    order=closure()
    flags=['gfortran','-std=f2018','-fimplicit-none','-ffree-line-length-none',
           '-O0','-J.','-I.']
    print(f'Public-module dependency closure: {len(order)} sources',flush=True)
    for index,source in enumerate(order,1):
        obj=folder/(source.stem+'.o')
        if obj.exists():
            continue
        subprocess.run([*flags,'-c',str(source),'-o',str(obj)],cwd=folder,check=True)
        if index%15 == 0 or index == len(order):
            print(f'Compiled {index}/{len(order)}',flush=True)
    objects=[str(folder/(source.stem+'.o')) for source in order]
    for test in PUBLIC_TESTS:
        exe=folder/Path(test).stem
        subprocess.run([*flags,str(ROOT/'test'/test),*objects,'-o',str(exe)],
                       cwd=folder,check=True)
        subprocess.run([str(exe)],cwd=folder,check=True)
    print(f'PASS {len(PUBLIC_TESTS)} public multiple/partial coherence tests',flush=True)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir',type=Path)
    args=parser.parse_args()
    if args.build_dir is None:
        with tempfile.TemporaryDirectory(prefix='scifort-multiple-public-') as tmp:
            check(Path(tmp))
    else:
        check(args.build_dir.resolve())


if __name__ == '__main__':
    main()
