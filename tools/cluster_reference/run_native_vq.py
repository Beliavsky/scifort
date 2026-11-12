#!/usr/bin/env python3
"""Compile and run the native vector-quantization regression without FPM."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def run(mode):
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-vq-native-') as temp:
        folder = Path(temp)
        compiler = ['gfortran', '-std=f2018', *flags, '-J', str(folder), '-I', str(folder)]
        objects = []
        for name in ('kinds','random','cluster_vq'):
            obj = folder / f'{name}.o'
            subprocess.run([*compiler, '-c', str(ROOT/'src'/f'{name}.f90'),
                            '-o', str(obj)], cwd=folder, check=True)
            objects.append(str(obj))
        for test_name in ('test_cluster_vq', 'test_cluster_vq_gaussian'):
            exe = folder/test_name
            subprocess.run([*compiler, str(ROOT/'test'/(test_name+'.f90')),
                            *objects, '-o', str(exe)], cwd=folder, check=True)
            result = subprocess.run([str(exe)], cwd=folder, text=True, check=True,
                                    capture_output=True)
            print(f'{mode}: {result.stdout.strip()}', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    args = parser.parse_args()
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        run(mode)


if __name__=='__main__':
    main()
