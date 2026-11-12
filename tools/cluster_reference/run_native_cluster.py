#!/usr/bin/env python3
"""Build and run the native hierarchical clustering tests without FPM."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def run(mode):
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-cluster-native-') as tmp:
        d = Path(tmp)
        command = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
                   '-J', str(d), '-I', str(d)]
        objs = []
        for name in ('kinds', 'spatial_distance', 'cluster_hierarchy'):
            obj = d / (name+'.o')
            subprocess.run([*command, '-c', str(ROOT/'src'/(name+'.f90')),
                            '-o', str(obj)], cwd=d, check=True)
            objs.append(str(obj))
        for test_name in ('test_cluster_hierarchy', 'test_cluster_hierarchy_ordering', 'test_cluster_hierarchy_flatdata', 'test_cluster_hierarchy_layout'):
            exe = d/test_name
            subprocess.run([*command, str(ROOT/'test'/(test_name+'.f90')),
                            *objs, '-o', str(exe)], cwd=d, check=True)
            result = subprocess.run([str(exe)], cwd=d, text=True, check=True,
                                    capture_output=True)
            print(f'{mode}: {result.stdout.strip()}', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    args = parser.parse_args()
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        run(mode)

if __name__=='__main__':
    main()
