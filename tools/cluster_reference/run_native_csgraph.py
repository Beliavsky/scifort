#!/usr/bin/env python3
"""Compile and run CSR graph and named hierarchical-linkage native tests."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def run(mode: str) -> None:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-csgraph-') as temp:
        build = Path(temp)
        cc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
              '-I', str(build), '-J', str(build)]
        for source in ('kinds', 'spatial_distance', 'cluster_hierarchy', 'sparse_csgraph'):
            subprocess.run([*cc, '-c', str(ROOT/'src'/f'{source}.f90'),
                            '-o', str(build/f'{source}.o')], cwd=build, check=True)
        for test, modules in (
            ('test_sparse_csgraph', ('kinds', 'sparse_csgraph')),
            ('test_sparse_csgraph_extensions', ('kinds', 'sparse_csgraph')),
            ('test_sparse_csgraph_advanced', ('kinds', 'sparse_csgraph')),
            ('test_sparse_csgraph_weighted', ('kinds', 'sparse_csgraph')),
            ('test_sparse_csgraph_yen', ('kinds', 'sparse_csgraph')),
            ('test_cluster_hierarchy_shortcuts', ('kinds', 'spatial_distance', 'cluster_hierarchy')),
        ):
            exe = build/test
            subprocess.run([*cc, str(ROOT/'test'/f'{test}.f90'),
                            *(str(build/f'{m}.o') for m in modules), '-o', str(exe)],
                           cwd=build, check=True)
            result = subprocess.run([str(exe)], cwd=build, text=True,
                                    capture_output=True, check=True)
            print(f'{mode}: {result.stdout.strip()}', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    args = parser.parse_args()
    for mode in (('checked', 'optimized') if args.mode == 'both' else (args.mode,)):
        run(mode)


if __name__ == '__main__':
    main()
