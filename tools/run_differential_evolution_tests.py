#!/usr/bin/env python3
"""Compile and validate SciFort differential evolution in both Fortran modes."""
from __future__ import annotations
import argparse
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--mode', choices=['checked', 'optimized', 'both'], default='both')
    args = p.parse_args()
    for mode in ('checked', 'optimized'):
        if args.mode not in ('both', mode):
            continue
        opts = ['-O0', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3', '-march=native']
        with tempfile.TemporaryDirectory(prefix='scifort_de_native_') as tmp:
            exe = pathlib.Path(tmp) / 'de_native'
            cmd = [args.fc, '-std=f2018', '-Wall', '-Wextra', '-Werror', '-J', tmp, '-I', tmp,
                   *opts, '-o', str(exe),
                   str(ROOT / 'src/kinds.f90'),
                   str(ROOT / 'src/optimize/differential_evolution.f90'),
                   str(ROOT / 'test/optimize_differential_evolution_reference.f90')]
            subprocess.run(cmd, check=True)
            print(f'--- {mode} native ---', flush=True)
            subprocess.run([str(exe)], check=True)
            ext = pathlib.Path(tmp) / 'de_constraints_native'
            cmd[-1] = str(ROOT / 'test/optimize_de_constraints_reference.f90')
            cmd[cmd.index(str(exe))] = str(ext)
            subprocess.run(cmd, check=True)
            print(f'--- {mode} integer/constraint native ---', flush=True)
            subprocess.run([str(ext)], check=True)
        print(f'--- {mode} SciPy reference ---', flush=True)
        subprocess.run([sys.executable, str(ROOT / 'test/compare_optimize_differential_evolution_scipy.py'),
                        '--fc', args.fc, '--mode', mode], check=True)
        print(f'--- {mode} mixed-integer/constraint SciPy reference ---', flush=True)
        subprocess.run([sys.executable, str(ROOT / 'test/compare_de_constraints_scipy.py'),
                        '--fc', args.fc, '--mode', mode], check=True)
        print(f'--- {mode} L-BFGS-B polish interface smoke test (test double) ---', flush=True)
        subprocess.run([sys.executable, str(ROOT / 'tools/check_de_polish_bridge.py'),
                        '--fc', args.fc, '--mode', mode], check=True)
    print('DIFFERENTIAL EVOLUTION VALIDATION PASSED')


if __name__ == '__main__':
    main()
