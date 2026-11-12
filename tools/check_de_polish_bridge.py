#!/usr/bin/env python3
"""Compile DE polish bridge using a TEST DOUBLE, not real L-BFGS-B.

The full numerical L-BFGS-B backend is part of the FPM dependency graph;
this smoke test checks only argument compatibility and optional callback wiring.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fc', default='gfortran')
    parser.add_argument('--mode', choices=('checked','optimized'), default='checked')
    args = parser.parse_args()
    opts = (['-O0','-fcheck=all','-fbacktrace'] if args.mode=='checked'
            else ['-O3','-march=native'])
    with tempfile.TemporaryDirectory(prefix='scifort_de_polish_bridge_') as tmp:
        exe = str(Path(tmp)/'bridge')
        cmd = [args.fc,'-std=f2018','-Wall','-Wextra','-Werror',
               '-J',tmp,'-I',tmp,*opts,'-o',exe,
               str(ROOT/'src/kinds.f90'),
               str(ROOT/'src/optimize/differential_evolution.f90'),
               str(ROOT/'tools/fixtures/mock_de_polish_minimize.f90'),
               str(ROOT/'src/optimize/differential_evolution_polish.f90'),
               str(ROOT/'tools/fixtures/test_de_polish_bridge.f90')]
        subprocess.run(cmd,check=True)
        subprocess.run([exe],check=True)
    print('DE L-BFGS-B POLISH BRIDGE INTERFACE CHECK PASSED (mock backend only)')

if __name__=='__main__':
    main()
