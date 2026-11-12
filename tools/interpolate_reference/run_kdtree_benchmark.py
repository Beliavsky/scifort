#!/usr/bin/env python3
"""Compile and run the independent SciFort KD-tree nearest-neighbor benchmark.

Usage: python tools/interpolate_reference/run_kdtree_benchmark.py [sites] [queries] [dimension]
"""
from __future__ import annotations
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]


def run(args,cwd):
    proc=subprocess.run(args,cwd=cwd,text=True,capture_output=True)
    if proc.returncode:
        raise RuntimeError(proc.stderr[-3000:]+proc.stdout[-3000:])
    return proc.stdout


def main():
    with tempfile.TemporaryDirectory() as td:
        directory=Path(td)
        files=['src/kinds.f90','src/interpolate_kdtree.f90',
               'tools/interpolate_reference/benchmark_kdtree.f90']
        command=['gfortran','-std=f2018','-ffree-line-length-none','-O3',
                 '-I'+td,'-J'+td,*(str(ROOT/f) for f in files),'-o',str(directory/'bench')]
        run(command,directory)
        print(run([str(directory/'bench'),*sys.argv[1:]],directory))


if __name__=='__main__':
    main()
