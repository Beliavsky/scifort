#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compile native tests and check SciFort matrix trig functions against SciPy."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

from compare_trigonometric_scipy import compare

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_matrix_functions.f90",
           "src/linalg_trigonometric.f90"]


def compile_executable(compiler, flags, directory, driver, name):
    executable = directory / name
    subprocess.run([compiler, "-std=f2018", *flags, "-J", str(directory),
                    "-I", str(directory), *[str(ROOT / src) for src in SOURCES],
                    str(ROOT / driver), "-o", str(executable)],
                   cwd=directory, check=True)
    return str(executable)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    parser.add_argument("--cases", type=int, default=600)
    parser.add_argument("--seed", type=int, default=20261107)
    parser.add_argument("--max-order", type=int, default=6)
    args = parser.parse_args()
    if args.cases < 1 or args.max_order < 1 or args.max_order > 40:
        parser.error("cases must be positive and max-order must be in [1, 40]")
    modes = ("checked", "optimized") if args.mode == "both" else (args.mode,)
    with tempfile.TemporaryDirectory(prefix="scifort-matrix-trig-") as folder:
        for mode in modes:
            directory = Path(folder) / mode
            directory.mkdir()
            flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked"
                     else ["-O3"])
            native = compile_executable(args.compiler, flags, directory,
                                        "test/test_linalg_trigonometric.f90", "native")
            subprocess.run([native], check=True)
            reference = compile_executable(args.compiler, flags, directory,
                                           "tools/linalg_reference/driver_trigonometric.f90", "reference")
            print(mode + ":", flush=True)
            compare(reference, args.cases, args.seed, args.max_order)


if __name__ == "__main__":
    main()
