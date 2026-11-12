#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compile and run SciFort expm_cond native and SciPy reference regressions."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

from compare_expm_condition_scipy import compare
from compare_expm_condition_stress import compare as compare_stress

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_expm_frechet.f90",
           "src/linalg_expm_condition.f90"]


def compile_executable(compiler, flags, directory, driver, name):
    executable = directory / name
    subprocess.run([compiler, "-std=f2018", *flags, "-J", str(directory),
                    "-I", str(directory), *[str(ROOT / src) for src in SOURCES],
                    str(ROOT / driver), "-o", str(executable)],
                   cwd=directory, check=True)
    return str(executable)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--compiler", default="gfortran")
    p.add_argument("--mode", choices=("checked", "optimized", "both"), default="both")
    p.add_argument("--cases", type=int, default=300)
    p.add_argument("--seed", type=int, default=20261107)
    p.add_argument("--max-order", type=int, default=6)
    p.add_argument("--skip-stress", action="store_true", help="skip analytic stress checks")
    args = p.parse_args()
    if args.cases < 1 or args.max_order < 1 or args.max_order > 48:
        p.error("cases must be positive and max-order in [1,48]")
    modes = ("checked", "optimized") if args.mode == "both" else (args.mode,)
    with tempfile.TemporaryDirectory(prefix="scifort-expm-cond-") as folder:
        for mode in modes:
            directory = Path(folder) / mode
            directory.mkdir()
            flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked"
                     else ["-O3"])
            native = compile_executable(args.compiler, flags, directory,
                                        "test/test_linalg_expm_condition.f90", "native")
            subprocess.run([native], check=True)
            reference = compile_executable(args.compiler, flags, directory,
                                           "tools/linalg_reference/driver_expm_condition.f90",
                                           "reference")
            print(mode + ":", flush=True)
            compare(reference, args.cases, args.seed, args.max_order)
            if not args.skip_stress:
                compare_stress(reference, args.seed)


if __name__ == "__main__":
    main()
