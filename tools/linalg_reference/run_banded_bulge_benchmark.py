#!/usr/bin/env python3
"""Build and time real/complex bulge reduction and real full-spectrum methods.

Single-process CPU-time results depend on hardware, compiler, load, and
matrix family. Use the independent SciPy tests for correctness assessments.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["src/kinds.f90", "src/linalg.f90", "src/linalg_banded.f90",
           "src/linalg_hessenberg.f90", "src/linalg_banded_tridiagonal.f90",
           "src/linalg_banded_eigen.f90"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="scifort_banded_benchmark_") as directory:
        work = Path(directory)
        exe = work / "bench_banded_bulge"
        subprocess.run([args.compiler, "-std=f2018", "-O3",
                        *[str(ROOT / name) for name in SOURCES],
                        str(ROOT / "tools/linalg_reference/benchmark_banded_bulge_real.f90"),
                        "-o", str(exe)], check=True, cwd=work)
        result = subprocess.run([str(exe)], text=True, capture_output=True,
                                check=True, timeout=180, cwd=work,
                                env={**os.environ,"OPENBLAS_NUM_THREADS":"1", "OMP_NUM_THREADS":"1"})
        print(result.stdout, end="")


if __name__ == "__main__":
    main()
