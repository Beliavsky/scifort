#!/usr/bin/env python3
"""Compile and run standalone native SciFort dense matrix-equation tests."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["checked", "optimized", "both"], default="both")
    args = parser.parse_args()
    modes = ["checked", "optimized"] if args.mode == "both" else [args.mode]
    common = ["src/kinds.f90", "src/linalg.f90", "src/linalg_matrix_functions.f90", "src/linalg_matrix_equations.f90"]
    cases = [
        ("matrix_equations", common + ["test/test_linalg_matrix_equations.f90"]),
        ("schur_equations", common + ["test/test_linalg_matrix_equations_schur.f90"]),
        ("riccati", common + ["src/linalg_riccati.f90", "test/test_linalg_riccati.f90"]),
    ]
    with tempfile.TemporaryDirectory(prefix="scifort-linalg-native-") as d:
        for mode in modes:
            flags = ["-O0", "-g", "-fcheck=all", "-fbacktrace"] if mode == "checked" else ["-O3"]
            for label, sources in cases:
                exe = Path(d) / ("linalg_" + label + "_" + mode + (".exe" if os.name == "nt" else ""))
                subprocess.run(["gfortran", *flags, "-ffree-line-length-none", 
                                *[str(ROOT / s) for s in sources], "-o", str(exe)],
                               cwd=d, check=True)
                result = subprocess.run([str(exe)], cwd=d, text=True, capture_output=True, check=True)
                print(f"{mode}: {label}: {result.stdout.strip()}")


if __name__ == "__main__":
    main()
