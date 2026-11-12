#!/usr/bin/env python3
"""Benchmark cached/indexed Clough-Tocher against an optional previous ZIP.

The same deterministic Fortran driver is compiled for both versions. No
SciPy, internet access, FPM, BLAS, or LAPACK is required. Construction and
query timings are measured separately; results are machine dependent.

Usage:
    python tools/interpolate_reference/benchmark_clough_tocher.py \
        --baseline-zip ../scifort-interpolate-clough-tocher-cubic-2026-10-31.zip
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("kinds", "interpolate_kdtree", "interpolate_scattered", "interpolate_clough_tocher")
DRIVER = ROOT / "tools/interpolate_reference/benchmark_clough_tocher_driver.f90"


def build_and_run(source_root: Path, temporary: Path) -> dict:
    temporary.mkdir(parents=True, exist_ok=True)
    cmd = ["gfortran", "-std=f2018", "-O3", "-ffree-line-length-none", "-I" + str(temporary), "-J" + str(temporary)]
    objects = []
    for name in SOURCES:
        obj = temporary / f"{name}.o"
        subprocess.run([*cmd, "-c", str(source_root / "src" / (name + ".f90")), "-o", str(obj)],
                       cwd=temporary, check=True, capture_output=True, text=True)
        objects.append(obj)
    exe = temporary / "benchmark"
    subprocess.run([*cmd, str(DRIVER), *map(str, objects), "-o", str(exe)], cwd=temporary,
                   check=True, capture_output=True, text=True)
    result = subprocess.run([str(exe)], check=True, capture_output=True, text=True, timeout=240)
    return dict(line.strip().split("=", 1) for line in result.stdout.splitlines())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-zip", type=Path, help="Prior complete SciFort ZIP to compare")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as td:
        temp = Path(td)
        current = build_and_run(ROOT, temp / "current_build")
        print("Indexed/cached implementation:", current, flush=True)
        if args.baseline_zip:
            with zipfile.ZipFile(args.baseline_zip) as archive:
                if not all(("scifort-main/src/" + name + ".f90") in archive.namelist() for name in SOURCES):
                    raise ValueError("The baseline ZIP is missing required modules")
                for name in SOURCES:
                    archive.extract("scifort-main/src/" + name + ".f90", temp / "baseline")
            baseline = build_and_run(temp / "baseline/scifort-main", temp / "baseline_build")
            print("Previous implementation:", baseline, flush=True)
            if int(baseline["triangles"]) != int(current["triangles"]):
                raise AssertionError("Different triangulations: benchmark is not comparable")
            a, b = float(current["checksum"]), float(baseline["checksum"])
            if abs(a-b) > 1e-9 * max(1.0, abs(b)):
                raise AssertionError(f"Benchmark checksum mismatch: {a} vs {b}")
            ratio = float(baseline["evaluate_sec"]) / float(current["evaluate_sec"])
            print(f"Evaluation speedup: {ratio:.2f}x (CPU times exclude construction)", flush=True)


if __name__ == "__main__":
    main()
