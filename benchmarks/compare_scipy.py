#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Validate and time matching native SciFort and SciPy workloads (not CI gates)."""
from __future__ import annotations

import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

# Set before NumPy/SciPy load their native libraries, including child processes.
THREAD_ENV = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS")
for variable in THREAD_ENV:
    os.environ[variable] = "1"

import numpy as np
import scipy
from scipy import fft, interpolate, linalg, special

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ["src/kinds.f90", "src/constants.f90", "src/math.f90",
           "src/stats/normal.f90", "src/linalg.f90", "src/linalg_dense.f90",
           "src/signal_windows_extra.f90", "src/signal_dpss.f90",
           "src/signal_spectral.f90", "src/interpolate_barycentric.f90"]
DRIVER = "benchmarks/scipy_speed_driver.f90"
# label, native mode, ordinary size, smoke-check size
CASES = [
    ("normal_cdf_array", "normal_cdf_array", 65536, 256),
    ("normal_cdf_scalar", "normal_cdf_scalar", 256, 16),
    ("normal_ppf_array", "normal_ppf_array", 16384, 128),
    ("fft_power2_reuse", "fft_reuse", 4096, 64),
    ("fft_prime_reuse", "fft_reuse", 4093, 67),
    ("fft_prime_setup", "fft_setup", 4093, 67),
    ("lu_factor_solve", "lu_factor_solve", 128, 16),
    ("lu_reuse", "lu_reuse", 128, 16),
    ("barycentric", "barycentric", 4096, 64),
]


def build(compiler: str, checked: bool) -> tuple[Path, dict]:
    resolved = shutil.which(compiler)
    if resolved is None:
        raise RuntimeError(f"Compiler not found: {compiler}")
    version = subprocess.check_output([resolved, "--version"], text=True).strip()
    flags = ["-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
             "-fbacktrace"] if checked else ["-O3"]
    flags += ["-std=f2018", "-fimplicit-none", "-Werror=implicit-interface"]
    digest = hashlib.sha256((resolved + version + repr(flags)).encode())
    for source in [*SOURCES, DRIVER]:
        digest.update((ROOT / source).read_bytes())
    directory = ROOT / "build" / "scipy-speed" / digest.hexdigest()[:16]
    directory.mkdir(parents=True, exist_ok=True)
    exe = directory / ("compare.exe" if os.name == "nt" else "compare")
    compile_command = [resolved, *flags, "-J", str(directory), "-I", str(directory),
                       "-c", *[str(ROOT / s) for s in SOURCES]]
    link_command = [resolved, *flags, "-I", str(directory), str(ROOT / DRIVER),
                    *[str(directory / (Path(s).stem + ".o")) for s in SOURCES],
                    "-o", str(exe)]
    if not exe.exists():
        print("Compiling benchmark modules (excluded from timings)...", flush=True)
        subprocess.run(compile_command, cwd=directory, check=True)
        subprocess.run(link_command, cwd=directory, check=True)
    return exe, {"compiler": version, "flags": flags, "source_sha256": digest.hexdigest(),
                 "commands": [compile_command, link_command], "profile": "checked" if checked else "optimized"}


def prepare(mode: str, n: int):
    """Allocate two deterministic input sets before timing; return one operation."""
    t = (np.arange(n, dtype=np.float64)[:, None] + np.array([0., .25])) / n
    values = np.empty(n)
    if mode.startswith("normal"):
        inputs = (1e-6 + (1 - 2e-6) * t) if mode == "normal_ppf_array" else -8 + 16*t
        inputs = np.asfortranarray(inputs)
        function = special.ndtri if mode == "normal_ppf_array" else special.ndtr

        def calculate(k):
            if mode == "normal_cdf_scalar":
                for i in range(n):
                    values[i] = special.ndtr(inputs[i, k])
            else:
                function(inputs[:, k], out=values)
            return values

        return calculate, inputs, None
    if mode.startswith("fft"):
        base = np.arange(n, dtype=np.float64) / n
        inputs = [np.sin(17*base) + .1*k + 1j*np.cos(11*base) for k in range(2)]

        def calculate(k):
            return fft.fft(inputs[k], norm="backward", overwrite_x=False, workers=1)

        return calculate, None, None
    if mode.startswith("lu"):
        index = np.arange(1, n+1, dtype=np.float64)
        matrix = np.asfortranarray(1/(index[:, None] + index[None, :]) + n*np.eye(n))
        exact = [np.sin(index) + .1*k for k in range(2)]
        rhs = [matrix @ x for x in exact]
        factors = linalg.lu_factor(matrix, overwrite_a=False, check_finite=True)

        def calculate(k):
            packed = (linalg.lu_factor(matrix, overwrite_a=False, check_finite=True)
                      if mode == "lu_factor_solve" else factors)
            return linalg.lu_solve(packed, rhs[k], overwrite_b=False, check_finite=True)

        return calculate, (matrix, rhs), exact
    if mode == "barycentric":
        nodes = np.cos(np.pi*np.arange(33)/32)
        samples = 1/(1+25*nodes**2)
        polynomial = interpolate.BarycentricInterpolator(nodes, samples, rng=20261122)
        queries = np.asfortranarray(-.99 + 1.98*t)

        def calculate(k):
            return polynomial(queries[:, k])

        return calculate, None, None
    raise ValueError(mode)


class Native:
    def __init__(self, exe: Path, mode: str, n: int, directory: Path):
        self.path = directory / "values.bin"
        self.process = subprocess.Popen([str(exe), mode, str(n), str(self.path)],
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        text=True)
        line = self.process.stdout.readline().split()
        if len(line) != 2 or line[0] != "READY":
            self.close()
            raise RuntimeError("Native driver failed during setup; see diagnostic above")
        self.clock_rate = int(line[1])

    def command(self, repetitions: int) -> str:
        self.process.stdin.write(f"{repetitions}\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline().strip()
        if not line:
            raise RuntimeError("Native driver exited unexpectedly; see diagnostic above")
        return line

    def validate(self, n: int):
        if self.command(0) != "VALIDATED":
            raise RuntimeError("Invalid validation response")
        return np.fromfile(self.path, dtype=np.complex128).reshape((n, 2), order="F")

    def batch(self, count: int):
        elapsed, checksum = map(float, self.command(count).split())
        if not np.isfinite(elapsed + checksum) or elapsed < 0:
            raise RuntimeError("Invalid native timing/checksum")
        return elapsed, checksum

    def close(self):
        if self.process.poll() is None:
            self.process.stdin.close()
        self.process.wait(timeout=30)
        self.process.stdout.close()


def python_batch(calculate, count: int):
    checksum = 0.
    start = time.perf_counter()
    for repeat in range(count):
        result = calculate(repeat % 2)
        checksum += float((result[0] + result[-1]).real)
    elapsed = time.perf_counter() - start
    if not np.isfinite(elapsed + checksum):
        raise RuntimeError("Invalid Python timing/checksum")
    return elapsed, checksum


def calibrate(batch, target: float) -> int:
    count = 1
    while True:
        elapsed, _ = batch(count)
        if elapsed >= target:
            return count
        if count >= 2**26:
            raise RuntimeError("Clock/calibration did not produce a measurable batch")
        count *= 2


def validate(actual, calculate, mode, inputs, exact):
    expected = np.column_stack([np.array(calculate(k), copy=True) for k in range(2)])
    if not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)):
        raise RuntimeError("Nonfinite validation result")
    error = float(np.max(np.abs(actual - expected)))
    # Relative-only CDF checking also checks the small lower-tail probabilities.
    if mode.startswith("normal_cdf"):
        np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=0)
    elif mode == "normal_ppf_array":
        np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-13)
    elif mode.startswith("fft"):
        if error > 3e-12 * max(1., float(np.max(np.abs(expected)))):
            raise RuntimeError(f"FFT validation error: {error}")
    else:
        np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)
    residual = None
    if mode.startswith("lu"):
        matrix, rhs = inputs
        residual = 0.
        for k in range(2):
            np.testing.assert_allclose(actual[:, k].real, exact[k], rtol=2e-12, atol=2e-13)
            r = np.linalg.norm(matrix @ actual[:, k].real - rhs[k], np.inf)
            scale = (np.linalg.norm(matrix, np.inf) * np.linalg.norm(actual[:, k], np.inf)
                     + np.linalg.norm(rhs[k], np.inf))
            residual = max(residual, float(r / scale))
        if residual > 5e-13:
            raise RuntimeError(f"Linear-solve residual: {residual}")
    return {"max_absolute_error": error, "scaled_solve_residual": residual,
            "checked_values": int(actual.size)}, expected


def check_checksum(actual: float, expected, count: int):
    terms = (expected[0, :] + expected[-1, :]).real
    reference = (count//2)*(terms[0] + terms[1]) + (count % 2)*terms[0]
    # Allow accumulation roundoff, while detecting missing/wrong timed work.
    tolerance = 1e-10 * count * max(1., float(np.max(np.abs(terms))))
    if abs(actual - reference) > tolerance:
        raise RuntimeError(f"Timed checksum mismatch: {actual} versus {reference}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--case", action="append", choices=[c[0] for c in CASES])
    parser.add_argument("--quick", action="store_true", help="small smoke workloads; not representative timings")
    parser.add_argument("--check-only", action="store_true", help="runtime-checked build and agreement checks, no timings")
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--target-seconds", type=float, default=.1, help="minimum calibrated batch duration")
    parser.add_argument("--output", type=Path, default=ROOT / "build/scipy-speed/results.json")
    args = parser.parse_args()
    if args.repeats < 3 or not np.isfinite(args.target_seconds) or args.target_seconds <= 0:
        parser.error("require repeats >= 3 and finite target-seconds > 0")
    exe, build_info = build(args.compiler, args.check_only)
    config = io.StringIO()
    with contextlib.redirect_stdout(config):
        np.show_config()
        scipy.show_config()
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True)
    report = {"build": build_info, "python": sys.version, "numpy": np.__version__,
              "scipy": scipy.__version__, "platform": platform.platform(),
              "machine": platform.machine(), "processor": platform.processor(),
              "git_revision": revision.stdout.strip(), "native_library_config": config.getvalue(),
              "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "started_utc": datetime.now(timezone.utc).isoformat(),
              "logical_cpu_count": os.cpu_count(), "requested_repeats": args.repeats,
              "target_batch_seconds": args.target_seconds,
              "thread_environment": {v: os.environ[v] for v in THREAD_ENV},
              "python_clock": vars(time.get_clock_info("perf_counter")),
              "quick": args.quick, "check_only": args.check_only, "cases": []}
    print("case                       n    SciFort ms     SciPy ms   SciPy/SciFort", flush=True)
    for label, mode, normal_size, quick_size in CASES:
        if args.case and label not in args.case:
            continue
        n = quick_size if args.quick else normal_size
        calculate, inputs, exact = prepare(mode, n)
        temporary = tempfile.TemporaryDirectory(dir=exe.parent)
        native = Native(exe, mode, n, Path(temporary.name))
        try:
            accuracy, expected = validate(native.validate(n), calculate, mode, inputs, exact)
            row = {"case": label, "n": n, "accuracy": accuracy, "native_clock_rate": native.clock_rate}
            if args.check_only:
                print(f"{label:25} {n:6}  agreement passed", flush=True)
            else:
                py_batch = lambda count: python_batch(calculate, count)
                # Warm up both input sets, then calibrate independently per backend.
                native.batch(2)
                py_batch(2)
                nf = calibrate(native.batch, args.target_seconds)
                npy = calibrate(py_batch, args.target_seconds)
                ftimes, ptimes = [], []
                for sample in range(args.repeats):
                    # Alternate backend order to reduce systematic drift bias.
                    if sample % 2:
                        pt, pc = py_batch(npy)
                        ft, fc = native.batch(nf)
                    else:
                        ft, fc = native.batch(nf)
                        pt, pc = py_batch(npy)
                    if ft <= 0 or pt <= 0:
                        raise RuntimeError("Clock resolution too coarse; increase target-seconds")
                    check_checksum(fc, expected, nf)
                    check_checksum(pc, expected, npy)
                    # Each checksum also consumes the result of every timed call.
                    ftimes.append(ft/nf)
                    ptimes.append(pt/npy)
                fmedian, pmedian = statistics.median(ftimes), statistics.median(ptimes)
                row.update(scifort_seconds=ftimes, scipy_seconds=ptimes,
                           scifort_batch_size=nf, scipy_batch_size=npy,
                           scifort_median_seconds=fmedian, scipy_median_seconds=pmedian,
                           scipy_over_scifort=pmedian/fmedian, last_checksums=[fc, pc])
                print(f"{label:25} {n:6} {fmedian*1e3:12.5f} {pmedian*1e3:12.5f} {pmedian/fmedian:15.3f}", flush=True)
            report["cases"].append(row)
        finally:
            native.close()
            temporary.cleanup()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Report: {args.output}")


if __name__ == "__main__":
    main()
