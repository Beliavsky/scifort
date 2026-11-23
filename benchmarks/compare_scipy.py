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
from scipy import fft, interpolate, linalg, special, stats

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ["src/kinds.f90", "src/constants.f90", "src/math.f90",
           "src/stats/normal.f90",
           *[f"src/stats/{name}.f90" for name in
             ("exponential", "logistic", "laplace", "cauchy", "weibull", "lognormal")],
           "src/linalg.f90", "src/linalg_dense.f90",
           "src/signal_windows_extra.f90", "src/signal_dpss.f90",
           "src/signal_spectral.f90", "src/interpolate_barycentric.f90"]
DRIVER = "benchmarks/scipy_speed_driver.f90"
# label, native mode, ordinary size, smoke-check size
FAMILIES = ("normal", "exponential", "logistic", "laplace", "cauchy", "weibull", "lognormal")
METHODS = ("pdf", "logpdf", "cdf", "sf", "logcdf", "logsf", "ppf", "isf")
CASES = [(f"{family}_{method}_array", f"{family}_{method}_array",
          65536 if family == "normal" and method == "cdf" else 16384, 128)
         for family in FAMILIES for method in METHODS]
CASES += [
    ("normal_cdf_scalar", "normal_cdf_scalar", 256, 16),
    *[(f"{direction}_{size}_{policy}", f"{direction}_{policy}", n, q)
      for direction in ("fft", "ifft")
      for size, n, q in (("power2", 4096, 64), ("mixed", 4050, 90), ("prime", 4093, 67))
      for policy in ("reuse", "setup")],
    *[(mode, mode, 128, 16) for mode in
      ("lu_factor_solve", "lu_reuse", "lu_transpose", "cho_factor_solve", "cho_reuse",
       "triangular_lower", "triangular_upper", "triangular_transpose")],
    *[(mode, mode, 4096, 64) for mode in
      ("barycentric", "barycentric_setup", "barycentric_derivative")],
]


def category(mode):
    if mode.split("_")[0] in FAMILIES:
        return "distribution"
    if mode.startswith(("fft", "ifft")):
        return "fft"
    if mode.startswith("barycentric"):
        return "interpolation"
    return "linear_solve"



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
    if category(mode) == "distribution":
        family, method, _ = mode.split("_")
        if method in ("ppf", "isf"):
            inputs = 1e-6 + (1 - 2e-6)*t
        elif family in ("exponential", "weibull", "lognormal"):
            inputs = .01 + 8*t
        else:
            inputs = -8 + 16*t
        inputs = np.asfortranarray(inputs)
        distribution = getattr(stats, {"normal": "norm", "exponential": "expon",
                                     "weibull": "weibull_min", "lognormal": "lognorm"}.get(family, family))
        shape = (1.7,) if family == "weibull" else ((.8,) if family == "lognormal" else ())
        function = getattr(distribution, method)
        # Retain the original direct-kernel normal throughput comparisons.
        direct = family == "normal" and method in ("cdf", "ppf")
        if direct:
            function = special.ndtr if method == "cdf" else special.ndtri

        def calculate(k):
            if mode == "normal_cdf_scalar":
                for i in range(n):
                    values[i] = special.ndtr(inputs[i, k])
                return values
            if direct:
                return function(inputs[:, k], out=values)
            return function(inputs[:, k], *shape)

        return calculate, inputs, None
    if mode.startswith(("fft", "ifft")):
        base = np.arange(n, dtype=np.float64) / n
        inputs = [np.sin(17*base) + .1*k + 1j*np.cos(11*base) for k in range(2)]

        def calculate(k):
            return (fft.ifft if mode.startswith("ifft") else fft.fft)(inputs[k], norm="backward", overwrite_x=False, workers=1)

        return calculate, None, None
    if category(mode) == "linear_solve":
        index = np.arange(1, n+1, dtype=np.float64)
        matrix = np.asfortranarray(1/(index[:, None] + index[None, :]) + n*np.eye(n))
        lower = mode == "triangular_lower"
        transpose = mode in ("lu_transpose", "triangular_transpose")
        if mode.startswith("triangular"):
            matrix = np.asfortranarray(np.tril(matrix) if lower else np.triu(matrix))
        operator = matrix.T if transpose else matrix
        exact = [np.sin(index) + .1*k for k in range(2)]
        rhs = [operator @ x for x in exact]
        factors = (linalg.cho_factor(matrix, lower=True) if mode.startswith("cho")
                   else linalg.lu_factor(matrix))

        def calculate(k):
            if mode.startswith("triangular"):
                return linalg.solve_triangular(matrix, rhs[k], lower=lower,
                                               trans="T" if transpose else "N", check_finite=True)
            if mode.startswith("cho"):
                packed = linalg.cho_factor(matrix, lower=True) if mode == "cho_factor_solve" else factors
                return linalg.cho_solve(packed, rhs[k], check_finite=True)
            packed = linalg.lu_factor(matrix) if mode == "lu_factor_solve" else factors
            return linalg.lu_solve(packed, rhs[k], trans=1 if transpose else 0, check_finite=True)

        return calculate, (operator, rhs), exact
    if mode.startswith("barycentric"):
        nodes = np.cos(np.pi*np.arange(33)/32)
        samples = 1/(1+25*nodes**2)
        polynomial = interpolate.BarycentricInterpolator(nodes, samples, rng=20261122)
        queries = np.asfortranarray(-.99 + 1.98*t)

        def calculate(k):
            fitted = (interpolate.BarycentricInterpolator(nodes, samples, rng=20261122)
                      if mode == "barycentric_setup" else polynomial)
            if mode == "barycentric_derivative":
                return fitted.derivative(queries[:, k], der=1)
            return fitted(queries[:, k])

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
        endpoints = result[0] + result[-1]
        checksum += float(endpoints.real + endpoints.imag)
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
    elif category(mode) == "distribution" and mode.split("_")[1] in ("pdf", "cdf", "sf"):
        np.testing.assert_allclose(actual, expected, rtol=3e-12, atol=0)
    elif category(mode) == "distribution":
        np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-13)
    elif mode.startswith(("fft", "ifft")):
        if error > 3e-12 * max(1., float(np.max(np.abs(expected)))):
            raise RuntimeError(f"FFT validation error: {error}")
    else:
        np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)
    residual = None
    if category(mode) == "linear_solve":
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
    endpoints = expected[0, :] + expected[-1, :]
    terms = endpoints.real + endpoints.imag
    reference = (count//2)*(terms[0] + terms[1]) + (count % 2)*terms[0]
    # Allow accumulation roundoff, while detecting missing/wrong timed work.
    tolerance = 1e-10 * count * max(1., float(np.max(np.abs(terms))))
    if abs(actual - reference) > tolerance:
        raise RuntimeError(f"Timed checksum mismatch: {actual} versus {reference}")


def markdown_report(report):
    """Render the same results as JSON, without inventing missing timing data."""
    lines = ["# SciFort versus SciPy timings", "", f"Run: {report['started_utc']}", "",
             f"Platform: {report['platform']}", "",
             f"Compiler: {report['build']['compiler'].splitlines()[0]}", "",
             f"Flags: `{' '.join(report['build']['flags'])}`", "",
             f"Python {report['python'].split()[0]}; NumPy {report['numpy']}; SciPy {report['scipy']}.", "",
             f"Quick workloads: {report['quick']}; validation only: {report['check_only']}.", "",
             "Times are milliseconds per full operation; ratios above one favor SciFort.",
             "These are application workloads, including documented allocation and API overhead.", "",
             "| Case | n | SciFort ms | SciPy ms | SciPy/SciFort | Max abs error |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for row in report["cases"]:
        timing = (f"{row['scifort_median_seconds']*1000:.6g} | "
                  f"{row['scipy_median_seconds']*1000:.6g} | {row['scipy_over_scifort']:.3f}"
                  if "scifort_seconds" in row else "- | - | -")
        lines.append(f"| {row['case']} | {row['n']} | {timing} | "
                     f"{row['accuracy']['max_absolute_error']:.3g} |")
    lines += ["", "## Individual samples", "",
              "Milliseconds per operation, in sampling order. Full metadata and checks are in the sibling JSON.", ""]
    for row in report["cases"]:
        if "scifort_seconds" in row:
            lines += [f"- `{row['case']}`: SciFort " +
                      ", ".join(f"{x*1000:.6g}" for x in row["scifort_seconds"]) +
                      "; SciPy " + ", ".join(f"{x*1000:.6g}" for x in row["scipy_seconds"]) + "."]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--case", action="append", choices=[c[0] for c in CASES])
    parser.add_argument("--category", action="append",
                        choices=("distribution", "fft", "linear_solve", "interpolation"))
    parser.add_argument("--list", action="store_true", help="list workloads without compiling")
    parser.add_argument("--quick", action="store_true", help="small smoke workloads; not representative timings")
    parser.add_argument("--check-only", action="store_true", help="runtime-checked build and agreement checks, no timings")
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--target-seconds", type=float, default=.1, help="minimum calibrated batch duration")
    parser.add_argument("--output", type=Path, default=ROOT / "build/scipy-speed/results.json")
    args = parser.parse_args()
    if args.repeats < 3 or not np.isfinite(args.target_seconds) or args.target_seconds <= 0:
        parser.error("require repeats >= 3 and finite target-seconds > 0")
    selected = [c for c in CASES if (not args.case or c[0] in args.case)
                and (not args.category or category(c[1]) in args.category)]
    if not selected:
        parser.error("no workloads match the selected cases/categories")
    if args.list:
        for label, mode, n, quick in selected:
            print(f"{label:32} {category(mode):14} n={n} quick={quick}")
        return
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
    for label, mode, normal_size, quick_size in selected:
        n = quick_size if args.quick else normal_size
        calculate, inputs, exact = prepare(mode, n)
        temporary = tempfile.TemporaryDirectory(dir=exe.parent)
        native = Native(exe, mode, n, Path(temporary.name))
        try:
            accuracy, expected = validate(native.validate(n), calculate, mode, inputs, exact)
            row = {"case": label, "category": category(mode), "n": n, "accuracy": accuracy, "native_clock_rate": native.clock_rate}
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
    markdown = args.output.with_suffix(".md")
    markdown.write_text(markdown_report(report), encoding="utf-8")
    print(f"Reports: {args.output} and {markdown}")


if __name__ == "__main__":
    main()
