#!/usr/bin/env python3
"""Compare SciFort bilinear transforms and elliptic prototypes with SciPy."""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy import signal
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ["kinds", "constants", "math", "special/bessel_i", "special/bessel_k",
           "special/bessel_jy", "special/bessel_order", "special/bessel_derivatives",
           "special/hankel", "special/bessel_complex", "signal_elliptic",
           "signal_analog_prototypes",
           "signal_frequency_transforms", "signal_bilinear", "signal_ellipap"]


def sorted_root_error(actual, expected):
    assert len(actual) == len(expected)
    if not len(actual):
        return 0.0
    costs = np.abs(actual[:, None] - expected[None, :])
    row, col = linear_sum_assignment(costs)
    return float(np.max(costs[row, col] / np.maximum(1.0, np.abs(expected[col]))))


def main():
    rng = np.random.default_rng(20261024)
    inputs = []
    cases = []

    def put_array(a):
        for v in a:
            inputs.append(f"{v.real:.17g} {v.imag:.17g}")

    # All-pass, strictly proper, complex, and high-order coefficient tests.
    for n in range(10):
        for d in range(0, n + 1):
            for fs in (0.3, 1.0, 7.5):
                b = rng.normal(size=d + 1) + 1j*rng.normal(size=d + 1)
                a = rng.normal(size=n + 1) + 1j*rng.normal(size=n + 1)
                a[0] += 3.0
                inputs.append("B")
                inputs.append(f"{len(b)} {len(a)} {fs:.17g}")
                put_array(b)
                put_array(a)
                cases.append(("B", b, a, fs))

    for n in range(9):
        for nz in range(n + 1):
            for fs in (0.2, 1.5, 24.0):
                z = -0.5-rng.random(nz) + 1j*rng.normal(size=nz)
                p = -1-rng.random(n) + 1j*rng.normal(size=n)
                k = 0.7 + 0.3j
                inputs.append("Z")
                inputs.append(f"{len(z)} {len(p)} {fs:.17g} {k.real:.17g} {k.imag:.17g}")
                put_array(z)
                put_array(p)
                cases.append(("Z", z, p, k, fs, False))

    # Real filters with conjugate-paired roots also match SciPy's gain rule.
    for n in (2, 4, 6, 8):
        for nz in range(0, n + 1, 2):
            for fs in (0.8, 8.0):
                pp = -1.0-rng.random(n//2) + 1j*(0.3+rng.random(n//2))
                zz = -0.2-rng.random(nz//2) + 1j*(0.4+rng.random(nz//2))
                p = np.concatenate((pp, pp.conj()))
                z = np.concatenate((zz, zz.conj()))
                k = 1.3 + 0.0j
                inputs.append("Z")
                inputs.append(f"{len(z)} {len(p)} {fs:.17g} {k.real:.17g} {k.imag:.17g}")
                put_array(z)
                put_array(p)
                cases.append(("Z", z, p, k, fs, True))

    for n in range(0, 21):
        for rp, rs in ((0.1, 30.), (0.5, 50.), (2., 65.), (4., 85.)):
            inputs.append("E")
            inputs.append(f"{n} {rp} {rs}")
            cases.append(("E", n, rp, rs))

    for n in (11, 12, 13, 14, 16, 18, 20):
        for norm in ("phase", "delay", "mag"):
            inputs.append("S")
            inputs.append(f"{n} {norm}")
            cases.append(("S", n, norm))

    with tempfile.TemporaryDirectory(prefix="scifort-bilinear-") as temp:
        work = Path(temp)
        objs = []
        for name in SOURCES:
            obj = work / f"{name.replace(chr(47), chr(95))}.o"
            subprocess.run(["gfortran", "-std=f2018", "-O0", "-Wimplicit-interface",
                            "-Werror=implicit-interface", "-J", str(work), "-I", str(work),
                            "-c", str(ROOT / "src" / (name + ".f90")), "-o", str(obj)],
                           check=True, cwd=work)
            objs.append(str(obj))
        exe = work / "reference"
        subprocess.run(["gfortran", "-std=f2018", "-O0", "-I", str(work),
                        str(ROOT / "tools/signal_reference/driver_signal_bilinear_ellipap.f90"),
                        *objs, "-o", str(exe)], check=True, cwd=work)
        proc = subprocess.run([str(exe)], input="\n".join(inputs) + "\n",
                              text=True, capture_output=True, check=True, cwd=work)

    lines = iter(proc.stdout.splitlines())
    maxdiff = {"B": 0.0, "Z": 0.0, "E": 0.0, "S": 0.0}
    elements = 0
    for case in cases:
        kind = case[0]
        name, status, nz, n_poles = next(lines).split()
        assert name == kind and int(status) == 0, (case, name, status)
        nz, n_poles = int(nz), int(n_poles)
        def get_array(count):
            return np.array([complex(*map(float, next(lines).split()))
                             for _ in range(count)], dtype=complex)
        if kind == "B":
            _, b, a, fs = case
            assert (nz, n_poles) == (max(len(b), len(a)),) * 2
            result_b, result_a = get_array(nz), get_array(n_poles)
            expected_b, expected_a = signal.bilinear(b, a, fs=fs)
            err = max(np.max(np.abs(result_b - expected_b) / np.maximum(1., np.abs(expected_b))),
                      np.max(np.abs(result_a - expected_a) / np.maximum(1., np.abs(expected_a))))
            elements += nz + n_poles
            tolerance = 1e-8
        elif kind == "Z":
            _, z, p, k, fs, real_system = case
            gain = complex(*map(float, next(lines).split()))
            zz, pp = get_array(nz), get_array(n_poles)
            zref, pref, kref = signal.bilinear_zpk(z, p, k, fs=fs)
            err = max(sorted_root_error(zz, zref), sorted_root_error(pp, pref),
                      abs(gain - k * np.prod(2*fs-z)/np.prod(2*fs-p)) /
                      max(1., abs(k * np.prod(2*fs-z)/np.prod(2*fs-p))))
            if real_system:
                err = max(err, abs(gain-kref)/max(1.0, abs(kref)))
            elements += nz + n_poles + 1
            tolerance = 4e-12
        elif kind in ("E", "S"):
            gain = float(next(lines))
            zz, pp = get_array(nz), get_array(n_poles)
            if kind == "E":
                _, n, rp, rs = case
                zref, pref, kref = signal.ellipap(n, rp, rs)
            else:
                _, n, norm = case
                zref, pref, kref = signal.besselap(n, norm)
            err = max(sorted_root_error(zz, np.atleast_1d(zref)),
                      sorted_root_error(pp, np.atleast_1d(pref)),
                      abs(gain-kref)/max(1., abs(kref)))
            elements += nz + n_poles + 1
            tolerance = 2e-7 if kind == "E" else 3e-8
        else:
            raise AssertionError(case)
        maxdiff[kind] = max(maxdiff[kind], err)
        assert err < tolerance, (case, err, tolerance)
    assert next(lines, None) is None
    print(f"PASS: {len(cases)} cases, {elements} coefficients/roots/gains")
    for name, err in maxdiff.items():
        print(f"{name}: maximum relative/scaled error = {err:.3e}")


if __name__ == "__main__":
    main()
