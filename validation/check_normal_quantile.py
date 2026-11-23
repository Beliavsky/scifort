#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check direct normal quantiles against SciPy and independent mpmath values."""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("speed", ROOT / "benchmarks/compare_scipy.py")
speed = importlib.util.module_from_spec(spec)
spec.loader.exec_module(speed)
import mpmath as mp
import numpy as np
from scipy import special
from scipy.stats import rankdata


def probabilities():
    rng = np.random.default_rng(20261122)
    bits = rng.integers(1, np.float64(1).view(np.uint64), size=10000, dtype=np.uint64)
    values = [np.linspace(0, 1, 20001), np.logspace(-323, -1, 5000), bits.view(np.float64),
              np.array([0., np.nextafter(0., 1.), np.finfo(float).tiny, np.nextafter(1., 0.), 1.])]
    for switch in [np.exp(-32.), 0.13533528323661269189, .5, 1-0.13533528323661269189]:
        neighbors = [switch]
        low, high = switch, switch
        for _ in range(32):
            low = np.nextafter(low, 0.)
            high = np.nextafter(high, 1.)
            neighbors += [low, high]
        values.append(np.array(neighbors))
    lower = np.concatenate(values)
    return np.unique(np.concatenate([lower, 1-lower]))


def check_fligner_reference():
    groups = [np.array(v) for v in [
        [8.88, 9.12, 9.04, 8.98, 9., 9.08, 9.01, 8.85, 9.06, 8.99],
        [8.88, 8.95, 9.29, 9.44, 9.15, 9.58, 8.36, 9.18, 8.67, 9.05],
        [8.95, 9.12, 8.95, 8.85, 9.03, 8.84, 9.07, 8.98, 8.86, 8.98]]]
    deviations = np.concatenate([abs(v-np.median(v)) for v in groups])
    p = .5 + rankdata(deviations)/(2*(deviations.size+1))
    with mp.workdps(100):
        scores = [mp.sqrt(2)*mp.erfinv(2*mp.mpf(float(v))-1) for v in p]
        mean = mp.fsum(scores)/len(scores)
        variance = mp.fsum((v-mean)**2 for v in scores)/(len(scores)-1)
        numerator = mp.fsum(10*(mp.fsum(scores[k:k+10])/10-mean)**2 for k in [0, 10, 20])
        reference = numerator/variance
        text = (ROOT / "test/test_hypothesis_extended.f90").read_text()
        stored = re.search(r"'fligner median statistic', flig%statistic, ([0-9.]+)_dp", text)
        assert stored is not None
        assert abs(mp.mpf(float(stored[1]))-reference) < mp.mpf("3e-15")
        print("Independent Fligner reference:", mp.nstr(reference, 45))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=["checked", "optimized"], default="checked")
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    _, metadata = speed.build(args.compiler, args.profile == "checked")
    directory = ROOT / "build/scipy-speed" / metadata["source_sha256"][:16]
    objects = [directory / (name + ".o") for name in ["kinds", "constants", "math", "normal"]]
    executables = []
    for source in ["validation/normal_quantile_driver.f90", "test/test_normal_quantile.f90"]:
        executable = directory / (Path(source).stem + (".exe" if sys.platform == "win32" else ""))
        subprocess.run([args.compiler, *metadata["flags"], "-I" + str(directory),
                        str(ROOT / source), *map(str, objects), "-o", str(executable)], check=True)
        executables.append(executable)
    subprocess.run([str(executables[1])], check=True)
    p = probabilities()
    payload = str(p.size) + "\n" + "\n".join(format(v, ".17e") for v in p) + "\n"
    result = subprocess.run([str(executables[0])], input=payload, text=True,
                            capture_output=True, check=True)
    actual = np.loadtxt(result.stdout.splitlines())
    target = special.ndtri(p)
    np.testing.assert_array_equal(actual[:, 1], -actual[:, 0])
    np.testing.assert_array_equal(np.isneginf(actual[:, 0]), p == 0)
    np.testing.assert_array_equal(np.isposinf(actual[:, 0]), p == 1)
    finite = (p > 0) & (p < 1)
    np.testing.assert_allclose(actual[finite, 0], target[finite], rtol=4e-15, atol=0)
    # Threshold-neighbor ties can differ by a few ulps; dense-grid strict
    # monotonicity is separately checked by the native regression.
    spacing = np.abs(np.spacing(target[finite]))
    maximum_ulps = float(np.max(np.abs(actual[finite, 0]-target[finite])/spacing))

    selected = np.unique(np.r_[np.flatnonzero(finite)[::max(1, sum(finite)//160)],
                                np.flatnonzero((p > .4999999999999) & (p < .5000000000001))])
    maximum_relative = mp.mpf(0)
    with mp.workdps(400):
        for index in selected:
            reference = mp.sqrt(2)*mp.erfinv(2*mp.mpf(float(p[index]))-1)
            value = mp.mpf(float(actual[index, 0]))
            if reference == 0:
                assert value == 0
            else:
                relative = abs((value-reference)/reference)
                maximum_relative = max(maximum_relative, relative)
                assert relative < mp.mpf("3e-15"), (float(p[index]), relative)
    print(f"{args.profile}: {p.size} SciPy comparisons, {len(selected)} independent "
          f"400-digit comparisons; max SciPy ulps={maximum_ulps:g}, "
          f"max independent relative error={float(maximum_relative):.6g}")
    check_fligner_reference()


if __name__ == "__main__":
    main()
