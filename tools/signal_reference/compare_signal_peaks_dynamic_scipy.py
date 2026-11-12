#!/usr/bin/env python3
"""Compare per-position SciFort peak intervals against SciPy 1.17.

Pass path to the compiled peak_dynamic_bounds_driver as the first argument.
The seeded driver covers all ten interval endpoints and mixed conditions.
"""
from __future__ import annotations

import sys
import warnings
import subprocess
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

FIELDS = (
    "peaks", "peak_heights", "left_thresholds", "right_thresholds",
    "prominences", "left_bases", "right_bases", "widths", "width_heights",
    "left_ips", "right_ips", "plateau_sizes", "left_edges", "right_edges",
)
INT_FIELDS = {"peaks", "left_bases", "right_bases", "plateau_sizes", "left_edges", "right_edges"}


def make_case(rng: np.random.Generator, case: int):
    n = int(rng.integers(3, 135)) if case % 29 else int(rng.integers(135, 201))
    kind = case % 6
    if kind == 0:
        x = rng.standard_normal(n) * 1.8
    elif kind == 1:
        x = np.sin(np.linspace(0, rng.uniform(4, 45), n)) * 2.5 + rng.normal(0, 0.15, n)
    elif kind == 2:
        x = rng.integers(-4, 5, n).astype(float)
    elif kind == 3:
        x = rng.uniform(-0.7, 0.7, n)
        for _ in range(max(1, n // 8)):
            x[int(rng.integers(n))] += float(rng.uniform(2, 6))
    elif kind == 4:
        x = rng.uniform(-2, 3, n)
        for _ in range(max(1, n // 14)):
            pos = int(rng.integers(1, n - 1))
            length = int(rng.integers(1, min(7, n - pos) + 1))
            x[pos:pos + length] = x[pos]
    else:
        x = np.full(n, float(rng.integers(-2, 3)))
        if n > 5:
            x[n // 3: n // 3 + 2] += 3.0
            x[2 * n // 3] += 5.0
    # SciPy's np.argsort may break equal-height *distinct* peaks differently
    # from SciFort's deterministic tie policy. Preserve flat-top runs, but
    # make heights of distinct runs unique before testing distance > 1.
    if case % 5 != 0:
        starts = np.r_[0, np.flatnonzero(np.diff(x) != 0) + 1]
        ends = np.r_[starts[1:], n]
        for start, end in zip(starts, ends):
            x[start:end] += 1e-8 * int(start)
    bounds = [
        rng.uniform(-3, 3, n), rng.uniform(-1, 6, n),
        rng.uniform(0, 2.5, n), rng.uniform(0, 5, n),
        rng.uniform(0, 5, n), rng.uniform(0, 6, n),
        rng.uniform(0, 3, n), rng.uniform(0, 7, n),
        rng.uniform(0, 3, n), rng.uniform(0, 5, n),
    ]
    if case % 7 == 0:
        # Nonrestrictive endpoints ensure surviving peaks and properties.
        for i in (0, 2, 4, 6, 8):
            bounds[i][:] = -4.0 if i == 0 else 0.0
        for i in (1, 3, 5, 7, 9):
            bounds[i][:] = 100.0
    if case % 20 == 0:
        # Inclusive boundaries, including plateau lengths and center heights.
        bounds[0] = x.copy()
        bounds[1] = x.copy()
    if case % 5 == 0:
        mask = (1 << 10) - 1
    elif case % 5 == 1:
        mask = 1 << ((case // 5) % 10)
    elif case % 5 == 2:
        mask = int(rng.integers(1, 1 << 10))
    elif case % 5 == 3:
        mask = 0
    else:
        mask = (3 << ((case % 5) * 2))
    distance = (1, 2, 3, 5, 9)[case % 5]
    window = (0, 2, 3, 6, 13, 32)[case % 6]
    rel_height = (0.0, 0.05, 0.5, 1.0, 1.6)[case % 5]
    return x, bounds, mask, distance, window, rel_height


def compare(exe: Path, cases: int = 1600) -> tuple[int, int, float]:
    rng = np.random.default_rng(20261025)
    inputs = []
    references = []
    for icase in range(cases):
        x, bounds, mask, distance, window, rel_height = make_case(rng, icase)
        inputs.append(f"{len(x)} {mask} {distance} {window} {rel_height:.17g}\n")
        inputs.append(" ".join(f"{v:.17g}" for v in x) + "\n")
        for b in bounds:
            inputs.append(" ".join(f"{v:.17g}" for v in b) + "\n")
        intervals = []
        for i in range(0, 10, 2):
            imin = bounds[i] if mask & (1 << i) else None
            imax = bounds[i + 1] if mask & (1 << (i + 1)) else None
            intervals.append((imin, imax) if imin is not None or imax is not None else None)
        args = {
            "height": intervals[0], "threshold": intervals[1],
            "prominence": intervals[2], "width": intervals[3],
            "plateau_size": intervals[4], "distance": distance,
            "rel_height": rel_height, "wlen": window or None,
        }
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            peaks, props = find_peaks(x, **args)
        references.append((peaks, props, mask, x))
    out = subprocess.run([str(exe.resolve())], input="".join(inputs), text=True,
                         capture_output=True, check=True)
    lines = out.stdout.splitlines()
    assert len(lines) == 15 * cases, (len(lines), cases)
    total_values = 0
    max_err = 0.0
    for icase, (peaks, props, mask, x) in enumerate(references):
        segment = lines[15 * icase:15 * (icase + 1)]
        status, n = [int(z) for z in segment[0].split()]
        assert status == 0, (icase, status, mask)
        assert n == len(peaks), (icase, n, len(peaks), mask)
        for j, name in enumerate(FIELDS):
            dtype = int if name in INT_FIELDS else float
            actual = np.fromstring(segment[j + 1], sep=" ", dtype=dtype)
            expected = peaks if name == "peaks" else props.get(name, np.array([], dtype=dtype))
            if name in {"peaks", "left_bases", "right_bases", "left_edges", "right_edges", "left_ips", "right_ips"}:
                expected = expected + 1
            if len(actual) != len(expected):
                raise AssertionError((icase, name, mask, len(actual), len(expected), x))
            if dtype is int:
                if not np.array_equal(actual, expected):
                    raise AssertionError((icase, name, mask, actual, expected, x))
            elif len(actual):
                err = float(np.max(np.abs(actual - expected)))
                max_err = max(max_err, err)
                if not np.allclose(actual, expected, rtol=1e-10, atol=1e-10):
                    raise AssertionError((icase, name, mask, actual, expected, err, x))
            total_values += len(actual)
    return cases, total_values, max_err


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./peak_dynamic_bounds_driver")
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 1600
    cases, values, max_err = compare(target, count)
    print(f"PASS: {cases} independent SciPy 1.17 peak cases; {values} values; "
          f"maximum absolute difference {max_err:.3g}")
