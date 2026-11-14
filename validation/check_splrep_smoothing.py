#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare native adaptive scalar FITPACK smoothing with SciPy."""
import argparse
from pathlib import Path
import subprocess
import re

import numpy as np
import scipy
from scipy.interpolate import splrep, splev, splint


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--upstream-tests", type=Path,
                        help="Pinned perazz fitpack_curve_tests.f90 for its fpknot crash data")
    parser.add_argument("--upstream-executable", type=Path,
                        help="Original perazz core linked with splrep_upstream_driver.f90")
    args = parser.parse_args()
    if args.upstream_tests and not args.upstream_executable:
        parser.error("The upstream crash fixture requires the upstream reference executable")
    rng = np.random.default_rng(20261113)
    cases = []
    for k in range(1, 6):
        for i in range(30):
            n = k+1 if i == 0 else 25+3*i
            x = np.cumsum(rng.uniform(0.4, 1.5, n))
            x = (x-x[0])/(x[-1]-x[0])*4-1
            y = np.sin((1+i % 5)*x)+0.12*rng.normal(size=n)
            w = np.exp(rng.uniform(-1.2, 1.2, n))
            bounds = (x[0]-(0.3 if i % 3 == 0 else 0),
                      x[-1]+(0.2 if i % 3 == 0 else 0))
            target = sum((w*y)**2)*(0.0001, 0.01, 0.1, 2, 0.001)[i % 5]
            q = np.linspace(*bounds, 47)
            cases.append((k, x, y, w, bounds, target, q))
    if args.upstream_tests:
        source = args.upstream_tests.read_text(encoding="utf-8")
        fixture = source.split("function test_fpknot_crash()", 1)[1]
        fixture = fixture.split("end function test_fpknot_crash", 1)[0]
        literals = re.search(r"\by\(\*\)\s*=\s*\[([^]]+)\]", fixture).group(1)
        # The inspected constructor uses default-real literals before binary64 assignment.
        y = np.array([float(v) for v in literals.replace("&", "").split(",")],
                     dtype=np.float32).astype(float)
        x = np.arange(109, dtype=float)
        assert len(y) == len(x)
        for target in (0.1, 1.0, 100.0):
            cases.append((1, x, y, np.ones(len(x)), (0.0, 108.0), target,
                          np.linspace(0, 108, 47)))
    lines = []
    for k, x, y, w, bounds, target, q in cases:
        lines.append(f"{len(x)} {k} {len(q)}")
        for values in ([target, *bounds], x, y, w, q):
            lines.append(" ".join(f"{v:.17g}" for v in values))
    result = subprocess.run([str(args.executable.resolve())],
                            input="\n".join(lines)+"\n", capture_output=True,
                            text=True, check=True)
    rows = iter(result.stdout.splitlines())
    reference_rows = None
    reference_name = f"SciPy {scipy.__version__}"
    if args.upstream_executable:
        reference = subprocess.run([str(args.upstream_executable.resolve())],
                                   input="\n".join(lines)+"\n", capture_output=True,
                                   text=True, check=True)
        reference_rows = iter(reference.stdout.splitlines())
        reference_name = "original perazz/FITPACK"
    count = adaptive = polynomial = 0
    maximum = residual_maximum = 0.0
    for index, (k, x, y, w, bounds, target, q) in enumerate(cases):
        status = int(next(rows))
        assert status == 0, (index, k, status)
        nt = int(next(rows))
        knots = np.fromstring(next(rows), sep=" ")
        coeff = np.fromstring(next(rows), sep=" ")
        actual = np.fromstring(next(rows), sep=" ")
        if reference_rows is None:
            ref = splrep(x, y, w=w, k=k, s=target, xb=bounds[0], xe=bounds[1])
            ref_actual = splev(x, ref)
            ref_queries = [splev(q, ref, der=d) for d in range(min(3, k)+1)]
            ref_integral = splint(*bounds, ref)
        else:
            assert int(next(reference_rows)) == 0, (index, "upstream fit failure")
            assert int(next(reference_rows)) == nt, (index, "upstream knot count")
            ref = (np.fromstring(next(reference_rows), sep=" "),
                   np.fromstring(next(reference_rows), sep=" "), k)
            ref_actual = np.fromstring(next(reference_rows), sep=" ")
            ref_queries = [np.fromstring(next(reference_rows), sep=" ")
                           for _ in range(min(3, k)+1)]
            ref_integral = float(next(reference_rows))
        assert nt == len(knots) == len(coeff) == len(ref[0]), (index, nt, len(ref[0]))
        np.testing.assert_allclose(knots, ref[0], rtol=2e-14, atol=2e-14)
        fp = sum((w*(actual-y))**2)
        assert fp <= target*1.0011+1e-24, (index, fp, target)
        if nt > 2*(k+1):
            adaptive += 1
            residual_maximum = max(residual_maximum, abs(fp-target)/target)
        else:
            polynomial += 1
        comparisons = [(coeff, ref[1]), (actual, ref_actual)]
        for order in range(min(3, k)+1):
            values = np.fromstring(next(rows), sep=" ")
            comparisons.append((values, ref_queries[order]))
        comparisons.append((np.array([float(next(rows))]),
                            np.array([ref_integral])))
        for actual, expected in comparisons:
            error = np.max(abs(actual-expected)/(1+abs(expected)))
            assert error < 3e-8, (index, k, error)
            maximum = max(maximum, error)
            count += actual.size
    assert next(rows, None) is None
    if reference_rows is not None:
        assert next(reference_rows, None) is None
    assert adaptive > 75 and polynomial > 20
    print(f"PASS: {len(cases)} fits, {adaptive} adaptive and {polynomial} polynomial; "
          f"{count} {reference_name} coefficient/value/derivative/integral comparisons")
    print(f"Maximum normalized difference {maximum:.6g}; "
          f"maximum adaptive residual relative error {residual_maximum:.6g}")


if __name__ == "__main__":
    main()
