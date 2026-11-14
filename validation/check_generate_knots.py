#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Compare native knot sequences with SciPy and pinned upstream selection logic.

SciPy 1.15.3 supplies the public nonperiodic oracle. For periodic selection,
execute the exact SciPy v1.18.0 Python control flow with an independent NumPy
cyclic LSQ residual evaluator, since 1.15.3 lacks the periodic public option.
The evaluator canonicalizes cyclic exterior knots before every solve; the
upstream loop otherwise leaves them stale between insertions within a batch.
Upstream source is downloaded only to the ignored build cache, hash checked,
and covered by LICENSES/scipy-generate-knots-BSD-3-Clause.txt.
"""
import ast
import hashlib
from pathlib import Path
import subprocess
import sys
import urllib.request

import numpy as np
import scipy
from scipy.interpolate import BSpline, generate_knots
from scipy.interpolate import _dierckx


def upstream_functions():
    cache = Path("build/generate-knots-upstream")
    cache.mkdir(parents=True, exist_ok=True)
    specs = [
        ("_fitpack_repro.py", "82449ae292bb361675439b422cec6c084f9fbfa6519eabec06b1009c51a818bd"),
        ("_bsplines.py", "b099333d1dc57d0d4ad5300c04cd2767350048cac2b9dee464e40ab469510caa"),
    ]
    namespace = dict(np=np, TOL=0.001, _get_residuals=residuals, add_knot=add_knot)
    for name, digest in specs:
        path = cache / name
        if not path.exists():
            url = "https://raw.githubusercontent.com/scipy/scipy/v1.18.0/scipy/interpolate/" + name
            path.write_bytes(urllib.request.urlopen(url, timeout=60).read())
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != digest:
            raise RuntimeError("Unexpected upstream source hash: " + str(path))
        tree = ast.parse(data)
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in {
                "_generate_knots_impl", "_not_a_knot", "_periodic_knots"
            }:
                node.decorator_list = []
                exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    return namespace


def residuals(x, y, t, k, w, periodic=False):
    if periodic:
        t = t.copy()
        independent = len(t)-2*k-1
        for j in range(1, k+1):
            t[k-j] = t[k-j+independent]-(x[-1]-x[0])
            t[-k+j-1] = t[k+j]+(x[-1]-x[0])
    design = BSpline.design_matrix(x, t, k).toarray()
    if periodic:
        independent = len(t) - 2*k - 1
        folded = np.zeros((len(x)-1, independent))
        for col in range(design.shape[1]):
            folded[:, col % independent] += design[:-1, col]
        design = folded
        active = len(x)-1
    else:
        active = len(x)
    a = design*w[:active, None]
    rhs = y[:active]*w[:active, None]
    c = np.linalg.lstsq(a, rhs, rcond=None)[0]
    err = np.sum((a@c-rhs)**2, axis=1)
    if periodic:
        err = np.r_[err, 0.0]
    return err, err.sum()


def add_knot(x, t, k, errors):
    value = _dierckx.fpknot(x, t, k, errors)
    return np.insert(t, np.searchsorted(t, value), value)


def native(exe, x, y, k, target, w, periodic, nest, bounds):
    head = [len(x), y.shape[1], k, int(periodic), nest or 0, int(w is not None),
            target, *bounds]
    weights = np.ones_like(x) if w is None else w
    rows = [" ".join(map(str, head)), " ".join(map(str, x))]
    rows += [" ".join(map(str, row)) for row in y]
    rows += [" ".join(map(str, weights))]
    result = subprocess.run([str(exe)], input="\n".join(rows)+"\n", text=True,
                            capture_output=True, check=True)
    lines = iter(result.stdout.splitlines())
    if int(next(lines)) != 0:
        raise AssertionError("native initialization failed: " + result.stdout)
    vectors = []
    for header in lines:
        code, size = map(int, header.split())
        if code != 0:
            raise AssertionError("native iteration failed: " + result.stdout)
        if not size:
            break
        vector = np.fromstring(next(lines), sep=" ")
        assert len(vector) == size
        vectors.append(vector)
    return vectors


def main():
    exe = Path(sys.argv[1]).resolve()
    upstream = upstream_functions()
    rng = np.random.default_rng(20261116)
    cases = comparisons = scalar_cases = periodic_cases = corrected_exteriors = 0
    residual_checks = site_scale_cases = 0
    max_error = 0.0
    for k in range(1, 6):
        for periodic in (False, True):
            for channels in (1, 3):
                for mode in range(7):
                    m = 47
                    x = np.r_[0.0, np.cumsum(rng.uniform(0.5, 1.5, m-1))]
                    x /= x[-1]
                    y = np.column_stack([
                        np.sin((12+j*5)*x)+0.1*rng.normal(size=m)
                        for j in range(channels)
                    ])
                    if periodic:
                        y[-1] = y[0]
                    w = rng.uniform(0.5, 1.5, m)
                    nest = None
                    bounds = (0.0, 1.0)
                    if mode == 0:
                        target = 0.0
                        w = None
                    elif mode == 1:
                        target = 1e6
                    elif mode == 2:
                        target = 0.8*channels
                    elif mode == 3:
                        target = 0.01*channels
                    elif mode == 4:
                        target = 1e-20
                    elif mode == 5:
                        target = 1e-15
                        nest = 2*(k+1)+5
                    else:
                        target = 0.5*channels
                        if not periodic:
                            bounds = (-0.2, 1.3)
                    actual = native(exe, x, y, k, target, w, periodic, nest, bounds)
                    weights = np.ones_like(x) if w is None else w
                    if target == 0:
                        name = "_periodic_knots" if periodic and k > 1 else "_not_a_knot"
                        expected = [upstream[name](x, k)]
                    elif not periodic:
                        ordinates = y[:, 0] if channels == 1 else y
                        expected = list(generate_knots(x, ordinates, k=k, w=w,
                                                      s=target, nest=nest,
                                                      xb=bounds[0], xe=bounds[1]))
                        scalar_cases += 1
                    else:
                        expected = list(upstream["_generate_knots_impl"](
                            x, y, weights, *bounds, k, target, nest, periodic
                        ))
                        # v1.18.0 yields the last capacity-limited vector before
                        # refreshing cyclic exterior knots. SciFort guarantees
                        # canonical cycles on every output; compare its core
                        # exactly and reconstruct only those stale exteriors.
                        if nest is not None:
                            final = expected[-1].copy()
                            independent = len(final)-2*k-1
                            for j in range(1, k+1):
                                final[k-j] = final[k-j+independent]-1.0
                                final[-k+j-1] = final[k+j]+1.0
                            if not np.array_equal(final, expected[-1]):
                                corrected_exteriors += 1
                            expected[-1] = final
                    label = (k, periodic, channels, mode)
                    assert len(actual) == len(expected), (label, [len(t) for t in actual],
                                                          [len(t) for t in expected])
                    for a, b in zip(actual, expected):
                        assert a.shape == b.shape, (label, a.shape, b.shape)
                        np.testing.assert_allclose(a, b, rtol=2e-14, atol=2e-14, err_msg=str(label))
                        max_error = max(max_error, float(np.max(np.abs(a-b))))
                        comparisons += len(a)
                    if target > 0:
                        last_fp = np.inf
                        for knots in actual:
                            _, fp = residuals(x, y, knots, k, weights, periodic)
                            assert fp <= last_fp+1e-10, ("residual increase", label, fp, last_fp)
                            last_fp = fp
                            residual_checks += 1
                        if nest is None:
                            assert last_fp < target*1.001+1e-20, ("target", label, last_fp)
                    # Balanced data/weight scaling preserves the objective and knots.
                    balanced = native(exe, x, y*1e100, k, target,
                                      weights*1e-100 if target else None, periodic, nest, bounds)
                    assert len(balanced) == len(actual), ("scale count", label)
                    for a, b in zip(actual, balanced):
                        np.testing.assert_array_equal(a, b, err_msg="balanced scaling " + str(label))
                    if mode == 2:
                        for scale in (1e-100, 1e100):
                            scaled = native(exe, x*scale, y, k, target, w, periodic, nest,
                                            tuple(v*scale for v in bounds))
                            assert len(scaled) == len(actual), ("site scale count", label)
                            for a, b in zip(actual, scaled):
                                np.testing.assert_allclose(a, b/scale, rtol=2e-14, atol=2e-14)
                            site_scale_cases += 1
                    cases += 1
                    periodic_cases += int(periodic)
    print(f"SciPy {scipy.__version__}; {cases} sequence cases + {cases} scaling cases; "
          f"{scalar_cases} direct public SciPy comparisons; {periodic_cases} periodic cases; "
          f"{comparisons} knot comparisons; max absolute error {max_error:.6g}; "
          f"{corrected_exteriors} stale upstream capacity exteriors canonicalized; "
          f"{residual_checks} residual checks; {site_scale_cases} site scale cases")


if __name__ == "__main__":
    main()
