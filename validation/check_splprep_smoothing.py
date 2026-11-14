#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare native shared-knot parametric smoothing with independent implementations."""
import argparse
from pathlib import Path
import re
import subprocess

import numpy as np
import scipy
from scipy.interpolate import BSpline, splprep, splev


def encode(cases):
    lines = []
    for k, x, w, u, bounds, target, q in cases:
        lines.append(f"{len(x)} {x.shape[1]} {k} {int(u is not None)} {len(q)}")
        arrays = [[target, *bounds], *x, w]
        if u is not None:
            arrays.append(u)
        arrays.append(q)
        lines.extend(" ".join(f"{v:.17g}" for v in row) for row in arrays)
    return "\n".join(lines)+"\n"


def run(executable, cases):
    result = subprocess.run([str(executable.resolve())], input=encode(cases),
                            text=True, capture_output=True, check=True)
    rows = iter(result.stdout.splitlines())
    output = []
    for index, (k, x, _, _, _, _, _) in enumerate(cases):
        status = int(next(rows))
        assert status == 0, (index, "fit status", status)
        nt = int(next(rows))
        u = np.fromstring(next(rows), sep=" ")
        t = np.fromstring(next(rows), sep=" ")
        c = np.array([np.fromstring(next(rows), sep=" ") for _ in range(x.shape[1])])
        samples = np.column_stack([np.fromstring(next(rows), sep=" ")
                                   for _ in range(x.shape[1])])
        queries = [np.column_stack([np.fromstring(next(rows), sep=" ")
                                    for _ in range(x.shape[1])])
                   for _ in range(min(3, k)+1)]
        assert len(t) == nt and c.shape == (x.shape[1], nt-k-1)
        output.append((u, t, c, samples, queries))
    assert next(rows, None) is None
    return output


def scipy_results(cases):
    results = []
    for index, (k, x, w, u, bounds, target, q) in enumerate(cases):
        (tck, ru), _, status, _ = splprep(x.T, u=u, w=w, k=k, s=target,
                                         ub=bounds[0], ue=bounds[1], full_output=True)
        assert status <= 0, (index, "SciPy status", status)
        results.append((ru, tck[0], np.asarray(tck[1]), np.array(splev(ru, tck)).T,
                        [np.array(splev(q, tck, der=d)).T for d in range(min(3, k)+1)]))
    return results


def stationarity(k, x, w, bounds, result):
    """Independently check the derivative-jump penalized normal equations."""
    u, t, c, _, _ = result
    span = bounds[1]-bounds[0]
    t = (t-bounds[0])/span
    u = (u-bounds[0])/span
    xs = np.max(abs(x)) or 1.0
    w = w/max(w)
    c = c.T/xs
    y = x/xs
    basis = BSpline(t, np.eye(len(c)), k)
    design = basis(u)
    gradient = design.T@(w[:, None]**2*(design@c-y))
    energy = max(1., np.linalg.norm(design.T@(w[:, None]**2*y)))
    interior = np.unique(t[k+1:-k-1])
    if not len(interior):
        error = np.linalg.norm(gradient)/energy
        assert error < 2e-9, ("polynomial normal equations", error)
        return error
    derivative = basis.derivative(k)
    jump = derivative(np.nextafter(interior, np.inf))-derivative(
        np.nextafter(interior, -np.inf))
    jump /= np.max(abs(jump))
    penalty = jump.T@(jump@c)
    multiplier = max(0., -np.sum(gradient*penalty)/np.sum(penalty**2))
    denominator = max(np.linalg.norm(gradient), multiplier*np.linalg.norm(penalty),
                      energy)
    error = np.linalg.norm(gradient+multiplier*penalty)/denominator
    assert error < 2e-9, ("penalized normal equations", error,
                          np.linalg.norm(gradient), np.linalg.norm(penalty),
                          energy, multiplier)
    return error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--upstream-executable", type=Path)
    parser.add_argument("--upstream-tests", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261113)
    cases = []
    for k in range(1, 6):
        for i in range(30):
            m = k+1 if i == 0 else 25+3*i
            ndim = 1+i % 5
            t = np.linspace(0, 1, m)
            x = np.column_stack([np.sin((3+j)*t)+0.13*rng.normal(size=m)+j*0.1
                                 for j in range(ndim)])
            w = np.exp(rng.uniform(-0.6, 0.6, m))
            u = None
            bounds = (0., 1.)
            if i % 2 == 0:
                u = np.cumsum(rng.uniform(0.6, 1.4, m))
                u = (u-u[0])/(u[-1]-u[0])*3-1
                bounds = (u[0]-(0.2 if i % 3 == 0 else 0),
                          u[-1]+(0.3 if i % 3 == 0 else 0))
            target = np.sum((x*w[:, None])**2)*(0.0001, 0.01, 0.1, 2, 0.005)[i % 5]
            cases.append((k, x, w, u, bounds, target, np.linspace(*bounds, 31)))
    if args.upstream_tests:
        source = args.upstream_tests.read_text(encoding="utf-8")
        fixture = source.split("function test_parametric_fit()", 1)[1]
        fixture = fixture.split("end function test_parametric_fit", 1)[0]
        u_text = re.search(r"\bu\s*=\s*\[real\(FP_REAL\)\s*::([^]]+)\]", fixture).group(1)
        u = np.array([float(v) for v in u_text.replace("&", "").split(",")])
        channels = []
        for j in (1, 2):
            literals = re.search(rf"x\({j},:\)\s*=\s*\[([^]]+)\]", fixture).group(1)
            channels.append(np.array([float(v) for v in literals.replace("&", "").split(",")],
                                     dtype=np.float32).astype(float))
        x = np.column_stack(channels)
        for target in (100., 1., 0.05, 0.25):
            cases.append((3, x, np.ones(len(x)), u, (u[0], u[-1]), target,
                          np.linspace(u[0], u[-1], 31)))
        for k in (3, 5):
            cases.append((k, x, np.ones(len(x)), None, (0., 1.), 0.25,
                          np.linspace(0, 1, 31)))
    actual = run(args.executable, cases)
    normalized = []
    scales = []
    for k, x, w, u, bounds, target, q in cases:
        xs = np.max(abs(x)) or 1.0
        ws = np.max(w)
        span = bounds[1]-bounds[0]
        norm_u = (u-bounds[0])/span if u is not None else None
        norm_q = (q-bounds[0])/span
        norm_s = np.exp(np.log(target)-2*np.log(xs)-2*np.log(ws))
        normalized.append((k, x/xs, w/ws, norm_u, (0., 1.), norm_s, norm_q))
        scales.append((xs, bounds[1]-bounds[0], bounds))
    if args.upstream_executable:
        references = run(args.upstream_executable, normalized)
        raw_references = run(args.upstream_executable, cases)
    else:
        references = scipy_results(normalized)
        raw_references = scipy_results(cases)
    maximum = residual_max = 0.0
    raw_maximum = prediction_budget_max = sample_distance_max = 0.0
    raw_worst = None
    stationarity_max = representation_max = 0.0
    count = adaptive = polynomial = 0
    for index, ((k, x, w, u, bounds, target, q), result) in enumerate(zip(cases, actual)):
        au, at, ac, av, aq = result
        try:
            stat = stationarity(k, x, w, bounds, result)
        except AssertionError as error:
            raise AssertionError((index, k, error.args)) from error
        stationarity_max = max(stationarity_max, stat)
        independent = BSpline(at, ac.T, k)
        for a, b in [(av, independent(au)),
                     *[(aq[d], independent(q, nu=d)) for d in range(min(3, k)+1)]]:
            error = np.max(abs(a-b)/(1+abs(b)))
            assert error < 3e-8, (index, "independent representation evaluation", error)
            representation_max = max(representation_max, error)
        ru, rt, rc, rv, rq = references[index]
        xs, span, bounds = scales[index]
        ru = (1-ru)*bounds[0]+ru*bounds[1]
        rt = (1-rt)*bounds[0]+rt*bounds[1]
        rc = rc*xs
        rv = rv*xs
        rq = [values*(xs/span**d) for d, values in enumerate(rq)]
        np.testing.assert_allclose(au, ru, rtol=3e-13, atol=3e-13)
        assert len(at) == len(rt), (index, "knot counts", len(at), len(rt))
        np.testing.assert_allclose(at, rt, rtol=3e-13, atol=3e-13)
        fp = np.sum((w[:, None]*(av-x))**2)
        assert fp <= target*1.0011+1e-24, (index, "residual", fp, target)
        if len(at) > 2*(k+1):
            adaptive += 1
            residual_max = max(residual_max, abs(fp-target)/target)
        else:
            polynomial += 1
        for a, b in [(ac, np.asarray(rc)), (av, rv), *zip(aq, rq)]:
            error = np.max(abs(a-b)/(1+abs(b)))
            if args.upstream_executable:
                assert error < 3e-8, (index, k, "normalized upstream error", error)
            maximum = max(maximum, error)
            count += a.size
        raw_u, raw_t, raw_c, raw_v, raw_q = raw_references[index]
        np.testing.assert_allclose(au, raw_u, rtol=3e-13, atol=3e-13)
        np.testing.assert_allclose(at, raw_t, rtol=3e-13, atol=3e-13)
        raw_fp = np.sum((w[:, None]*(raw_v-x))**2)
        assert raw_fp <= target*1.0011+1e-24, (index, "raw reference residual")
        # With fixed knots, residual components in the penalty eigenbasis are
        # monotone in the smoothing parameter. Thus squared prediction distance
        # is bounded by the difference in squared residuals, plus roundoff.
        distance = np.sum((w[:, None]*(av-raw_v))**2)
        energy = max(1., np.sum((w[:, None]*x)**2))
        budget = abs(fp-raw_fp)+1e-10*energy
        assert distance <= budget, (index, "smoothing-path prediction budget", distance, budget)
        prediction_budget_max = max(prediction_budget_max, distance/budget)
        sample_distance_max = max(sample_distance_max, np.sqrt(distance/target))
        comparisons = [("coefficients", ac, raw_c), ("samples", av, raw_v)]
        comparisons.extend((f"query derivative {d}", a, b)
                           for d, (a, b) in enumerate(zip(aq, raw_q)))
        for label, a, b in comparisons:
            error = np.max(abs(a-b)/(1+abs(b)))
            if error > raw_maximum:
                raw_maximum = error
                raw_worst = (index, k, label, float(np.max(abs(a-b))))
    assert adaptive > 75 and polynomial > 20
    name = "original perazz/FITPACK" if args.upstream_executable else f"SciPy {scipy.__version__}"
    print(f"PASS: {len(cases)} fits, {adaptive} adaptive, {polynomial} polynomial; "
          f"{count} coefficient/value/derivative diagnostics with {name}")
    if args.upstream_executable:
        print("All normalized upstream coefficient/value/derivative errors are below 3e-8")
    print(f"Normalized-reference maximum difference {maximum:.6g}; "
          f"maximum adaptive residual relative error {residual_max:.6g}")
    print(f"Raw-reference maximum normalized difference {raw_maximum:.6g}; "
          f"maximum squared prediction distance / residual budget {prediction_budget_max:.6g}")
    print(f"Independent stationarity maximum {stationarity_max:.6g}; "
          f"independent representation evaluation maximum {representation_max:.6g}")
    print(f"Maximum raw weighted sample distance / sqrt(target) {sample_distance_max:.6g}; "
          f"worst raw diagnostic (case, degree, quantity, max absolute difference) {raw_worst}")


if __name__ == "__main__":
    main()
