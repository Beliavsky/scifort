#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare scattered spherical smoothing with SciPy and independent polynomial LSQ."""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.interpolate import BSpline, SmoothSphereBivariateSpline


def native_batch(exe, cases):
    lines = []
    for case in cases:
        t, p, r, w, s, eps, nt, np_, qt, qp = case
        lines += [f"{len(t)} {len(qt)} {nt} {np_}", f"{s:.17g} {eps:.17g}"]
        for row in (t, p, r, w):
            lines.append(" ".join(f"{x:.17g}" for x in row))
        lines += [f"{x:.17g} {y:.17g}" for x, y in zip(qt, qp)]
    result = subprocess.run([str(exe.resolve())], input="\n".join(lines)+"\n",
                            text=True, capture_output=True, check=True)
    rows = iter(result.stdout.splitlines())
    results = []
    for case in cases:
        meta = np.fromstring(next(rows), sep=" ")
        if int(meta[0]) not in (0, 2):
            results.append((meta, None, None, None))
            continue
        knots = [np.fromstring(next(rows), sep=" ") for _ in range(2)]
        coef = np.fromstring(next(rows), sep=" ")
        values = np.array([np.fromstring(next(rows), sep=" ") for _ in case[-2]])
        results.append((meta, knots, coef, values))
    assert next(rows, None) is None
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261116)
    cases, refs = [], []
    for m in (60, 120, 240):
        t = rng.uniform(0.05, np.pi-0.05, m)
        p = rng.uniform(0, 2*np.pi, m)
        t[:4] = [0, np.pi, np.pi/2, np.pi/2]
        p[:4] = [0, 2*np.pi, 0, 2*np.pi]
        qt = np.r_[t, 0, np.pi, rng.uniform(0, np.pi, 25)]
        qp = np.r_[p, 0.3, 0.7, rng.uniform(0, 2*np.pi, 25)]
        w = rng.uniform(0.7, 1.8, m)
        z = t/np.pi
        for kind in ("cubic", "harmonic", "noisy"):
            r = 1+2*z*z*(3-2*z)
            if kind != "cubic":
                r = np.cos(t)+0.7*np.sin(t)*np.cos(p)+0.2*np.sin(t)**2*np.sin(2*p)
            if kind == "noisy":
                r += rng.normal(0, 0.04, m)
            for s in (0.002*m, 0.02*m, 10.0*m):
                ref = SmoothSphereBivariateSpline(t, p, r, w=w, s=s)
                cap = 8+int(np.sqrt(m//2))
                cases.append((t, p, r, w, s, 1e-16, cap, cap, qt, qp))
                refs.append(ref)
    outputs = native_batch(args.executable, cases)
    max_error = max_coef = max_residual = 0.0
    comparisons = polynomials = 0
    for case, ref, (meta, knots, coef, values) in zip(cases, refs, outputs):
        t, p, r, w, s, _, _, _, qt, qp = case
        assert meta[0] in (0, 2), (len(t), s, meta)
        for a, b in zip(knots, ref.get_knots()):
            assert a.shape == b.shape and np.allclose(a, b, rtol=0, atol=3e-13)
        ce = np.max(abs(coef-ref.get_coeffs()))/(1+np.max(abs(ref.get_coeffs())))
        expected = np.column_stack([ref.ev(qt, qp, dtheta=a, dphi=b)
            for a, b in ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1))])
        error = np.max(abs(values-expected)/(1+abs(expected)))
        re = abs(meta[3]-ref.get_residual())/(1+ref.get_residual())
        assert ce < 3e-8 and error < 3e-7 and re < 3e-9, (len(t), s, ce, error, re)
        fp = np.sum((w*(r-values[:len(t), 0]))**2)
        assert abs(fp-meta[3]) < 3e-10*(1+fp), (fp, meta)
        assert fp <= 1.0011*s
        if all(len(a) == 8 for a in knots):
            x = t/np.pi
            h = x*x*(3-2*x)
            basis = np.column_stack((1-h, h))
            beta = np.linalg.lstsq(w[:, None]*basis, w*r, rcond=None)[0]
            qx = qt/np.pi
            qh = qx*qx*(3-2*qx)
            independent = beta[0]*(1-qh)+beta[1]*qh
            assert np.max(abs(values[:, 0]-independent)) < 3e-11
            assert meta[2] == 2 and meta[1] == -1
            polynomials += 1
        max_error = max(max_error, error)
        max_coef = max(max_coef, ce)
        max_residual = max(max_residual, re)
        comparisons += values.size
    print(f"PASS: {len(cases)} scattered fits, {comparisons} value/derivative references, "
          f"{polynomials} independent polynomial solves; SciPy {scipy.__version__}")
    print(f"Maximum normalized errors: value {max_error:.6g}, coefficients {max_coef:.6g}, "
          f"residual {max_residual:.6g}")
    check_rank_deficiency(args.executable)


def check_rank_deficiency(executable):
    """Build the constrained tensor space explicitly, independently of Givens QR."""
    tt = np.r_[[0.0]*4, np.pi/2, [np.pi]*4]
    nodes = np.arange(4)*np.pi/2
    tp = np.arange(-3, 8)*np.pi/2
    longitude = BSpline(tp, np.eye(7), 3)(nodes)
    fold = np.zeros((7, 4))
    fold[np.arange(7), np.arange(7) % 4] = 1
    cs = np.linalg.solve(longitude @ fold,
                         np.column_stack((np.cos(nodes), np.sin(nodes))))
    cs = fold @ cs
    mapping = np.zeros((35, 10))
    for j in range(7):
        mapping[j, 0] = 1
        mapping[7+j, :3] = [1, *cs[j]]
        mapping[14+j, 3+j % 4] = 1
        mapping[21+j, 7:] = [*cs[j], 1]
        mapping[28+j, 9] = 1

    def design(t, p, dt=0, dp=0):
        bt = BSpline(tt, np.eye(5), 3)(t, nu=dt)
        bp = BSpline(tp, np.eye(7), 3)(p, nu=dp)
        return (bt[:, :, None]*bp[:, None, :]).reshape(-1, 35) @ mapping

    p = 2*np.pi*np.arange(40)/40
    t = np.full(40, np.pi/2)
    r = 2+np.sin(p)
    a = design(t, p)
    beta, _, rank, _ = np.linalg.lstsq(a, r, rcond=1e-8)
    assert rank == 4
    residual = np.sum((r-a @ beta)**2)
    qt = np.r_[t, 0.3, np.pi/2, 2.6]
    qp = np.r_[p, 0.2, 0.7, 4.1]
    case = (t, p, r, np.ones(40), residual*1.00001, 1e-16, 12, 12, qt, qp)
    meta, knots, coefficients, values = native_batch(executable, [case])[0]
    assert tuple(meta[:3]) == (2, 4, 10), meta
    assert abs(meta[3]-residual) < 3e-12
    expected = np.column_stack([design(qt, qp, dt, dp) @ beta
        for dt, dp in ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1))])
    error = np.max(abs(values-expected))
    assert error < 3e-10, error
    assert np.max(abs(coefficients-mapping @ beta)) < 3e-10
    assert np.allclose(knots[0], tt, rtol=0, atol=3e-15)
    assert np.allclose(knots[1], tp, rtol=0, atol=3e-15)
    print(f"PASS: rank-4 minimum-norm fit, {values.size} dense-SVD references; "
          f"maximum absolute error {error:.6g}; residual {residual:.16g}")


if __name__ == "__main__":
    main()
