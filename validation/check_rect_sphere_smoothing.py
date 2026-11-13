#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check native adaptive rectangular spherical splines against SciPy."""
import argparse
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import RectSphereBivariateSpline, BSpline

from check_rect_sphere import native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--upstream", type=Path, help="Pinned FITPACK driver with documented fixes")
    args = parser.parse_args()
    rng = np.random.default_rng(20261115)
    configs = [((False, False), (None, None), (False, False), (False, False)),
        ((True, True), (None, None), (False, False), (False, False)),
        ((True, False), (None, None), (False, False), (False, False)),
        ((False, True), (None, None), (False, False), (False, False)),
        ((True, True), (None, None), (False, False), (True, True)),
        ((True, False), (1.2, None), (True, False), (True, False)),
        ((False, True), (None, -0.8), (False, True), (False, True)),
        ((True, True), (1.2, -0.8), (False, False), (False, False))]
    count = comparisons = scipy_count = upstream_count = polynomial_count = 0
    max_value = max_coef = max_residual = 0.0
    for mu, mv, origin in ((6, 8, -1.2), (9, 12, 0.65), (12, 16, 0.0)):
        u = np.linspace(0.12, np.pi-0.17, mu)
        v = origin+np.arange(mv)*2*np.pi/mv
        if mu == 12:
            u[1:-1] += rng.uniform(-0.03, 0.03, mu-2)
            v[1:] += rng.uniform(-0.02, 0.02, mv-1)
        nt, np_ = np.meshgrid(u, v, indexing="ij")
        r = (np.cos(nt)+0.7*np.sin(nt)*np.cos(np_)
             +0.2*np.sin(nt)**2*np.sin(2*np_)+rng.normal(0, 0.04, nt.shape))
        qt = np.r_[nt.ravel(), 0, np.pi, rng.uniform(0, np.pi, 25)]
        qp = np.r_[np_.ravel(), origin, origin+2*np.pi,
                   rng.uniform(origin, origin+2*np.pi, 25)]
        for continuity, poles, exact, flat in configs:
            for s in (0.01, 0.2, 2.0, 1000.0):
                knots, coef, values, meta = native(args.executable, u, v, r,
                    continuity, poles, exact, flat, qt, qp, s=s, metadata=True)
                label = (mu, mv, continuity, poles, exact, flat, s)
                references = []
                # Legacy SciPy 1.15.3 C1 pole fits can be nonfinite or depend on
                # prior calls. Restrict strict SciPy comparison to C0 poles.
                if not any(continuity):
                    ref = RectSphereBivariateSpline(u, v, r, s=s)
                    expected = np.column_stack([ref.ev(qt, qp, dtheta=a, dphi=b)
                        for a, b in ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1))])
                    references.append((ref.get_knots(), ref.get_coeffs(), expected,
                                       ref.get_residual()))
                    scipy_count += 1
                if args.upstream:
                    rk, rc, ry, rm = native(args.upstream, u, v, r, continuity,
                        poles, exact, flat, qt, qp, s=s, metadata=True)
                    references.append((rk, rc, ry, rm[3]))
                    upstream_count += 1
                for rk, rc, expected, residual in references:
                    for a, b in zip(knots, rk):
                        assert a.shape == b.shape and np.allclose(a, b, rtol=0, atol=3e-15), label
                    ce = np.max(abs(coef-rc))/(1+np.max(abs(rc)))
                    error = np.max(abs(values-expected)/(1+abs(expected)))
                    re = abs(meta[3]-residual)/(1+residual)
                    assert ce < 3e-8 and error < 3e-7 and re < 3e-9, (label, ce, error, re)
                    max_value = max(max_value, error)
                    max_coef = max(max_coef, ce)
                    max_residual = max(max_residual, re)
                    comparisons += values.size
                if s == 1000.0:
                    # Independent constrained cubic LSQ oracle for the no-knot
                    # limit. C1 periodic longitude-constant fits have flat poles.
                    groups = [0, 0 if continuity[0] else 1,
                              3 if continuity[1] else 2, 3]
                    fixed = np.zeros(4)
                    fixed_groups = {}
                    for pole, value, is_exact in zip((0, 3), poles, exact):
                        if value is not None and is_exact:
                            fixed_groups[groups[pole]] = value
                    free = sorted(set(groups)-fixed_groups.keys())
                    mapping = np.zeros((4, len(free)))
                    for i, group in enumerate(groups):
                        if group in fixed_groups:
                            fixed[i] = fixed_groups[group]
                        else:
                            mapping[i, free.index(group)] = 1
                    tk = np.r_[[0.0]*4, [np.pi]*4]
                    basis = BSpline.design_matrix(u, tk, 3).toarray()
                    design = np.repeat(basis @ mapping, mv, axis=0)
                    rhs = r.ravel()-np.repeat(basis @ fixed, mv)
                    for pole, value, is_exact in zip((0, 3), poles, exact):
                        if value is not None and not is_exact:
                            design = np.vstack((design, mapping[pole]))
                            rhs = np.r_[rhs, value-fixed[pole]]
                    coefficients = fixed+mapping @ np.linalg.lstsq(design, rhs, rcond=None)[0]
                    poly = BSpline(tk, coefficients, 3)
                    for column, derivative in ((0, 0), (1, 1), (3, 2)):
                        assert np.max(abs(values[:, column]-poly(qt, nu=derivative))) < 3e-10, label
                    assert np.max(abs(values[:, [2, 4, 5]])) < 3e-11, label
                    polynomial_count += 1
                # FITPACK's public residual excludes soft pole-observation residuals.
                grid_residual = np.sum((values[:mu*mv, 0]-r.ravel())**2)
                assert abs(grid_residual-meta[3]) < 3e-11*(1+meta[3]), label
                assert meta[3] <= s*1.0011, label
                assert meta[1] == -1, label
                count += 1
    print(f"PASS: {count} native fits, {scipy_count} SciPy and {upstream_count} upstream fits, "
          f"{polynomial_count} independent polynomial fits; {comparisons} references; "
          f"SciPy {scipy.__version__}")
    print(f"Maximum normalized errors: value {max_value:.6g}, "
          f"coefficients {max_coef:.6g}, residual {max_residual:.6g}")


if __name__ == "__main__":
    main()
