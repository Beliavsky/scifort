#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare native fixed-knot spherical splines with SciPy and analytic data."""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.interpolate import LSQSphereBivariateSpline


def native(exe, theta, phi, data, weights, tt, tp, qt, qp, eps=1e-16):
    payload = [f"{len(theta)} {len(tt)} {len(tp)} {len(qt)} {eps:.17g}"]
    payload += [" ".join(f"{v:.17g}" for v in row)
                for row in zip(theta, phi, data, weights)]
    if len(tt):
        payload.append(" ".join(f"{v:.17g}" for v in tt))
    payload.append(" ".join(f"{v:.17g}" for v in tp))
    payload += [f"{t:.17g} {p:.17g}" for t, p in zip(qt, qp)]
    result = subprocess.run([str(exe.resolve())], input="\n".join(payload)+"\n",
                            capture_output=True, text=True, check=True)
    lines = result.stdout.splitlines()
    status, rank, dimension, residual = map(float, lines[0].split())
    assert status in (0, 2), (status, result.stdout, result.stderr)
    knots = [np.fromstring(line, sep=" ") for line in lines[1:3]]
    coef = np.fromstring(lines[3], sep=" ")
    values = np.array([np.fromstring(line, sep=" ") for line in lines[4:]])
    assert values.shape == (len(qt), 6)
    return int(status), int(rank), int(dimension), residual, knots, coef, values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(20261113)
    qt = np.r_[0, np.pi, np.linspace(0.02, np.pi-0.02, 35)]
    qp = np.r_[0, 2*np.pi, rng.uniform(0, 2*np.pi, 35)]
    max_eval_error = max_coef_error = max_residual_error = 0.0
    cases = 0
    for nlat, nlon in ((0, 3), (1, 3), (3, 5), (4, 7), (2, 4)):
        tt = np.linspace(0, np.pi, nlat+2)[1:-1]
        tp = np.linspace(0, 2*np.pi, nlon+2)[1:-1]
        if (nlat,nlon) == (2,4):
            tt = np.array([0.35, 2.3])
            tp = np.array([0.3, 1.1, 3.4, 5.6])
        for kind in ("constant", "axis", "harmonic", "noisy"):
            theta = rng.uniform(0, np.pi, 240)
            phi = rng.uniform(0, 2*np.pi, 240)
            # Repeated pole observations, both seam representations, unsorted data.
            theta[:8] = [0, 0, np.pi, np.pi, 0.4, 0.4, 2.5, 2.5]
            phi[:8] = [0, 2, 1, 2*np.pi, 0, 2*np.pi, 0, 2*np.pi]
            data = np.ones_like(theta)*2.5
            if kind != "constant":
                # A global cubic with zero polar slope is exactly in every space.
                u = theta/np.pi
                data = 1+2*(3*u*u-2*u*u*u)
            if kind in ("harmonic", "noisy"):
                data += np.sin(theta)*np.cos(phi)+0.3*np.sin(theta)**2*np.sin(2*phi)
            if kind == "noisy":
                data += rng.normal(0, 0.02, len(data))
            weights = rng.uniform(0.3, 3, len(data))
            got = native(args.executable,theta,phi,data,weights,tt,tp,qt,qp)
            status, rank, dimension, residual, knots, coef, values = got
            assert status == 0 and rank == dimension == 6+(nlon+1)*nlat
            ref = LSQSphereBivariateSpline(theta,phi,data,tt,tp,w=weights)
            for a, b in zip(knots, ref.get_knots()):
                assert np.allclose(a, b, rtol=0, atol=2e-15)
            coef_error = np.max(np.abs(coef-ref.get_coeffs()))/(1+np.max(np.abs(coef)))
            assert coef_error < 2e-10, (nlat,nlon,kind,coef_error)
            max_coef_error = max(max_coef_error,coef_error)
            expected = np.column_stack([ref.ev(qt,qp,dtheta=a,dphi=b)
                for a,b in ((0,0),(1,0),(0,1),(2,0),(0,2),(1,1))])
            error = np.max(np.abs(values-expected)/(1+np.abs(expected)))
            assert error < 3e-9, (nlat,nlon,kind,error)
            max_eval_error = max(max_eval_error,error)
            residual_error = abs(residual-ref.get_residual())/(1+ref.get_residual())
            assert residual_error < 2e-10, (nlat,nlon,kind,residual_error)
            max_residual_error = max(max_residual_error,residual_error)
            if kind in ("constant", "axis"):
                u = qt/np.pi
                analytic = np.full_like(qt,2.5) if kind == "constant" else 1+2*(3*u*u-2*u*u*u)
                assert np.max(abs(values[:,0]-analytic)) < 1e-11
            cases += 1
    print(f"PASS: {cases} fits, {cases*len(qt)*6} value/derivative references; SciPy {scipy.__version__}")
    print(f"Maximum normalized errors: evaluation {max_eval_error:.6g}, "
          f"coefficients {max_coef_error:.6g}, residual {max_residual_error:.6g}")


if __name__ == "__main__":
    main()
