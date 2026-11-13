# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Independent solve_bvp API comparisons; requires NumPy and SciPy.

Pass an already compiled validation/bvp_reference_driver.f90 executable.
"""
import argparse
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.integrate import solve_bvp


def reference(family, q, tol, branch):
    x = np.linspace(0, 1, 5)
    y = np.zeros((2, 5), dtype=complex if family >= 5 else float)
    options = dict(tol=tol, max_nodes=3000)
    if family == 0:
        fun = lambda t, v: np.array([q*v[0], -q*v[1]])
        bc = lambda a, b: np.array([a[0]+b[1]-1-2*np.exp(-q),
                                    b[0]+2*a[1]-np.exp(q)-4])
    elif family == 1:
        if branch:
            y[0] = 3
        fun = lambda t, v: np.array([v[1], -q*np.exp(v[0])])
        bc = lambda a, b: np.array([a[0], b[0]])
    elif family == 2:
        fun = lambda t, v: np.array([v[1], -q*v[1]])
        bc = lambda a, b: np.array([a[0], b[0]-1+np.exp(-q)])
    elif family == 3:
        y[0] = 1
        fun = lambda t, v: np.array([v[1], -q*q*v[0]])
        bc = lambda a, b: np.array([a[1], b[0]-np.sin(q)/q])
        options['S'] = np.diag([0, -2])
    elif family == 4:
        y[0] = np.exp(.8*q*x)
        y[1] = 2*np.exp(-1.5*q*x)
        fun = lambda t, v, p: np.array([p[0]*v[0], p[1]*v[1]])
        bc = lambda a, b, p: np.array([a[0]-1, a[1]-2, b[0]-np.exp(q),
                                       b[1]-2*np.exp(-2*q)])
        options['p'] = [.8*q, -1.5*q]
    elif family == 5:
        y[:] = 1
        fun = lambda t, v: np.array([1j*q*v[0], -q*v[1]])
        bc = lambda a, b: np.array([a[0]+b[1]-1-2*np.exp(-q),
                                    b[0]+2*a[1]-np.exp(1j*q)-4])
    else:
        y[0] = 1
        cq = q*(1+.3j)
        fun = lambda t, v: np.array([v[1], -cq*cq*v[0]])
        bc = lambda a, b: np.array([a[1], b[0]-np.sin(cq)/cq])
        options['S'] = np.diag([0, -2]).astype(complex)
        def jac(t, v):
            a = np.zeros((2, 2, len(t)), dtype=complex)
            a[0, 1] = 1
            a[1, 0] = -cq*cq
            return a
        options['fun_jac'] = jac
        options['bc_jac'] = lambda a, b: (np.array([[0, 1], [0, 0]], dtype=complex),
                                          np.array([[0, 0], [1, 0]], dtype=complex))
    return solve_bvp(fun, bc, x, y, **options)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    cases = [(0, q, 0) for q in (-5, -1, 0, 1, 5)]
    cases += [(1, q, branch) for q in (1, 2, 3) for branch in (0, 1)]
    cases += [(2, q, 0) for q in (10, 100, 1000)]
    cases += [(f, q, 0) for f in (3, 4, 5, 6) for q in (.2, 1, 3)]
    checks = 0
    reference_failures = []
    worst = np.zeros(3)
    for tol in (1e-4, 1e-7):
        for family, q, branch in cases:
            result = subprocess.run([str(args.executable.resolve()), str(family), str(q),
                                     str(tol), str(branch)], check=True, capture_output=True,
                                    text=True)
            lines = result.stdout.splitlines()
            status, nodes, iterations, rms, bc = map(float, lines[0].split())
            assert status == 0, (family, q, tol, result.stdout[:200])
            assert rms <= tol and bc <= tol, (family, q, tol, rms, bc)
            parameters = np.fromstring(lines[1], sep=' ')[1:]
            values = np.loadtxt(lines[2:]).reshape(3, 81, 5)
            actual = values[:, :, 1:3].transpose(0, 2, 1)
            if family >= 5:
                actual = actual + 1j*values[:, :, 3:5].transpose(0, 2, 1)
            ref = reference(family, q, tol, branch)
            if not ref.success:
                assert family == 6, (family, q, tol, ref.message)
                reference_failures.append((family, q, tol, ref.message))
            if len(parameters):
                np.testing.assert_allclose(parameters, ref.p, atol=10*tol, rtol=10*tol)
                checks += len(parameters)
            for derivative in range(3):
                expected = ref.sol(values[derivative, :, 0], derivative)
                if not ref.success:
                    # Entire regular Lane-Emden solution, independent of SciPy convergence.
                    t = values[derivative, :, 0]
                    cq = q*(1+.3j)
                    u = np.ones_like(t, dtype=complex)
                    up = np.zeros_like(u)
                    u[1:] = np.sin(cq*t[1:])/(cq*t[1:])
                    up[1:] = (cq*t[1:]*np.cos(cq*t[1:])-np.sin(cq*t[1:]))/(cq*t[1:]**2)
                    upp = np.full_like(u, -cq*cq/3)
                    upp[1:] = -2*up[1:]/t[1:]-cq*cq*u[1:]
                    uppp = np.zeros_like(u)
                    uppp[1:] = -2*upp[1:]/t[1:]+2*up[1:]/t[1:]**2-cq*cq*up[1:]
                    expected = np.array(((u, up), (up, upp), (upp, uppp))[derivative])
                # Piecewise-cubic second derivatives differ with the chosen meshes.
                limit = (tol, tol**.75, tol**.5)[derivative]
                error = np.abs(actual[derivative]-expected)/(1+np.abs(expected))
                worst[derivative] = max(worst[derivative], np.max(error)/tol**(1, .75, .5)[derivative])
                assert np.max(error) <= limit, (family, q, tol, derivative, np.max(error))
                checks += error.size
            # Analytic oracles independent of either collocation implementation.
            t = values[0, :, 0]
            if family == 0:
                exact = np.array([np.exp(q*t), 2*np.exp(-q*t)])
            elif family == 2:
                exact = np.array([1-np.exp(-q*t), q*np.exp(-q*t)])
            elif family == 3:
                exact = np.array([np.sinc(q*t/np.pi), np.zeros_like(t)])
                exact[1, 1:] = (q*t[1:]*np.cos(q*t[1:])-np.sin(q*t[1:]))/(q*t[1:]**2)
            elif family == 4:
                exact = np.array([np.exp(q*t), 2*np.exp(-2*q*t)])
            elif family == 5:
                exact = np.array([np.exp(1j*q*t), 2*np.exp(-q*t)])
            else:
                exact = None
            if exact is not None:
                oracle_tol = 5*tol
                if family == 0:
                    # Coupled endpoint constraints amplify the local defect at large |q|.
                    boundary = np.array([[1, np.exp(-q)], [np.exp(q), 2]])
                    oracle_tol *= np.sqrt(np.linalg.cond(boundary))
                np.testing.assert_allclose(actual[0], exact, atol=oracle_tol, rtol=oracle_tol,
                                           err_msg=str((family, q, tol)))
                checks += exact.size
    print(f'SciPy {scipy.__version__}: {2*len(cases)} BVP cases, {checks} checked values passed.')
    print('Maximum error / tol**(1, .75, .5):', worst)
    for failure in reference_failures:
        print('SciPy reference did not converge; analytic oracle used:', failure)


if __name__ == '__main__':
    main()
