# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Compare native AAA fits with SciPy and analytic rational spectra."""
import argparse
from pathlib import Path
import subprocess
import warnings

import numpy as np
import scipy
from scipy.interpolate import AAA
from scipy.special import gamma
from scipy.optimize import linear_sum_assignment


def samples(family, a):
    x = -2+4*np.arange(201)/200
    if family == 0:
        f = (x+.2)/(x*x+.4*x+1+a*a)
    elif family == 1:
        f = np.exp(a*x)
    elif family == 2:
        f = np.tanh(a*x)
    elif family == 3:
        f = np.sin(a*x)
    elif family == 4:
        x = np.exp(-a*(2-x)/4)
        f = np.log(x)
    elif family == 5:
        x = 2*np.exp(1j*np.pi*(x+2)/2.01)
        p, q = a*(.3+.4j), a*(-.2+.6j)
        f = (x+.2)/((x-p)*(x-q))
    elif family == 6:
        x = .75*x+.003
        f = gamma(x)
    elif family == 7:
        f = 1/(x+3)+a/(x-.313)
    else:
        f = 1/(x+3)+1e-8*np.sin(a*x)
    return x.astype(complex), f.astype(complex)


def match(actual, expected, tolerance):
    assert len(actual) == len(expected), (actual, expected)
    if not len(actual):
        return 0
    expected = np.asarray(expected)
    error = abs(actual[:, None]-expected[None, :])/(1+abs(expected[None, :]))
    i, j = linear_sum_assignment(error)
    assert np.max(error[i, j]) <= tolerance, (actual, expected, np.max(error[i, j]))
    return len(actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    cases = [(family, a, tol, 48, 1)
             for family, coefficients in [(0, (.5, 1., 2.)), (1, (.5, 1., 2.)),
                                           (2, (1., 3.)), (3, (1., 5.)),
                                           (4, (3., 10.)), (5, (.5, 1.)), (6, (1.,))]
             for a in coefficients for tol in (1e-8, 1e-11)]
    cases += [(7, 1e-10, 1e-13, 12, c) for c in (0, 1, 2)]
    cases += [(8, 127., 1e-12, 12, 0)]
    checks = spectral = 0
    diagnostics = []
    worst = 0.
    for family, a, tol, budget, clean in cases:
        output = subprocess.run([str(args.executable.resolve()), str(family), str(a),
                                 str(tol), str(budget), str(clean)], capture_output=True,
                                text=True, check=True).stdout.splitlines()
        status, fit, cleanup, terms, sample_error = map(float, output[0].split())
        assert terms > 0, (family, a, tol, status)
        assert fit in (0, 2), (family, a, tol, status, fit)
        data = np.loadtxt(output[1:122])
        assert np.all(data[:, 4] == 0), (family, a, tol, data[data[:, 4] != 0])
        query, actual = data[:, 0]+1j*data[:, 1], data[:, 2]+1j*data[:, 3]
        x, f = samples(family, a)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            ref = AAA(x, f, rtol=tol, max_terms=budget, clean_up=bool(clean),
                      clean_up_tol=1e-6 if clean == 2 else 1e-13)
        expected = ref(query)
        # AAA's tolerance controls error on the samples, scaled by the largest value.
        # Different SVDs can choose different supports and rational representations.
        bound = 8*max(tol, np.finfo(float).eps**.75)*(np.max(abs(f))+abs(expected))
        if family == 7 and clean == 2:
            bound = 1e-7*(1+abs(expected))
        if family == 8:
            bound = 1e-5*(1+abs(expected))
        error = abs(actual-expected)
        assert np.all(error <= bound), (family, a, tol, np.max(error/bound))
        worst = max(worst, np.max(error/bound))
        checks += len(query)
        ps, rs, ds, npoles, nroots, nres = map(int, output[122].split())
        pos = 123
        arrays = []
        for count in (npoles, nroots, nres):
            values = np.loadtxt(output[pos:pos+count]).reshape(count, 2) if count else np.empty((0, 2))
            arrays.append(values[:, 0]+1j*values[:, 1])
            pos += count
        poles, roots, residues = arrays
        assert pos == len(output)
        if cleanup or ps or rs or ds:
            diagnostics.append((family, a, tol, int(cleanup), ps, rs, ds))
        if family in (0, 5):
            assert (cleanup, ps, rs, ds) == (0, 0, 0, 0), (family, a, tol, diagnostics[-1:])
            if family == 0:
                exact_poles = np.array([-.2+1j*np.sqrt(a*a+.96), -.2-1j*np.sqrt(a*a+.96)])
                exact_residues = .5+0*exact_poles
            else:
                exact_poles = a*np.array([.3+.4j, -.2+.6j])
                exact_residues = (exact_poles+.2)/(exact_poles-exact_poles[::-1])
            spectral += match(poles, exact_poles, 1e-8)
            spectral += match(roots, [-.2], 1e-8)
            rows, cols = linear_sum_assignment(abs(poles[:, None]-exact_poles[None, :]))
            np.testing.assert_allclose(residues[rows], exact_residues[cols], atol=1e-8, rtol=1e-8)
            spectral += len(residues)
        if family == 7:
            assert ps == 0 and ds == 0, (family, clean, ps, ds)
            target = [-3.] if clean == 2 else [-3., .313]
            spectral += match(poles, target, 1e-6)
    print(f'SciPy {scipy.__version__}: {len(cases)} cases, {checks} query values, '
          f'{spectral} analytic spectral values passed.')
    print(f'Maximum error/bound: {worst:.4g}')
    for diagnostic in diagnostics:
        print('Reported cleanup/pole/root/residue diagnostic:', diagnostic)


if __name__ == '__main__':
    main()
