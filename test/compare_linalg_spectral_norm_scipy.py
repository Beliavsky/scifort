#!/usr/bin/env python3
"""Independent spectral-norm comparisons with SciPy and NumPy.

The Fortran program is given an initial vector to make its power iteration
bitwise-reproducible between checked and optimized builds. SciPy draws a
DIFFERENT initial vector, so its randomized estimate need not be identical.
This script independently recomputes the fixed-vector iteration, compares
its result to Fortran, and checks both estimates against the exact SVD norm.

Usage: python test/compare_linalg_spectral_norm_scipy.py ./driver
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import scipy
from scipy.linalg import interpolative


def cases():
    rng = np.random.default_rng(20261109)
    # Two dtypes * two operations * multiple rectangular dimensions.
    for kind in range(4):
        complex_type = kind >= 2
        for m in (0, 1, 2, 3, 6, 11, 21):
            for n in (0, 1, 2, 4, 8, 17):
                for variant in range(5):
                    a = rng.normal(size=(m, n))
                    b = 0.2*rng.normal(size=(m, n))
                    if complex_type:
                        a = a + 1j*rng.normal(size=(m, n))
                        b = b + 1j*0.2*rng.normal(size=(m, n))
                    if variant == 1:  # exact rank-one
                        u = rng.normal(size=(m, 1))
                        v = rng.normal(size=(1, n))
                        if complex_type:
                            u = u + 1j*rng.normal(size=(m, 1))
                            v = v + 1j*rng.normal(size=(1, n))
                        a = u @ v
                    if variant == 2:  # spectrally separated diagonal
                        a = np.zeros((m, n),dtype=a.dtype)
                        for i in range(min(m,n)):
                            a[i,i] = 5.0**(-i)
                    if variant == 3:
                        b = a.copy()  # exact zero difference
                    if variant == 4:
                        scale = 1e70 if (m+n) % 2 else 1e-70
                        a = a * scale
                        b = b * scale
                    start = rng.normal(size=n)
                    if complex_type:
                        start = start + 1j*rng.normal(size=n)
                    # Fortran and NumPy perform 30 complete power steps.
                    yield kind, m, n, 30, a, b, start


def numeric_value(z):
    return f"{float(np.real(z)):.18e} {float(np.imag(z)):.18e}"


def prepare(testcases):
    items = [str(len(testcases))]
    for kind, m, n, iterations, a, b, start in testcases:
        items.append(f"{kind} {m} {n} {iterations}")
        items.extend(numeric_value(z) for z in a.flat)
        items.extend(numeric_value(z) for z in b.flat)
        items.extend(numeric_value(z) for z in start.flat)
    return "\n".join(items) + "\n"


def direct_power(a, start, iterations):
    # Scale before products for large/tiny norms and complex matrices.
    if not a.size:
        return 0.0
    scale = float(np.max(np.abs(a)))
    if scale == 0.0:
        return 0.0
    a = a/scale
    x = start/np.linalg.norm(start)
    for _ in range(iterations):
        y = a @ x
        ny = np.linalg.norm(y)
        if ny == 0.0:
            return 0.0
        x = a.conj().T @ (y/ny)
        nx = np.linalg.norm(x)
        if nx == 0.0:
            return 0.0
        x = x/nx
    return float(np.linalg.norm(a @ x)) * scale


def main(driver):
    testcases = list(cases())
    completed = subprocess.run(
        [driver], input=prepare(testcases), text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True, timeout=180,
    )
    output = completed.stdout.splitlines()
    if len(output) != len(testcases):
        raise AssertionError((len(output),len(testcases),completed.stderr[:600]))
    max_power_diff = 0.0
    max_true_error = 0.0
    max_scipy_error = 0.0
    total_nontrivial = 0
    by_kind = [0]*4
    max_operator_diff = 0.0
    for icase, (args,line) in enumerate(zip(testcases,output),1):
        kind, m, n, iterations, a, b, start = args
        fields = line.split()
        if fields[:4] != ['CASE',str(icase),'0','0']:
            raise AssertionError(('bad Fortran result',icase,fields))
        fd, fo = (float(x) for x in fields[4:])
        mat = a if kind % 2 == 0 else a - b
        if not mat.size:
            expected = exact = scipy_result = 0.0
        else:
            expected = direct_power(mat,start,iterations)
            scale = max(float(np.max(np.abs(mat))),1e-150)
            exact = float(np.linalg.svd(mat/scale,compute_uv=False)[0])*scale
            scipy_result = float(interpolative.estimate_spectral_norm(
                mat,its=iterations,rng=np.random.default_rng(31+icase)))
        normscale = max(abs(exact),1e-150)
        rel_power = abs(fd-expected)/normscale
        rel_op = abs(fo-expected)/normscale
        rel_exact = abs(fd-exact)/normscale
        if rel_power > 6e-12 or rel_op > 6e-12:
            raise AssertionError((icase,kind,(m,n),'power',rel_power,rel_op,fd,fo,expected))
        # No randomized method produces a rigorous accuracy guarantee, but
        # a 30-iteration estimate must remain a lower bound (up to roundoff).
        if fd > exact+2e-12*normscale:
            raise AssertionError((icase,kind,'estimated norm exceeds true 2-norm',fd,exact))
        if exact > 0.0:
            # These small test matrices have enough power iterations for a
            # close approximation except for rare near-degenerate spectra.
            if rel_exact > 0.075:
                raise AssertionError((icase,kind,'poor power convergence',rel_exact))
            if not np.isfinite(scipy_result) or scipy_result > exact+1e-10*normscale:
                raise AssertionError((icase,kind,'unexpected SciPy estimate',scipy_result,exact))
            max_scipy_error = max(max_scipy_error,abs(scipy_result-exact)/normscale)
            total_nontrivial += 1
        max_power_diff = max(max_power_diff,rel_power)
        max_operator_diff = max(max_operator_diff,rel_op)
        max_true_error = max(max_true_error,rel_exact)
        by_kind[kind] += 1
    print(f'SciPy {scipy.__version__}: {len(testcases)} spectral norm comparison cases passed')
    print('Cases by category (real, real difference, complex, complex difference):',by_kind)
    print('Nonzero reference norms:',total_nontrivial)
    print(f'Max relative error vs independently computed power iteration: {max_power_diff:.3e}')
    print(f'Max relative error for matrix-free operator: {max_operator_diff:.3e}')
    print(f'Max relative error vs exact SVD norm: {max_true_error:.3e}')
    print(f'Max relative error of SciPy randomized reference vs exact: {max_scipy_error:.3e}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: compare_linalg_spectral_norm_scipy.py ./driver')
    main(sys.argv[1])
