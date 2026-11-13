#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check complex Y0/Y1/Y1-prime zeros against SciPy and 80-digit mpmath."""
import argparse
from pathlib import Path
import subprocess

import mpmath as mp
import numpy as np
import scipy
from scipy.special import y0_zeros, y1_zeros, y1p_zeros


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--count", type=int, default=160)
    args = parser.parse_args()
    assert 16 <= args.count <= 1000
    mp.mp.dps = 80
    oracles = (y0_zeros, y1_zeros, y1p_zeros)
    root_error = companion_error = residual = mp_root_error = mp_value_error = 0.0
    checked_mp = 0
    scipy_outliers = 0
    for family, oracle in enumerate(oracles):
        result = subprocess.run([str(args.executable.resolve())],
            input=f"{family} {args.count}\n", text=True, capture_output=True, check=True)
        lines = result.stdout.splitlines()
        assert list(map(int, lines[0].split())) == [0, args.count], (family, lines[0])
        table = np.array([list(map(float, row.split())) for row in lines[1:]])
        assert table.shape == (args.count, 4)
        got = table[:, 0]+1j*table[:, 1]
        values = table[:, 2]+1j*table[:, 3]
        assert np.all(np.isfinite(got)) and np.all(np.isfinite(values))
        assert np.all(got.imag > 0) and np.all(np.diff(got.real) < 0)
        if family == 2:
            assert got[0].real > 0 and np.all(got[1:].real < 0)
        else:
            assert np.all(got.real < 0)
        # SciPy 1.15.3 returns nonfinite roots beyond index 140. Compare its
        # supported prefix explicitly; independent mpmath checks cover the tail.
        scipy_count = min(args.count, 140)
        ref, ref_values = oracle(scipy_count, complex=True)
        # SciPy returns +Y1 for Y0; SciFort retains the actual derivative -Y1.
        if family == 0:
            ref_values = -ref_values
        rerr = np.abs(got[:scipy_count]-ref)/(1+np.abs(ref))
        verr = np.abs(values[:scipy_count]-ref_values)/np.abs(ref_values)
        suspects = set(np.flatnonzero((rerr > 5e-12) | (verr > 2e-10)).tolist())
        root_error = max(root_error, float(np.max(rerr)))
        companion_error = max(companion_error, float(np.max(verr)))
        # Each high-precision solve starts from the mathematical phase seed,
        # independently of either library's computed root.
        selected = sorted(set(range(16)) | suspects | {args.count//4-1, args.count//2-1,
                                                      3*args.count//4-1, args.count-1})
        for j in selected:
            k = j+1
            offset = (mp.mpf("0.25"), mp.mpf("0.75"), mp.mpf("1.25"))[family]
            seed = -(k-offset)*mp.pi+mp.j*mp.log(3)/2
            if family == 1 and k == 1:
                seed = mp.mpc("-0.5", "0.5")
            if family == 2 and k == 1:
                seed = mp.mpc("0.5", "1.0")
            def equation(z):
                if family < 2:
                    return mp.bessely(family, z)
                return mp.bessely(0, z)-mp.bessely(1, z)/z
            expected = mp.findroot(equation, seed, tol=mp.mpf("1e-65"))
            z = mp.mpc(float(got[j].real), float(got[j].imag))
            if family == 0:
                value = -mp.bessely(1, expected)
            elif family == 1:
                value = mp.bessely(0, expected)-mp.bessely(1, expected)/expected
            else:
                value = mp.bessely(1, expected)
            re = float(abs(z-expected)/(1+abs(expected)))
            ve = float(abs(mp.mpc(complex(values[j]))-value)/abs(value))
            res = float(abs(equation(z))*mp.sqrt(abs(z)))
            assert re <= 5e-13, (family, k, re)
            assert ve <= 3e-11, (family, k, ve)
            assert res <= 2e-12*(1+abs(z)), (family, k, res)
            if j in suspects:
                # Any disagreement with SciPy must be adjudicated using the
                # independent root, never by relaxing the native accuracy bound.
                scipy_re = float(abs(mp.mpc(complex(ref[j]))-expected)/(1+abs(expected)))
                scipy_ve = float(abs(mp.mpc(complex(ref_values[j]))-value)/abs(value))
                assert scipy_re > 5e-12 or scipy_ve > 2e-10, (family, k, re, ve)
                scipy_outliers += 1
            mp_root_error = max(mp_root_error, re)
            mp_value_error = max(mp_value_error, ve)
            residual = max(residual, res)
            # The lower conjugate is also a root on the principal branch.
            assert abs(equation(mp.conj(z))) <= mp.mpf("2e-12")*(1+abs(z))/mp.sqrt(abs(z))
            checked_mp += 1
        print(f"family {family}: {args.count} native roots, {scipy_count} SciPy pairs, "
              f"{len(selected)} mpmath pairs PASS",
              flush=True)
    print(f"SciPy {scipy.__version__}, mpmath {mp.__version__} (80 digits)")
    print(f"SciPy: max normalized root error {root_error:.6g}; "
          f"max relative companion error {companion_error:.6g}")
    print(f"SciPy discrepancies independently resolved by mpmath: {scipy_outliers}")
    print(f"mpmath: {checked_mp} pairs; max normalized root error {mp_root_error:.6g}; "
          f"max relative companion error {mp_value_error:.6g}; "
          f"max |residual|*sqrt(|z|) {residual:.6g}")


if __name__ == "__main__":
    main()
