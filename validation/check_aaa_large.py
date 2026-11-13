#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check public AAA spectra and residues from larger retained representations."""
import argparse
from pathlib import Path
import subprocess

import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import linear_sum_assignment

from check_aaa_spectrum import matching, roots


def secular_roots(s, u):
    n = len(s)+1
    a = np.zeros((n, n), complex)
    b = a.copy()
    a[0, 1:] = u/max(abs(u))
    a[1:, 0] = 1
    a[np.arange(1, n), np.arange(1, n)] = s
    b[np.arange(1, n), np.arange(1, n)] = 1
    return roots(*eigvals(a, b, homogeneous_eigvals=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    cases = []
    for n in (41, 65, 81, 97):
        t = np.linspace(-1, 1, n)
        for mode in (0, 1):
            x = t.astype(complex)
            f = (1.1+np.sin(1000*t)+0.3*np.cos(747*t)).astype(complex)
            if mode:
                x += 0.2j*np.sin(3*t)
                f = np.sin(239*np.arange(1, n+1))+1j*np.cos(137*np.arange(1, n+1))
            for clean in (0, 1, 2):
                cases.append((x, f, mode, clean))
    lines = []
    for x, f, mode, clean in cases:
        lines += [f"{len(x)} {mode} {clean}", "0.1"]
        lines += [f"{a.real:.17g} {a.imag:.17g} {b.real:.17g} {b.imag:.17g}"
                  for a, b in zip(x, f)]
    result = subprocess.run([str(args.executable.resolve())],input="\n".join(lines)+"\n",
                            capture_output=True,text=True,check=True)
    rows = iter(result.stdout.splitlines())
    maximum = residue_max = 0.0
    spectral = residues = removed_total = large = 0
    for x, _, mode, clean in cases:
        status, removed, m, ps, zs, rs, np_, nz, nr = map(int,next(rows).split())
        assert status == ps == zs == rs == 0, (len(x), mode, clean, status, ps, zs, rs)
        data = np.array([np.fromstring(next(rows),sep=" ") for _ in range(m)])
        s = data[:,0]+1j*data[:,1]
        v = data[:,2]+1j*data[:,3]
        w = data[:,4]+1j*data[:,5]
        arrays = []
        for count in (np_, nz, nr):
            a = np.array([np.fromstring(next(rows),sep=" ") for _ in range(count)])
            arrays.append(a[:,0]+1j*a[:,1] if count else np.array([],complex))
        p, z, r = arrays
        expected_p = secular_roots(s,w)
        expected_z = secular_roots(s,w*v)
        pe, ze = matching(p,expected_p),matching(z,expected_z)
        assert pe < 3e-8 and ze < 3e-8, (len(x),mode,clean,pe,ze)
        d = expected_p[:,None]-s
        expected_r = np.sum(w*v/d,axis=1)/(-np.sum(w/d**2,axis=1))
        i,j = linear_sum_assignment(abs(p[:,None]-expected_p[None,:]))
        re = np.max(abs(r[i]-expected_r[j])/(1+abs(expected_r[j])))
        assert re < 3e-8, (len(x),mode,clean,re)
        maximum = max(maximum,pe,ze)
        residue_max = max(residue_max,re)
        spectral += len(p)+len(z)
        residues += len(r)
        removed_total += removed
        if m > 39:
            large += 1
    assert next(rows,None) is None
    assert large >= 16
    assert removed_total > 0
    print(f"PASS: {len(cases)} public real/complex fits, {large} above 39 supports; "
          f"{spectral} spectral values, {residues} residues; {removed_total} supports removed")
    print(f"Maximum normalized errors: spectra {maximum:.6g}, residues {residue_max:.6g}")


if __name__ == "__main__":
    main()
