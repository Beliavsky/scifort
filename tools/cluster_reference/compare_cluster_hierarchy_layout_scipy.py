#!/usr/bin/env python3
"""Numerically cross-check SciFort hierarchy layout and predicates vs SciPy.

Run from anywhere; builds the standalone Fortran reference driver with gfortran,
using a temporary directory. Requires NumPy, SciPy, and gfortran, not FPM.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import dendrogram, is_isomorphic, linkage

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    ROOT / 'src' / 'kinds.f90',
    ROOT / 'src' / 'spatial_distance.f90',
    ROOT / 'src' / 'cluster_hierarchy.f90',
    ROOT / 'tools' / 'cluster_reference' / 'cluster_hierarchy_layout_reference_driver.f90',
)
METHODS = ('single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward')


def run(ncases: int, optimization: str) -> None:
    rng = np.random.default_rng(312865)
    data = []
    headers = []
    ref = []
    scipy_label_mismatches = 0
    for case in range(ncases):
        n = int(rng.integers(2, 29))
        mode = case % 5
        # Include non-monotone trees, ties, and negative observations.
        x = rng.normal(size=(n, int(rng.integers(1, 5))))
        if case % 11 == 0:
            x = np.round(x, 1)
        z0 = linkage(x, method=METHODS[case % len(METHODS)])
        z = z0.copy()
        z[:, :2] += 1.0
        count = int(rng.integers(1, min(n, 9) + 1))
        t1 = rng.integers(-count, count + 1, size=n, dtype=np.int32)
        if case % 3 == 0:
            names = np.unique(t1)
            replacements = rng.choice(np.arange(-1000, 1000), len(names), replace=False)
            t2 = np.array([replacements[np.searchsorted(names, label)] for label in t1], dtype=np.int32)
        else:
            t2 = rng.integers(-count, count + 1, size=n, dtype=np.int32)
        # Compare against the mathematical equivalence relation, not just
        # SciPy's current adjacent-boundary test: SciPy 1.17 can yield false
        # positives when a label repeats nonadjacently after sorting by T1.
        expected_iso = int(np.array_equal(t1[:, None] == t1[None, :],
                                          t2[:, None] == t2[None, :]))
        if int(is_isomorphic(t1, t2)) != expected_iso:
            scipy_label_mismatches += 1
        kwargs = (
            {'count_sort': 'ascending'} if mode == 1 else
            {'count_sort': 'descending'} if mode == 2 else
            {'distance_sort': 'ascending'} if mode == 3 else
            {'distance_sort': 'descending'} if mode == 4 else {}
        )
        expected = dendrogram(z0, no_plot=True, **kwargs)
        lines = [f'{n} {mode}']
        lines.extend(' '.join(f'{v:.17g}' for v in row) for row in z)
        lines.append(' '.join(str(i) for i in t1))
        lines.append(' '.join(str(i) for i in t2))
        data.extend(lines)
        headers.append(n)
        ref.append((case, expected, expected_iso))

    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp) / 'cluster_hierarchy_layout_reference_driver'
        compiler = shutil.which('gfortran')
        if compiler is None:
            raise RuntimeError('gfortran is required for this comparison')
        cmd = [compiler, '-std=f2018', optimization, '-ffree-line-length-none',
               '-J', tmp, *map(str, SOURCES), '-o', str(exe)]
        if optimization == '-O0':
            cmd += ['-g', '-fcheck=all', '-fbacktrace']
        subprocess.run(cmd, cwd=tmp, check=True, capture_output=True, text=True)
        result = subprocess.run([str(exe)], input='\n'.join(data)+'\n',
                                text=True, capture_output=True, check=True)
        output = iter(result.stdout.splitlines())
    for n, (case, expected, expected_iso) in zip(headers, ref):
        status, iso, converted = map(int, next(output).split())
        assert (status, iso, converted) == (0, expected_iso, 1), (case,status,iso,converted,expected_iso)
        leaves = np.fromstring(next(output), dtype=int, sep=' ')
        np.testing.assert_array_equal(leaves, np.asarray(expected['leaves']) + 1,
                                      err_msg=f'case {case} leaves')
        coords = np.stack([np.fromstring(next(output), dtype=float, sep=' ')
                           for _ in range(n-1)])
        np.testing.assert_allclose(coords[:, :4], expected['icoord'], rtol=1e-14,
                                   atol=1e-12, err_msg=f'case {case} x-coordinates')
        np.testing.assert_allclose(coords[:, 4:], expected['dcoord'], rtol=1e-14,
                                   atol=1e-12, err_msg=f'case {case} y-coordinates')
    assert next(output, None) is None, 'unexpected output from reference driver'
    print(f'PASS: {ncases} SciPy dendrogram + rigorous is_isomorphic + MATLAB roundtrips ({optimization}); 'f'{scipy_label_mismatches} SciPy 1.17 false-positive partition equivalences')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', type=int, default=250)
    parser.add_argument('--optimization', default='-O0', choices=['-O0','-O3'])
    args = parser.parse_args()
    run(args.cases, args.optimization)
