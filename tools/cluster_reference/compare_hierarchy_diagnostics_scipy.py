#!/usr/bin/env python3
"""Randomized SciPy 1.17 reference for hierarchical diagnostics and tree cuts.

Uses only the public SciPy Python API as a numerical oracle, not its source.
The Fortran program consumes a SciPy-generated linkage matrix with node IDs
shifted by one, ensuring comparisons isolate the new routines rather than the
agglomerative merge constructor.  Cluster labels are compared by partition.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.cluster import hierarchy as h

ROOT = Path(__file__).resolve().parents[2]
METHODS = ("single", "complete", "average", "weighted", "centroid", "median", "ward")


def compile_driver(folder: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    base = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
            '-J', str(folder), '-I', str(folder)]
    objects = []
    for src in ('kinds', 'spatial_distance', 'cluster_hierarchy'):
        obj = folder / (src + '.o')
        subprocess.run([*base, '-c', str(ROOT / 'src' / (src + '.f90')),
                        '-o', str(obj)], cwd=folder, check=True, capture_output=True)
        objects.append(str(obj))
    exe = folder / 'cluster_diagnostics_driver'
    subprocess.run([*base, str(ROOT / 'tools/cluster_reference' /
                                     'cluster_hierarchy_diagnostics_driver.f90'),
                    *objects, '-o', str(exe)], cwd=folder, check=True, capture_output=True)
    return exe


def same_partition(first, second):
    first = np.asarray(first)
    second = np.asarray(second)
    return np.array_equal(first[:, None] == first[None, :],
                          second[:, None] == second[None, :])


def same_cut_columns(first, second, description):
    assert first.shape == second.shape, (description, first.shape, second.shape)
    for i in range(first.shape[1]):
        if not same_partition(first[:, i], second[:, i]):
            raise AssertionError(f'{description}: column {i}, {first[:,i]} vs {second[:,i]}')


def run_case(exe, z_py, depth, incon_cut, mono_cut, k, heights):
    n = len(z_py) + 1
    z = z_py.copy()
    z[:, :2] += 1
    lines = [f'{n} {depth} {incon_cut:.17e} {mono_cut:.17e} {k}']
    lines += [' '.join(f'{v:.17e}' for v in row) for row in z]
    lines += [' '.join(f'{v:.17e}' for v in heights)]
    proc = subprocess.run([str(exe)], input='\n'.join(lines) + '\n', text=True,
                          capture_output=True, check=True)
    output = proc.stdout.splitlines()
    cursor = 0

    def get_float(rows):
        nonlocal cursor
        values = np.array([np.fromstring(s, sep=' ') for s in output[cursor:cursor + rows]])
        cursor += rows
        return values

    def get_int(rows):
        nonlocal cursor
        values = np.array([np.fromstring(s, sep=' ', dtype=np.int64)
                           for s in output[cursor:cursor + rows]], dtype=np.int64)
        cursor += rows
        return values

    expected_r = h.inconsistent(z_py, d=depth)
    obtained_r = get_float(n - 1)
    np.testing.assert_allclose(obtained_r, expected_r, rtol=5e-9, atol=2e-9,
                               err_msg='inconsistent')
    dist = h.maxdists(z_py)
    np.testing.assert_allclose(get_float(1)[0], dist, rtol=1e-13, atol=1e-13,
                               err_msg='maxdists')
    for j in range(4):
        expected_stat = h.maxRstat(z_py, expected_r, j)
        np.testing.assert_allclose(get_float(1)[0], expected_stat,
                                   rtol=5e-9, atol=2e-9,
                                   err_msg=f'maxRstat col {j}')
    np.testing.assert_allclose(get_float(1)[0], h.maxinconsts(z_py, expected_r),
                               rtol=5e-9, atol=2e-9, err_msg='maxinconsts')
    expected_inc = h.fcluster(z_py, incon_cut, criterion='inconsistent',
                              depth=depth, R=expected_r)
    assert same_partition(get_int(1)[0], expected_inc), 'fcluster inconsistent'
    assert same_partition(get_int(1)[0], expected_inc), 'fcluster inconsistent with R'
    expected_mono = h.fcluster(z_py, mono_cut, criterion='monocrit', monocrit=dist)
    assert same_partition(get_int(1)[0], expected_mono), 'fcluster monocrit'
    expected_k = h.fcluster(z_py, k, criterion='maxclust_monocrit', monocrit=dist)
    assert same_partition(get_int(1)[0], expected_k), 'fcluster maxclust_monocrit'
    full_cuts = get_int(n)
    count_cuts = get_int(n)
    # SciPy sorts nonmonotone centroid/median merges by *distance* when
    # constructing cut_tree. SciFort follows linkage-row order instead,
    # preserving the hierarchy's actual agglomeration sequence.
    if h.is_monotonic(z_py):
        same_cut_columns(full_cuts, h.cut_tree(z_py), 'cut_tree full')
        expected_counts = np.column_stack([h.cut_tree(z_py, n_clusters=[count])[:, 0]
                                           for count in (1, k, n, k)])
        same_cut_columns(count_cuts, expected_counts,
                         'cut_tree selected cluster counts')
        # Comparing scalar cuts separately avoids an upstream SciPy bug
        # where unsorted selections can leave array columns uninitialized.
        expected_heights = np.column_stack([h.cut_tree(z_py, height=[height])[:, 0]
                                            for height in heights])
        same_cut_columns(get_int(n), expected_heights, 'cut_tree height')
    else:
        assert output[cursor].strip() == 'NONMONOTONE', output[cursor]
        cursor += 1
    assert cursor == len(output), f'unexpected output trailing lines: {output[cursor:]}'
    return (n - 1) * 4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', type=int, default=140)
    parser.add_argument('--seed', type=int, default=20261103)
    parser.add_argument('--mode', choices=['checked', 'optimized', 'both'], default='both')
    args = parser.parse_args()
    modes = ('checked', 'optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        rng = np.random.default_rng(args.seed)
        total_stats = 0
        with tempfile.TemporaryDirectory(prefix='scifort-cluster-diag-') as folder_name:
            exe = compile_driver(Path(folder_name), mode)
            for case in range(args.cases):
                n = int(rng.integers(6, 22))
                ndim = int(rng.integers(1, 6))
                x = rng.normal(size=(n, ndim)) * rng.uniform(0.25, 3.0)
                method = METHODS[case % len(METHODS)]
                z = h.linkage(x, method=method)
                depth = [0, 1, 2, 3, 5, n + 2][case % 6]
                r = h.inconsistent(z, d=depth)
                dist = h.maxdists(z)
                # Move strictly away from a candidate threshold to avoid
                # rounding-dependent ties across compilers.
                incon_cut = float(np.quantile(r[:, 3], 0.3 + 0.05 * (case % 5)) + 1e-8)
                mono_cut = float(np.quantile(dist, 0.25 + 0.1 * (case % 5)) + 1e-8)
                k = int(rng.integers(2, n - 1))
                first_height, last_height = np.min(z[:, 2]), np.max(z[:, 2])
                heights = [float(first_height - 1e-7),
                           float(first_height + 1e-8),
                           float(np.median(z[:, 2]) + 1e-8),
                           float(last_height + 1e-7)]
                try:
                    total_stats += run_case(exe, z, depth, incon_cut, mono_cut, k, heights)
                except (AssertionError, subprocess.CalledProcessError) as exc:
                    raise AssertionError(f'case={case} method={method} depth={depth} '
                                         f'n={n}: {exc}') from exc
        print(f'{mode}: {args.cases} randomized SciPy hierarchy-diagnostics '
              f'comparisons PASS; {total_stats} inconsistency entries', flush=True)


if __name__ == '__main__':
    main()
