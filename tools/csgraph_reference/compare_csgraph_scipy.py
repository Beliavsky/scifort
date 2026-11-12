#!/usr/bin/env python3
"""Randomized reference checks against SciPy 1.17's sparse.csgraph algorithms.

The comparison checks partitions (rather than arbitrary numeric labels) and
minimum-spanning-forest cost (rather than an arbitrary choice of tied edges).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components, dijkstra, minimum_spanning_tree

ROOT = Path(__file__).resolve().parents[2]


def same_partition(a: np.ndarray, b: np.ndarray) -> bool:
    return bool(np.array_equal(a[:, None] == a[None, :],
                               b[:, None] == b[None, :]))


def cases() -> list[np.ndarray]:
    rng = np.random.default_rng(20261103)
    result = []
    for n in (1, 2, 3, 5, 8, 12, 18):
        for trial in range(40):
            if trial in (0, 1, 2):
                a = np.zeros((n, n), dtype=np.float64)
            else:
                weights = rng.integers(-20 if trial % 10 == 0 else 1,
                                       35, size=(n, n))
                weights[weights == 0] = 1
                density = (0.07, 0.2, 0.55, 0.9)[trial % 4]
                a = np.where(rng.random((n, n)) < density, weights, 0).astype(float)
                if trial % 3 == 0:
                    a = np.triu(a, 1)
                    a += a.T
            np.fill_diagonal(a, 0)
            result.append(a)
    return result


def compare(mode: str, graphs: list[np.ndarray]) -> None:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-csgraph-reference-') as temp:
        build = Path(temp)
        cc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
              '-J', str(build), '-I', str(build)]
        objs = []
        for name in ('kinds', 'sparse_csgraph'):
            obj = build/(name+'.o')
            subprocess.run([*cc, '-c', str(ROOT/'src'/(name+'.f90')),
                            '-o', str(obj)], cwd=build, check=True)
            objs.append(str(obj))
        exe = build/'csgraph_driver'
        subprocess.run([*cc, str(Path(__file__).with_name('csgraph_reference_driver.f90')),
                        *objs, '-o', str(exe)], cwd=build, check=True)
        lines = [str(len(graphs))]
        for a in graphs:
            lines.append(str(a.shape[0]))
            lines += [' '.join(str(v) for v in row) for row in a]
        output = subprocess.run([str(exe)], cwd=build, input='\n'.join(lines)+'\n',
                                check=True, text=True, capture_output=True).stdout.splitlines()
        idx = 0
        checked_distance_rows = 0
        for trial, a in enumerate(graphs):
            graph = csr_array(a)
            n = a.shape[0]
            nr, nw, ns = map(int, output[idx].split())
            idx += 1
            assert nr == n
            weak = np.fromstring(output[idx], dtype=int, sep=' ')
            idx += 1
            strong = np.fromstring(output[idx], dtype=int, sep=' ')
            idx += 1
            tree_edges, tree_weight = output[idx].split()
            idx += 1
            tree_edges, tree_weight = int(tree_edges), float(tree_weight)
            ref_nw, ref_weak = connected_components(graph, directed=True, connection='weak')
            ref_ns, ref_strong = connected_components(graph, directed=True, connection='strong')
            assert (nw, ns) == (ref_nw, ref_ns), (mode, trial, 'components count')
            assert same_partition(weak, ref_weak), (mode, trial, 'weak partitions')
            assert same_partition(strong, ref_strong), (mode, trial, 'strong partitions')
            tree = minimum_spanning_tree(graph)
            assert tree_edges == tree.nnz, (mode, trial, 'tree edge count')
            assert abs(tree_weight-float(tree.sum())) < 1.e-8, (mode, trial, 'tree cost')
            mode_str = output[idx].strip()
            idx += 1
            if np.any(a < 0):
                assert mode_str == 'NEGATIVE'
                continue
            assert mode_str == 'POSITIVE'
            direct = np.array([np.fromstring(line, sep=' ') for line in output[idx:idx+n]])
            idx += n
            undirect = np.array([np.fromstring(line, sep=' ') for line in output[idx:idx+n]])
            idx += n
            ref_direct = dijkstra(graph, directed=True)
            assert np.allclose(direct, ref_direct, rtol=1.e-11, atol=1.e-11), (mode,trial,'dijkstra')
            checked_distance_rows += n
            # SciPy warns of undefined behavior for asymmetric reciprocal
            # arc weights when directed=False; compare only symmetric graphs.
            if np.array_equal(a, a.T):
                ref_undirect = dijkstra(graph, directed=False)
                assert np.allclose(undirect, ref_undirect, rtol=1.e-11, atol=1.e-11), (mode,trial,'undirected')
        assert idx == len(output), (idx,len(output))
    print(f'{mode}: {len(graphs)} graph cases, {checked_distance_rows} directed source-path rows PASS')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked', 'optimized', 'both'), default='both')
    args = parser.parse_args()
    for mode in (('checked', 'optimized') if args.mode == 'both' else (args.mode,)):
        compare(mode, cases())


if __name__ == '__main__':
    main()
