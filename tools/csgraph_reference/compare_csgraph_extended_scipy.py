#!/usr/bin/env python3
"""Compare SciFort APSP, graph Laplacians, and DFS with SciPy 1.17.

Random cases include asymmetric, disconnected, signed, dense/sparse,
cyclic and acyclic graphs, and self-loops.  Compare predecessor path
*validity*, not the arbitrary choice between tied optimal paths.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import warnings

warnings.filterwarnings("ignore", message="Graph has negative weights: dijkstra.*", category=UserWarning)

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import (
    NegativeCycleError, depth_first_order, floyd_warshall,
    johnson, laplacian,
)

ROOT = Path(__file__).resolve().parents[2]


def cases() -> list[np.ndarray]:
    rng = np.random.default_rng(20261103)
    graphs = []
    for n in (1, 2, 3, 5, 8, 12, 18):
        for trial in range(30):
            a = np.zeros((n, n))
            if trial > 1:
                density = (.06, .2, .45, .8)[trial % 4]
                weights = rng.integers(-6 if trial % 3 == 0 else 1, 15, size=(n, n))
                weights[weights == 0] = 1
                a = np.where(rng.random((n, n)) < density, weights, 0).astype(float)
                if trial % 5 == 0:  # DAG, allowing negative edges but no cycles
                    a = np.triu(a, 1)
                if trial % 5 == 1:  # undirected graph
                    a = np.triu(abs(a), 1)
                    a += a.T
                np.fill_diagonal(a, 0)
                if trial % 10 == 7 and n > 1:  # planted negative 2-cycle
                    a[0, 1] = -4.
                    a[1, 0] = 1.
            graphs.append(a)
    return graphs


def check_paths(a: np.ndarray, distances: np.ndarray, pred: np.ndarray,
                *, directed: bool, unweighted: bool) -> None:
    n = len(a)
    assert pred.shape == (n, n)
    for source in range(n):
        assert pred[source, source] == 0
        for target in range(n):
            reachable = np.isfinite(distances[source, target])
            if not reachable or target == source:
                assert pred[source, target] == 0, (source, target, pred[source, target])
                continue
            seen = set()
            cost = 0.
            node = target
            while node != source:
                assert node not in seen, ('predecessor loop', source, target)
                seen.add(node)
                parent = int(pred[source, node]) - 1
                assert 0 <= parent < n
                edge_weight = a[parent, node]
                if not directed:
                    other = a[node, parent]
                    if other != 0 and (edge_weight == 0 or other < edge_weight):
                        edge_weight = other
                assert edge_weight != 0, ('missing edge', source, target, parent, node)
                cost += 1. if unweighted else edge_weight
                node = parent
            assert np.isclose(cost, distances[source, target], atol=1e-10, rtol=1e-10), (
                source, target, cost, distances[source, target]
            )


def check_depth_first(a: np.ndarray, order: np.ndarray, pred: np.ndarray,
                      directed: bool) -> None:
    n = len(a)
    assert len(pred) == n and order[0] == 1
    assert len(order) == len(set(order))
    ref_order, ref_pred = depth_first_order(csr_array(a), 0, directed=directed)
    assert set(order) == set(ref_order + 1)
    assert pred[0] == 0
    for node in range(n):
        if node == 0 or node + 1 not in order:
            assert pred[node] == 0
        else:
            parent = int(pred[node]) - 1
            assert parent + 1 in order
            assert a[parent, node] != 0 or (not directed and a[node, parent] != 0)
    if directed:
        # For row-sorted CSR both implementations choose the first unvisited arc.
        assert np.array_equal(order, ref_order + 1), (order, ref_order + 1)


def compare(mode: str, graphs: list[np.ndarray]) -> None:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-csgraph-extended-') as temp:
        build = Path(temp)
        cc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
              '-J', str(build), '-I', str(build)]
        objects = []
        for name in ('kinds', 'sparse_csgraph'):
            obj = build / (name + '.o')
            subprocess.run([*cc, '-c', str(ROOT/'src'/f'{name}.f90'),
                            '-o', str(obj)], cwd=build, check=True)
            objects.append(str(obj))
        exe = build/'compare_csgraph_extended'
        subprocess.run([*cc, str(Path(__file__).with_name('csgraph_extended_driver.f90')),
                        *objects, '-o', str(exe)], cwd=build, check=True)
        lines = [str(len(graphs))]
        for a in graphs:
            n = len(a)
            lines.append(str(n))
            lines.extend(' '.join(map(str, row)) for row in a)
        result = subprocess.run([str(exe)], cwd=build, input='\n'.join(lines)+'\n',
                                capture_output=True, text=True, check=True)
        output = result.stdout.splitlines()
        cursor = 0
        apsp_success = cycles = lap_pass = dfs_pass = 0
        for case_number, a in enumerate(graphs):
            n = int(output[cursor]); cursor += 1
            assert n == len(a)
            sparse = csr_array(a)
            for directed in (False, True):
                for unweighted in (False, True):
                    for algo in (floyd_warshall, johnson):
                        status = int(output[cursor]); cursor += 1
                        try:
                            ref = algo(sparse, directed=directed, unweighted=unweighted)
                        except NegativeCycleError:
                            assert status == 2, (mode, case_number, algo.__name__, 'missed negative cycle')
                            cycles += 1
                            continue
                        assert status == 0, (mode, case_number, algo.__name__, 'false negative cycle')
                        d = np.array([np.fromstring(line, sep=' ') for line in output[cursor:cursor+n]])
                        cursor += n
                        p = np.array([np.fromstring(line, sep=' ', dtype=int) for line in output[cursor:cursor+n]])
                        cursor += n
                        assert np.allclose(d, ref, atol=1e-10, rtol=1e-10), (
                            mode, case_number, algo.__name__, directed, unweighted, d, ref
                        )
                        if not (np.any(np.diag(a) != 0)):
                            check_paths(a, d, p, directed=directed, unweighted=unweighted)
                        apsp_success += 1
                size = int(output[cursor]); cursor += 1
                order = np.fromstring(output[cursor], dtype=int, sep=' '); cursor += 1
                pred = np.fromstring(output[cursor], dtype=int, sep=' '); cursor += 1
                assert len(order) == size
                check_depth_first(a, order, pred, directed)
                dfs_pass += 1
            for sym in (False, True):
                for outgoing in (False, True):
                    for normalize in (False, True):
                        status = int(output[cursor]); cursor += 1
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore', RuntimeWarning)
                            ref_lap, ref_diag = laplacian(a, normed=normalize,
                                return_diag=True, use_out_degree=outgoing, symmetrized=sym)
                        if np.any(~np.isfinite(ref_lap)) or np.any(~np.isfinite(ref_diag)):
                            assert status == 1
                            continue
                        assert status == 0, (mode, case_number, sym, outgoing, normalize)
                        diag = np.fromstring(output[cursor], sep=' '); cursor += 1
                        lap = np.array([np.fromstring(line, sep=' ') for line in output[cursor:cursor+n]])
                        cursor += n
                        assert np.allclose(lap, ref_lap, rtol=1e-10, atol=1e-10), (
                            mode, case_number, 'laplacian', sym, outgoing, normalize, lap, ref_lap
                        )
                        assert np.allclose(diag, ref_diag, rtol=1e-10, atol=1e-10), (
                            mode, case_number, 'degree', sym, outgoing, normalize, diag, ref_diag
                        )
                        lap_pass += 1
        assert cursor == len(output), (cursor, len(output))
    print(f'{mode}: {len(graphs)} graphs; {apsp_success} all-pairs comparisons; '
          f'{cycles} negative cycles; {lap_pass} Laplacians; {dfs_pass} DFS traversals PASS', flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['checked','optimized','both'], default='both')
    args = parser.parse_args()
    graphs = cases()
    for mode in (('checked','optimized') if args.mode == 'both' else (args.mode,)):
        compare(mode, graphs)


if __name__ == '__main__':
    main()
