#!/usr/bin/env python3
"""SciPy comparisons for Bellman-Ford and breadth-first graph traversals."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import (bellman_ford, breadth_first_order,
                                  NegativeCycleError)
from compare_csgraph_scipy import cases as base_cases

ROOT = Path(__file__).resolve().parents[2]


def cases() -> list[np.ndarray]:
    random = np.random.default_rng(92317)
    graphs = base_cases()
    for n in (2, 3, 7, 14, 20):
        for iteration in range(25):
            w = random.integers(-15, 17, size=(n, n)).astype(float)
            w[w == 0] = 1
            # Directed DAG: negative edges, but no directed cycles.
            w = np.where(random.random((n, n)) < (.1+0.7*(iteration%4)/3), w, 0)
            graphs.append(np.triu(w, k=1))
    return graphs


def check_bfs(a: np.ndarray, order: np.ndarray, pred: np.ndarray,
              directed: bool) -> None:
    n = a.shape[0]
    graph = csr_array(a)
    ref_order, ref_pred = breadth_first_order(graph, 0, directed=directed)
    assert set(order) == set(ref_order+1), ('BFS visited',a,order,ref_order)
    assert len(order) == len(set(order)) and order[0] == 1
    depth = np.full(n, np.inf)
    depth[0] = 0
    assert pred[0] == 0
    for v in order[1:]:
        previous = pred[v-1]
        assert previous != 0
        assert (a[previous-1,v-1] != 0 or
                (not directed and a[v-1,previous-1] != 0)), ('BFS edge',a,v,previous)
        assert np.isfinite(depth[previous-1]), ('BFS ancestor order',a)
        depth[v-1] = depth[previous-1]+1
    ref_depth = np.full(n, np.inf)
    ref_depth[0] = 0
    for v in ref_order[1:]:
        ref_depth[v] = ref_depth[ref_pred[v]] + 1
    assert np.array_equal(depth,ref_depth), ('BFS path lengths',a,depth,ref_depth)


def compare(mode: str, graphs: list[np.ndarray]) -> None:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-csgraph-paths-') as temp:
        build = Path(temp)
        cc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
              '-I', str(build), '-J', str(build)]
        objs = []
        for name in ('kinds','sparse_csgraph'):
            obj = build/(name+'.o')
            subprocess.run([*cc, '-c', str(ROOT/'src'/(name+'.f90')),
                            '-o', str(obj)], cwd=build, check=True)
            objs.append(str(obj))
        exe = build/'paths_driver'
        subprocess.run([*cc, str(Path(__file__).with_name('csgraph_paths_driver.f90')),
                        *objs, '-o', str(exe)], cwd=build, check=True)
        input_lines = [str(len(graphs))]
        for a in graphs:
            input_lines.append(str(len(a)))
            input_lines.extend(' '.join(str(x) for x in row) for row in a)
        output = subprocess.run([str(exe)], cwd=build, input='\n'.join(input_lines)+'\n',
                                check=True, text=True, capture_output=True).stdout.splitlines()
        cursor = 0
        checked_sources = 0
        negative_cycles = 0
        for number,a in enumerate(graphs):
            n = int(output[cursor]); cursor += 1
            assert n == len(a)
            graph = csr_array(a)
            for i in range(n):
                line = output[cursor].strip()
                cursor += 1
                try:
                    reference = bellman_ford(graph,indices=i)
                except NegativeCycleError:
                    assert line == 'NEGATIVE_CYCLE', (number,i,'cycle not detected')
                    negative_cycles += 1
                else:
                    actual = np.fromstring(line,sep=' ')
                    assert np.allclose(actual,reference,atol=1.e-11,rtol=1.e-11), (number,i,actual,reference)
                    checked_sources += 1
            for directed in (True,False):
                nr = int(output[cursor]); cursor += 1
                order = np.fromstring(output[cursor],dtype=int,sep=' '); cursor += 1
                pred = np.fromstring(output[cursor],dtype=int,sep=' '); cursor += 1
                assert len(order) == nr and len(pred) == n
                check_bfs(a,order,pred,directed)
        assert cursor == len(output)
    print(f'{mode}: {len(graphs)} graphs, {checked_sources} Bellman-Ford sources, '
          f'{negative_cycles} negative-cycle detections, {len(graphs)*2} BFS traversals PASS')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('checked','optimized','both'),default='both')
    args = parser.parse_args()
    for mode in (('checked','optimized') if args.mode == 'both' else (args.mode,)):
        compare(mode,cases())


if __name__ == '__main__':
    main()
