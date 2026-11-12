#!/usr/bin/env python3
"""Randomized SciPy comparisons for RCM, max flow, and bipartite matching.

The Fortran source is independent of SciPy; this file is test-only.
Matching permutations and max-flow decompositions may differ; compare their
mathematical invariants and optimal objective values, not arbitrary tie breaks.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_array
from scipy.sparse.csgraph import (maximum_flow, maximum_bipartite_matching,
                                   reverse_cuthill_mckee, structural_rank)

ROOT = Path(__file__).resolve().parents[2]


def solve_bruteforce_matching(nr, nc, pairs):
    adj = [set() for _ in range(nr)]
    for r, c in pairs:
        adj[r].add(c)
    states = {0}
    for opts in adj:
        states = states | {mask | (1 << c) for mask in states for c in opts if not mask & (1 << c)}
    return max((mask.bit_count() for mask in states), default=0)


def solve_bruteforce_cut(n, edges, s, t):
    # Independent exact min-cut enumeration for small graphs.
    mid = [v for v in range(n) if v not in (s, t)]
    best = 10**30
    for mask in range(1 << len(mid)):
        on_source = {s} | {v for k, v in enumerate(mid) if (mask >> k) & 1}
        cut = sum(w for u, v, w in edges if u in on_source and v not in on_source)
        best = min(best, cut)
    return best


def bandwidth(n, edges, perm):
    inv = np.empty(n, dtype=int)
    inv[np.asarray(perm, dtype=int)] = np.arange(n)
    return max((abs(int(inv[u]) - int(inv[v])) for u, v, _ in edges), default=0)


def cases_random(seed=41924):
    rng = np.random.default_rng(seed)
    cases = []
    for i in range(170):
        n = int(rng.integers(2, 25))
        m = int(rng.integers(0, min(4 * n + 1, n * n + 25)))
        edges = [(int(rng.integers(n)), int(rng.integers(n)), int(rng.integers(0, 19)))
                 for _ in range(m)]
        s, t = rng.choice(n, size=2, replace=False)
        cases.append(('F', n, int(s), int(t), edges))
    for i in range(230):
        nr = int(rng.integers(0, 40))
        nc = int(rng.integers(0, 42))
        m = int(rng.integers(0, (nr * nc // 4) + 1)) if nr * nc else 0
        pairs = {(int(rng.integers(nr)), int(rng.integers(nc))) for _ in range(m)}
        edges = [(r, c, 1) for r, c in sorted(pairs)]
        cases.append(('M', nr, nc, 0, edges))
    for i in range(45):
        nr = int(rng.integers(0, 8))
        nc = int(rng.integers(0, 8))
        edges = [(r, c, 0) for r in range(nr) for c in range(nc) if rng.random() < .4]
        cases.append(('M', nr, nc, 1, edges))
    for i in range(175):
        n = int(rng.integers(0, 81))
        symmetric = i % 3 == 0
        edges = []
        if n:
            # Randomly relabel a sparse banded graph, then possibly add loops.
            labels = rng.permutation(n)
            for u in range(n):
                for offset in (1, 2, 3):
                    if u + offset >= n or rng.random() > (.70 if offset == 1 else .28):
                        continue
                    v = u + offset
                    w = int(rng.integers(-9, 10))
                    edges.append((int(labels[u]), int(labels[v]), w))
                    if symmetric:
                        edges.append((int(labels[v]), int(labels[u]), w))
                if rng.random() < .15:
                    edges.append((int(labels[u]), int(labels[u]), 0))
        cases.append(('R', n, int(symmetric), 0, edges))
    rng.shuffle(cases)
    return cases


def matrix_from_edges(nr, nc, edges, dtype=np.int64):
    if edges:
        rows, cols, vals = zip(*edges)
    else:
        rows, cols, vals = (), (), ()
    return coo_array((np.asarray(vals, dtype=dtype),
                      (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32))),
                     shape=(nr, nc)).tocsr()


def run(mode, seed=41924):
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    with tempfile.TemporaryDirectory(prefix='scifort-advanced-') as dirname:
        build = Path(dirname)
        compiler = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
                    '-I', str(build), '-J', str(build)]
        for module in ('kinds', 'sparse_csgraph'):
            subprocess.run([*compiler, '-c', str(ROOT/'src'/f'{module}.f90'),
                            '-o', str(build/f'{module}.o')], check=True, cwd=build)
        exe = build/'csgraph_advanced_driver'
        subprocess.run([*compiler,
                        str(ROOT/'tools/csgraph_reference/csgraph_advanced_driver.f90'),
                        *(str(build/f'{m}.o') for m in ('kinds', 'sparse_csgraph')),
                        '-o', str(exe)], check=True, cwd=build)
        cases = cases_random(seed)
        lines = [str(len(cases))]
        for task, n, a1, a2, edges in cases:
            first, second = ((a1 + 1, a2 + 1) if task == 'F' else (a1, a2))
            lines.append(f'{task} {n} {len(edges)} {first} {second}')
            lines += [f'{u+1} {v+1} {w}' for u, v, w in edges]
        result = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                                text=True, capture_output=True, check=True)
        output = iter(result.stdout.splitlines())
        counts = {'F': 0, 'M': 0, 'R': 0, 'brute_flow': 0, 'brute_match': 0}
        for task, n, a1, a2, edges in cases:
            matrix = matrix_from_edges(n, a1 if task == 'M' else n, edges)
            if task == 'F':
                status, got = map(int, next(output).split())
                assert status == 0, (task, status, edges)
                flow = np.asarray([[int(z) for z in next(output).split()] for _ in range(n)])
                assert flow.shape == (n, n)
                ref = maximum_flow(matrix, a1, a2)
                assert got == int(ref.flow_value), (got, ref.flow_value, edges)
                assert np.array_equal(flow, -flow.T)
                assert np.all(flow <= matrix.toarray())
                assert np.all(flow >= -matrix.toarray().T)
                net = flow.sum(axis=1)
                expected = np.zeros(n, dtype=np.int64)
                expected[a1], expected[a2] = got, -got
                assert np.array_equal(net, expected), (net, expected, edges)
                if n <= 9:
                    assert got == solve_bruteforce_cut(n, edges, a1, a2)
                    counts['brute_flow'] += 1
            elif task == 'M':
                one = list(map(int, next(output).split()))
                two = list(map(int, next(output).split()))
                status_rank, rank = map(int, next(output).split())
                assert one[0] == two[0] == status_rank == 0
                col_to_row = np.array(one[1:]) - 1
                row_to_col = np.array(two[1:]) - 1
                assert len(col_to_row) == a1 and len(row_to_col) == n
                assert np.count_nonzero(col_to_row >= 0) == rank
                assert np.count_nonzero(row_to_col >= 0) == rank
                assert rank == structural_rank(matrix)
                # Structural equality of the two inverse maps, including empty shapes.
                for c, r in enumerate(col_to_row):
                    if r >= 0:
                        assert 0 <= r < n and row_to_col[r] == c
                        assert (r, c) in {(u, v) for u, v, _ in edges}
                for r, c in enumerate(row_to_col):
                    if c >= 0:
                        assert 0 <= c < a1 and col_to_row[c] == r
                # The matching itself can differ; its size must be SciPy-optimal.
                scipy_match = maximum_bipartite_matching(matrix, perm_type='row')
                assert rank == np.count_nonzero(scipy_match >= 0)
                if a2 == 1:
                    assert rank == solve_bruteforce_matching(n,a1,[(r,c) for r,c,_ in edges])
                    counts['brute_match'] += 1
            else:
                line = list(map(int, next(output).split()))
                assert line[0] == 0, (task, line)
                perm = np.array(line[1:]) - 1
                assert sorted(perm.tolist()) == list(range(n))
                if n:
                    reference = reverse_cuthill_mckee(matrix, symmetric_mode=bool(a1))
                    assert sorted(reference.tolist()) == list(range(n))
                # The pseudo-peripheral selection and tie breaks are not part of
                # SciPy's public API; both produce valid RCM permutations.
                if n > 0 and len(edges):
                    # A loose quality guard on our ordering, not equality to SciPy.
                    ours = bandwidth(n,edges,perm)
                    original = bandwidth(n,edges,np.arange(n))
                    assert ours <= max(original, 2 * bandwidth(n,edges,reference) + 1), (ours, original)
            counts[task] += 1
        assert next(output,None) is None
        print(f'{mode}: '+' '.join(f'{k}={v}' for k,v in counts.items()), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    p.add_argument('--seed', type=int, default=41924)
    opts = p.parse_args()
    for mode in ('checked','optimized') if opts.mode=='both' else (opts.mode,):
        run(mode, seed=opts.seed)


if __name__ == '__main__':
    main()
