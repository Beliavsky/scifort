#!/usr/bin/env python3
"""Randomized SciPy and brute-force checks of sparse assignment and CSR Laplacians."""
from __future__ import annotations

import argparse
from itertools import permutations
from pathlib import Path
import random
import subprocess
import tempfile

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import laplacian as scipy_laplacian
from scipy.sparse.csgraph import min_weight_full_bipartite_matching

ROOT = Path(__file__).resolve().parents[2]


def cases_for_seed(seed: int, count: int):
    rng = random.Random(seed)
    cases = []
    for _ in range(count):
        if rng.random() < 0.08:
            nrows, ncols = rng.randrange(0, 12), rng.randrange(0, 12)
        else:
            nrows, ncols = rng.randrange(1, 10), rng.randrange(1, 10)
        prob = rng.uniform(0.16, 1.0)
        edges = [(i, j, float(rng.randint(-20, 20)) + (0.125 if rng.random() < 0.4 else 0.0))
                 for i in range(nrows) for j in range(ncols) if rng.random() < prob]
        if rng.random() < 0.4:
            # Force a full matching, while retaining the original weights.
            for i, j in zip(rng.sample(range(nrows), min(nrows, ncols)),
                            rng.sample(range(ncols), min(nrows, ncols))):
                if not any(r == i and c == j for r, c, _ in edges):
                    edges.append((i, j, float(rng.randint(-20, 20))))
        if edges and rng.random() < 0.3:
            i,j,_ = rng.choice(edges)
            edges.append((i,j,float(rng.randint(-20,20)) + 0.25))
        rng.shuffle(edges)
        cases.append(('W', nrows, ncols, rng.randrange(2), 0, 0, edges))

    for _ in range(count):
        n = rng.randrange(0, 25)
        num_edges = rng.randrange(0, 4*n+1)
        # Mixed signed edge weights; normalization rejects negative degrees.
        positive = rng.random() < 0.7
        edges = [(rng.randrange(n), rng.randrange(n), float(rng.randint(0 if positive else -5, 12)))
                 for _ in range(num_edges)] if n else []
        cases.append(('L', n, n, rng.randrange(2), rng.randrange(2), rng.randrange(2), edges))
    return cases


def oracle_small(nr, nc, edges, maximize):
    """Exhaustive matching optimum, independent of SciPy's matching algorithm."""
    lookup = {}
    for r, c, w in edges:
        lookup[(r, c)] = min(lookup.get((r, c), w), w) if not maximize else max(lookup.get((r, c), w), w)
    k = min(nr, nc)
    if k == 0:
        return 0.0
    if nr <= nc:
        candidate_sets = ([(i, p[i]) for i in range(k)] for p in permutations(range(nc), k))
    else:
        candidate_sets = ([(p[i], i) for i in range(k)] for p in permutations(range(nr), k))
    values = [sum(lookup[x] for x in matching) for matching in candidate_sets
              if all(x in lookup for x in matching)]
    if not values:
        return None
    return max(values) if maximize else min(values)


def check_cases(exe: Path, cases):
    inp = [str(len(cases))]
    for task, nr, nc, a, b, c, edges in cases:
        inp.append(f'{task} {nr} {nc} {len(edges)} {a}')
        inp.extend(f'{r+1} {col+1} {w:.17g}' for r, col, w in edges)
        if task == 'L':
            inp.append(f'{b} {c}')
    result = subprocess.run([str(exe)], input='\n'.join(inp)+'\n', capture_output=True,
                            text=True, check=True)
    output = iter(result.stdout.splitlines())
    counts = {'weighted': 0, 'infeasible': 0, 'exhaustive': 0,
              'laplacian': 0, 'negative_degree': 0}
    for task, nr, nc, a, b, c, edges in cases:
        if task == 'W':
            status_str, cost_str, n_str = next(output).split()
            status, got_cost, got_size = int(status_str), float(cost_str), int(n_str)
            # SciPy's CSR constructor sums duplicates; collapse to the cheapest
            # (or dearest) parallel edge before passing weights to its solver.
            lookup = {}
            for r, j, w in edges:
                key = (r,j)
                if key not in lookup or ((w > lookup[key]) if a else (w < lookup[key])):
                    lookup[key] = w
            if nr and nc and lookup:
                shifted = np.asarray([w+100.0 for w in lookup.values()], dtype=float)
                matrix = csr_matrix((shifted, ([r for r, _ in lookup],
                                               [j for _, j in lookup])), shape=(nr, nc))
            else:
                matrix = csr_matrix((nr, nc), dtype=float)
            scipy_possible = True
            try:
                ref_r, ref_c = min_weight_full_bipartite_matching(matrix, maximize=bool(a))
                reference_cost = float(matrix[ref_r, ref_c].sum()) - 100.0*len(ref_r)
            except ValueError:
                scipy_possible = False
                reference_cost = None
            if not scipy_possible:
                assert status == 3, (nr, nc, edges, status)
                counts['infeasible'] += 1
            else:
                assert status == 0, (nr, nc, a, edges, status)
                assert got_size == min(nr, nc)
                pairs = [tuple(map(int, next(output).split())) for _ in range(got_size)]
                assert [r for r, _ in pairs] == sorted(r for r, _ in pairs)
                assert len({r for r, _ in pairs}) == len({col for _, col in pairs}) == got_size
                expected_cost = 0
                for pair in pairs:
                    key = (pair[0]-1, pair[1]-1)
                    assert key in lookup, (pair, edges)
                    expected_cost += lookup[key]
                assert np.isclose(got_cost, expected_cost, atol=2e-9), (cost_str, expected_cost)
                assert np.isclose(got_cost, reference_cost, atol=2e-9), (got_cost, reference_cost, pairs, edges)
                counts['weighted'] += 1
            if nr <= 5 and nc <= 5:
                exact = oracle_small(nr,nc,edges,bool(a))
                assert (exact is None) == (status==3), (nr,nc,edges,exact,status)
                if exact is not None:
                    assert np.isclose(got_cost, exact, atol=2e-9), (got_cost,exact,edges)
                counts['exhaustive'] += 1
        else:
            status, n = map(int, next(output).split())
            assert n == nr
            rows = np.array([e[0] for e in edges], dtype=int)
            cols = np.array([e[1] for e in edges], dtype=int)
            weights = np.array([e[2] for e in edges], dtype=float)
            matrix = coo_matrix((weights, (rows, cols)), shape=(nr, nr)).tocsr()
            # SciPy may produce nonfinite normalized entries on signed degree.
            if a:
                _, unnormalized_degrees = scipy_laplacian(matrix,normed=False,return_diag=True,
                    use_out_degree=bool(b),symmetrized=bool(c))
                if np.any(unnormalized_degrees < 0):
                    assert status == 1, (edges, a,b,c,status)
                    counts['negative_degree'] += 1
                    continue
            assert status == 0, (status, a,b,c,edges)
            degrees = np.array([float(next(output)) for _ in range(nr)])
            dense = np.array([[float(x) for x in next(output).split()] for _ in range(nr)]).reshape(nr,nr)
            ref, ref_degree = scipy_laplacian(matrix, normed=bool(a), return_diag=True,
                                              use_out_degree=bool(b),symmetrized=bool(c))
            assert np.allclose(degrees,ref_degree,rtol=1e-10,atol=1e-10), (degrees, ref_degree, edges)
            assert np.allclose(dense,ref.toarray(),rtol=1e-10,atol=1e-10), (dense,ref.toarray(),edges)
            counts['laplacian'] += 1
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=33047)
    parser.add_argument('--cases', type=int, default=220)
    parser.add_argument('--mode', choices=('checked','optimized','both'), default='both')
    args = parser.parse_args()
    cases = cases_for_seed(args.seed,args.cases)
    for mode in ('checked','optimized') if args.mode=='both' else (args.mode,):
        flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3']
        with tempfile.TemporaryDirectory(prefix='scifort-weighted-') as d:
            path = Path(d)
            cc = ['gfortran','-std=f2018','-ffree-line-length-none',*flags,'-I',str(path),'-J',str(path)]
            for name in ('kinds','sparse_csgraph'):
                subprocess.run([*cc,'-c',str(ROOT/'src'/f'{name}.f90'),'-o',str(path/f'{name}.o')],
                               cwd=path, check=True)
            exe = path/'weighted_driver'
            subprocess.run([*cc,str(ROOT/'tools/csgraph_reference/csgraph_weighted_driver.f90'),
                            str(path/'kinds.o'),str(path/'sparse_csgraph.o'),'-o',str(exe)], cwd=path,check=True)
            print(f'{mode}: {check_cases(exe,cases)}',flush=True)


if __name__ == '__main__':
    main()
