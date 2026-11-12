#!/usr/bin/env python3
"""Yen K-shortest simple paths: SciPy 1.17 and exhaustive DFS comparisons.

Run with --mode checked / optimized / both. Uses only stdlib, NumPy, SciPy and
GNU Fortran; no extra SciFort library or preexisting FPM build is required.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import random
import subprocess
import tempfile

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import NegativeCycleError, yen as scipy_yen

ROOT = Path(__file__).resolve().parents[2]


def create_cases(seed: int, count: int):
    rng = random.Random(seed)
    cases = []
    for iteration in range(count):
        n = rng.randrange(1, 10)
        source, target = rng.randrange(n), rng.randrange(n)
        k = rng.randrange(0, 15)
        directed = rng.choice([True, False])
        unweighted = rng.random() < 0.2
        family = iteration % 6
        if family in (1, 2, 5):
            directed = True
        if family == 4:
            n = max(n,3)
            source, target = rng.randrange(n), rng.randrange(n)
        edges = []
        # Cycles only have nonnegative weights in these two families.
        if family in (0, 3):
            candidates = [(u,v) for u in range(n) for v in range(n) if rng.random() < 0.36]
            weights = [float(rng.randint(0, 16))/4 for _ in candidates]
        elif family == 1:
            # DAG with signed weights, no negative cycles.
            candidates = [(u,v) for u in range(n) for v in range(u+1,n)
                          if rng.random() < 0.65]
            weights = [float(rng.randrange(-20,24))/4 for _ in candidates]
        elif family in (2,5):
            # Potential differences guarantee no negative cycles, despite
            # possibly negative arcs and strongly connected directed graphs.
            h = [rng.randrange(-12,13)/4 for _ in range(n)]
            candidates = [(u,v) for u in range(n) for v in range(n)
                          if rng.random() < 0.46]
            weights = [float(rng.randrange(0,14))/4 + h[v]-h[u] for u,v in candidates]
        else:
            # A known negative cycle, including when unreachable from source.
            candidates = [(0,1),(1,2),(2,0)]
            weights = [-3.0,1.0,1.0]
            candidates += [(u,v) for u in range(3,n) for v in range(3,n)
                           if u != v and rng.random()<0.18]
            weights += [float(rng.randrange(0,16))/4 for _ in candidates[3:]]
            unweighted = False
            k = max(1,k)
        edges.extend((u,v,w) for (u,v),w in zip(candidates,weights))
        if edges and rng.random() < 0.4:
            # Extra parallel edges of different weights, including zeros.
            for _ in range(rng.randrange(1,4)):
                u,v,w = rng.choice(edges)
                if family == 4:
                    alt = w+float(rng.randrange(0,8))/4
                elif family in (1,2,5):
                    alt = w+float(rng.randrange(0,8))/4
                else:
                    alt = float(rng.randint(0,16))/4
                edges.append((u,v,alt))
        rng.shuffle(edges)
        cases.append((n,source,target,k,directed,unweighted,edges))
    return cases


def edge_lookup(n, edges, directed, unweighted):
    """Store only the least expensive arc for a vertex-to-vertex transition."""
    lookup = {}
    for u,v,w in edges:
        w = 1.0 if unweighted else w
        if (u,v) not in lookup or w < lookup[(u,v)]:
            lookup[(u,v)] = w
        if not directed and ((v,u) not in lookup or w < lookup[(v,u)]):
            lookup[(v,u)] = w
    return lookup


def enumerate_simple_paths(n, source, target, lookup):
    """Independent exhaustive oracle, valid for negative arcs without cycles."""
    outgoing = defaultdict(list)
    for (u,v),w in lookup.items():
        outgoing[u].append((v,w))
    paths = []
    def visit(node, path, visited, cost):
        if node == target:
            paths.append((cost,tuple(path)))
            return
        for nxt,weight in outgoing[node]:
            if nxt in visited:
                continue
            visit(nxt,path+[nxt],visited | {nxt},cost+weight)
    visit(source,[source],{source},0.0)
    paths.sort(key=lambda e:(e[0],e[1]))
    return paths


def scipy_reference(n, source, target, k, directed, unweighted, edges):
    # Coalesce parallel edges by min rather than scipy.sparse COO's default SUM.
    lookup = edge_lookup(n,edges,True,False)
    rows = [a for a,b in lookup]
    cols = [b for a,b in lookup]
    values = [lookup[(a,b)] for a,b in lookup]
    graph = csr_array((np.asarray(values,dtype=float),(np.asarray(rows,dtype=np.int32),
                   np.asarray(cols,dtype=np.int32))),shape=(n,n))
    graph.indices = graph.indices.astype(np.int32)
    graph.indptr = graph.indptr.astype(np.int32)
    # SciPy 1.17 dereferences an invalid zero-length workspace in this case;
    # do not invoke the buggy routine at all.
    if k == 0 or (graph.nnz == 0 and source == target):
        return None,None,None
    try:
        scipy_dist, scipy_pred = scipy_yen(graph,source,target,k,
            directed=directed,unweighted=unweighted,return_predecessors=True)
    except NegativeCycleError:
        return 2,None,None
    except ValueError as error:
        # SciPy 1.17 attempts a zero-sized memoryview for empty CSR and
        # source == sink; the independent exhaustive check still applies.
        if graph.nnz == 0 and source == target and "Invalid shape" in str(error):
            return None,None,None
        raise
    return 0,scipy_dist,scipy_pred


def decode_path(row, source, target):
    if source == target:
        assert not any(row), row
        return (source,)
    route = [target]
    seen = {target}
    for _ in range(len(row)):
        u = int(row[route[-1]])-1
        if u < 0:
            raise AssertionError(f'broken predecessor chain: {row}')
        if u in seen:
            raise AssertionError(f'cyclic predecessor chain: {row}')
        route.append(u)
        seen.add(u)
        if u == source:
            return tuple(reversed(route))
    raise AssertionError(f'unbounded predecessor chain: {row}')


def compare(exe: Path, cases):
    lines = [str(len(cases))]
    for n,source,target,k,directed,unweighted,edges in cases:
        lines.append(f'{n} {len(edges)} {source+1} {target+1} {k} {int(directed)} {int(unweighted)}')
        lines.extend(f'{u+1} {v+1} {w:.17g}' for u,v,w in edges)
    result = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',text=True,
                            capture_output=True,check=True)
    out = iter(result.stdout.splitlines())
    counts = dict(cases=0, paths=0, exhaustive=0, negative_cycles=0, zero_k=0, scipy_skipped=0)
    for case in cases:
        n,source,target,k,directed,unweighted,edges = case
        status, npaths = map(int,next(out).split())
        got_costs = []
        got_paths = []
        for _ in range(npaths):
            row = next(out).split()
            assert len(row)==n+1, (row,n)
            got_costs.append(float(row[0]))
            got_paths.append(decode_path([int(x) for x in row[1:]],source,target))
        ref_status,ref_costs,ref_pred = scipy_reference(*case)
        if ref_status is not None:
            assert status==ref_status,(case,status,ref_status,got_costs,ref_costs)
        else:
            counts["scipy_skipped"] += 1
        if status == 2:
            assert npaths==0,(case,npaths)
            counts['negative_cycles'] += 1
            continue
        assert npaths <= k
        if not k:
            counts['zero_k'] += 1
        lookup = edge_lookup(n,edges,directed,unweighted)
        assert len(set(got_paths)) == npaths,(case,got_paths)
        for nodes,cost in zip(got_paths,got_costs):
            assert nodes[0]==source and nodes[-1]==target
            assert len(set(nodes))==len(nodes)
            expected_cost = sum(lookup[(a,b)] for a,b in zip(nodes,nodes[1:]))
            assert np.isclose(expected_cost,cost,rtol=1e-11,atol=1e-11), (nodes,cost,expected_cost)
        assert np.all(np.diff(got_costs)>=-1e-10),(case,got_costs)
        # Path tie-breaking can differ, but the sorted optimal costs must agree.
        if ref_status is not None:
            assert len(got_costs)==len(ref_costs),(case,got_costs,ref_costs)
            assert np.allclose(got_costs,ref_costs,rtol=1e-10,atol=1e-10), (case,got_costs,ref_costs)
        # Independently enumerate all routes on small test graphs.
        if n <= 8:
            exact = enumerate_simple_paths(n,source,target,lookup)
            target_costs = [cost for cost,path in exact[:k]]
            assert np.allclose(got_costs,target_costs,rtol=1e-11,atol=1e-11), (case,got_costs,target_costs,got_paths)
            counts['exhaustive'] += 1
        counts['paths'] += npaths
        counts['cases'] += 1
    assert next(out,None) is None
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=44187)
    parser.add_argument('--cases',type=int,default=240)
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    args = parser.parse_args()
    cases = create_cases(args.seed,args.cases)
    for mode in ('checked','optimized') if args.mode=='both' else (args.mode,):
        flags = (['-O0','-g','-fcheck=all','-fbacktrace'] if mode=='checked' else ['-O3'])
        with tempfile.TemporaryDirectory(prefix='scifort-yen-') as temp:
            build = Path(temp)
            cc = ['gfortran','-std=f2018','-ffree-line-length-none',*flags,
                '-I',str(build),'-J',str(build)]
            for name in ('kinds','sparse_csgraph'):
                subprocess.run([*cc,'-c',str(ROOT/'src'/f'{name}.f90'),
                    '-o',str(build/f'{name}.o')],cwd=build,check=True)
            exe = build/'csgraph_yen_driver'
            subprocess.run([*cc,str(ROOT/'tools/csgraph_reference/csgraph_yen_driver.f90'),
                str(build/'kinds.o'),str(build/'sparse_csgraph.o'),'-o',str(exe)],
                cwd=build,check=True)
            print(f'{mode}: {compare(exe,cases)}',flush=True)


if __name__=='__main__':
    main()
