#!/usr/bin/env python3
"""Randomized SciPy tests for shortest_path dispatcher and graph trees."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import (
    NegativeCycleError, breadth_first_tree, depth_first_tree, shortest_path,
)
from compare_csgraph_extended_scipy import cases

ROOT = Path(__file__).resolve().parents[2]


def compare(mode: str) -> None:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    graphs = cases()[::3]  # 70 deterministic graphs, many with negative edges
    with tempfile.TemporaryDirectory(prefix='scifort-csgraph-dispatch-') as temp:
        build = Path(temp)
        cc = ['gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
              '-J', str(build), '-I', str(build)]
        objects = []
        for name in ('kinds', 'sparse_csgraph'):
            obj = build/f'{name}.o'
            subprocess.run([*cc, '-c', str(ROOT/'src'/f'{name}.f90'), '-o', str(obj)],
                           check=True, cwd=build)
            objects.append(str(obj))
        exe = build/'csgraph_dispatch'
        subprocess.run([*cc, str(Path(__file__).with_name('csgraph_dispatch_driver.f90')),
                        *objects, '-o', str(exe)], cwd=build, check=True)
        inp = [str(len(graphs))]
        for a in graphs:
            inp.append(str(len(a)))
            inp += [' '.join(map(str,row)) for row in a]
        result = subprocess.run([str(exe)], input='\n'.join(inp)+'\n',
                                text=True, capture_output=True, cwd=build, check=True)
        output = result.stdout.splitlines()
        cursor = checked = cycles = trees = 0
        for number, a in enumerate(graphs):
            n = int(output[cursor]); cursor+=1
            assert n == len(a)
            graph = csr_array(a)
            sources = np.array([0,n-1])
            for directed in (False,True):
                for unweighted in (False,True):
                    for algorithm in ('auto','FW','J','BF','D'):
                        status = int(output[cursor]); cursor += 1
                        if algorithm=='D' and not unweighted and np.any(a < 0):
                            assert status == 1  # SciFort rejects negative weights for Dijkstra
                            continue
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore', UserWarning)
                            try:
                                # SciPy cannot take indices with FW; slice full APSP.
                                if algorithm == 'FW':
                                    ref = shortest_path(graph,method=algorithm,directed=directed,
                                                        unweighted=unweighted)[sources]
                                else:
                                    ref = shortest_path(graph,method=algorithm,directed=directed,
                                                        unweighted=unweighted,indices=sources)
                            except NegativeCycleError:
                                assert status == 2, (mode,number,algorithm,directed,unweighted,status)
                                cycles += 1
                                continue
                        if status==2 and algorithm=='BF':
                            # BF is source-local: SciPy raises if any queried source
                            # reaches a negative cycle. This must agree for our chosen roots.
                            raise AssertionError((mode,number,algorithm,'false negative cycle'))
                        assert status==0, (mode,number,algorithm,directed,unweighted,status)
                        dist=np.array([np.fromstring(line,sep=' ') for line in output[cursor:cursor+2]])
                        cursor+=2
                        pred=np.array([np.fromstring(line,sep=' ',dtype=int) for line in output[cursor:cursor+2]])
                        cursor+=2
                        assert dist.shape==(2,n) and pred.shape==(2,n)
                        assert np.allclose(dist,ref,atol=1e-10,rtol=1e-10), (
                            mode,number,algorithm,directed,unweighted,dist,ref)
                        checked+=1
                for function in (breadth_first_tree,depth_first_tree):
                    status=int(output[cursor]);cursor+=1
                    assert status==0
                    tree=np.array([np.fromstring(line,sep=' ') for line in output[cursor:cursor+n]])
                    cursor+=n
                    reference=function(graph,0,directed=directed).toarray()
                    assert np.count_nonzero(tree)==np.count_nonzero(reference), (
                        mode,number,function.__name__,directed,tree,reference)
                    if directed or np.array_equal(a,a.T):
                        assert np.array_equal(tree,reference), (
                            mode,number,function.__name__,directed,tree,reference)
                    trees+=1
        assert cursor==len(output), (cursor,len(output))
    print(f'{mode}: {len(graphs)} graphs; {checked} shortest_path comparisons; '
          f'{cycles} negative cycles; {trees} search trees PASS',flush=True)


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['both','checked','optimized'],default='both')
    args=parser.parse_args()
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        compare(mode)


if __name__=='__main__':
    main()
