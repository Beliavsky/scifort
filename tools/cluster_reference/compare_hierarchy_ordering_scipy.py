#!/usr/bin/env python3
"""Compare SciFort hierarchy leaf-order optimization with SciPy and brute force.

SciPy's orientation is not treated as the one true sequence: equal optima
may differ, and the installed SciPy implementation sometimes returns a
slightly longer path. We independently enumerate all allowed tree flips
for small cases, and require SciFort's path cost to reach that optimum.
"""
from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from numpy.testing import assert_allclose
from scipy.cluster import hierarchy as sch
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parents[2]


def compile_native(workdir: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    command = ['gfortran','-std=f2018','-ffree-line-length-none',*flags,
               '-J',str(workdir),'-I',str(workdir)]
    objects = []
    for module in ('kinds','spatial_distance','cluster_hierarchy'):
        obj = workdir/(module+'.o')
        subprocess.run([*command, '-c', str(ROOT/'src'/(module+'.f90')),
                        '-o',str(obj)],cwd=workdir,check=True)
        objects.append(str(obj))
    exe = workdir/'ordering_driver'
    subprocess.run([*command,str(ROOT/'tools/cluster_reference'/'cluster_hierarchy_ordering_driver.f90'),
                    *objects,'-o',str(exe)],cwd=workdir,check=True)
    return exe


def native_optimal(exe: Path, z: np.ndarray, y: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    n = len(z)+1
    data = [str(n)]
    data.extend(f'{v:.17e}' for v in y)
    for row in z:
        data.append(' '.join(f'{v:.17e}' for v in row))
    result = subprocess.run([str(exe)],input='\n'.join(data)+'\n',
                            capture_output=True,text=True,check=True)
    lines = result.stdout.splitlines()
    status, count = map(int,lines[0].split())
    if status != 0 or count != n-1:
        raise AssertionError(f'Native invalid output: {result.stdout} {result.stderr}')
    output = np.array([[float(v) for v in ln.split()] for ln in lines[1:n]])
    leaves = np.array([int(v) for v in lines[n].split()],dtype=int)
    return output,leaves


def tree_valid(original: np.ndarray, ordered: np.ndarray, leaves: np.ndarray):
    np.testing.assert_array_equal(np.sort(ordered[:,:2],axis=1),
                                  np.sort(original[:,:2],axis=1))
    np.testing.assert_array_equal(ordered[:,2:4],original[:,2:4])
    n = len(leaves)
    np.testing.assert_array_equal(np.sort(leaves),np.arange(1,n+1))
    check = ordered.copy()
    check[:,:2] -= 1
    if not sch.is_valid_linkage(check):
        raise AssertionError('native tree invalid according to SciPy')
    np.testing.assert_array_equal(sch.leaves_list(check)+1,leaves)


def cost(order: np.ndarray, dissimilarity: np.ndarray) -> float:
    o = np.asarray(order,dtype=int)-1
    return float(np.sum(dissimilarity[o[:-1],o[1:]]))


def all_allowed_orders(z: np.ndarray):
    n = len(z)+1
    def descend(node):
        if node < n: return [(node+1,)]
        a,b = map(int,z[node-n,:2])
        left = descend(a-1)
        right = descend(b-1)
        return [a_order+b_order for a_order,b_order in product(left,right)] + [
                b_order+a_order for a_order,b_order in product(left,right)]
    return descend(2*n-2)


def test_one(exe: Path, y: np.ndarray, z: np.ndarray,
             exhaustive: bool) -> tuple[bool,bool]:
    n = len(z)+1
    input_z = z.copy()
    input_z[:,:2] += 1
    ordered, leaves = native_optimal(exe,input_z,y)
    tree_valid(input_z,ordered,leaves)
    D = squareform(y)
    native_score = cost(leaves,D)
    base_score = cost(sch.leaves_list(z)+1,D)
    if native_score > base_score+2e-11*max(1.0,base_score):
        raise AssertionError(f'Reordering worsened original: {native_score} > {base_score}')
    scipy_score = cost(sch.leaves_list(sch.optimal_leaf_ordering(z,y))+1,D)
    if native_score > scipy_score+2e-11*max(1.0,scipy_score):
        raise AssertionError(f'SciFort score {native_score:.16g} exceeds SciPy {scipy_score:.16g}')
    if exhaustive:
        exhaustive_best = min(cost(o,D) for o in all_allowed_orders(input_z))
        assert_allclose(native_score,exhaustive_best,rtol=2e-12,atol=2e-12,
                        err_msg='not the global optimal ordering over all tree flips')
    return (native_score+1e-10 < scipy_score, native_score+1e-10 < base_score)


def run_one(mode: str, cases: int, seed: int):
    rng = np.random.default_rng(seed)
    improved_over_scipy = 0
    improved_over_original = 0
    exhaustive_count = 0
    methods = ('single','complete','average','weighted','centroid','median','ward')
    with tempfile.TemporaryDirectory(prefix='scifort-ordering-') as tmp:
        exe = compile_native(Path(tmp),mode)
        for c in range(cases):
            n = int(rng.integers(2,13 if c%4==0 else 25))
            p = int(rng.integers(1,6))
            method = methods[c%len(methods)]
            X = rng.standard_normal((n,p))
            if c%11 == 0:
                X = np.round(X,1)
            if c%13 == 0:
                X = X*1000.0
            if c%17 == 0:
                X[0] = X[1]
            y = pdist(X)
            z = sch.linkage(y,method=method)
            exhaustive = n<=8
            better_scipy,better_original = test_one(exe,y,z,exhaustive)
            improved_over_scipy += better_scipy
            improved_over_original += better_original
            exhaustive_count += exhaustive
        # Zero dissimilarities, ties, singletons and non-Euclidean distances.
        for n in (2,3,4,6,8):
            y = np.zeros(n*(n-1)//2)
            z = sch.linkage(y,method='complete')
            test_one(exe,y,z,True)
        for n in (5,8,11):
            y = rng.uniform(0.1,5.0,size=n*(n-1)//2)
            z = sch.linkage(y,method='centroid')
            test_one(exe,y,z,n<=8)
    print(f'{mode}: {cases} random trees + 8 edge cases PASS; '
          f'{exhaustive_count} random cases independently brute-force verified; '
          f'{improved_over_original} strictly shorter than original; '
          f'{improved_over_scipy} shorter than SciPy 1.17',flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=200)
    parser.add_argument('--seed',type=int,default=26909)
    args = parser.parse_args()
    for mode in (('checked','optimized') if args.mode=='both' else (args.mode,)):
        run_one(mode,args.cases,args.seed)

if __name__ == '__main__':
    main()
