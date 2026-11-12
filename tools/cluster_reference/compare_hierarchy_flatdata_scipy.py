#!/usr/bin/env python3
"""Compare Fortran fclusterdata/fcluster and cophenetic correlation with SciPy.

Partitions are compared modulo label relabeling. SciPy is used only as an
independent numerical oracle, not as a source-code dependency.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

ROOT = Path(__file__).resolve().parents[2]
METHODS = ('single','complete','average','weighted','centroid','median','ward')
METRICS = ('euclidean','cityblock','chebyshev','sqeuclidean','cosine','minkowski')
CRITERIA = ('inconsistent','distance','maxclust')


def compile_native(workdir: Path, mode: str) -> Path:
    flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if mode == 'checked' else ['-O3']
    base = ['gfortran','-std=f2018','-ffree-line-length-none',*flags,
            '-J',str(workdir),'-I',str(workdir)]
    objects = []
    for module in ('kinds','spatial_distance','cluster_hierarchy'):
        obj = workdir / (module + '.o')
        subprocess.run([*base,'-c',str(ROOT/'src'/(module+'.f90')),
                        '-o',str(obj)],cwd=workdir,check=True)
        objects.append(str(obj))
    exe = workdir / 'flatdata_driver'
    subprocess.run([*base,str(ROOT/'tools/cluster_reference'/'cluster_hierarchy_flatdata_driver.f90'),
                    *objects,'-o',str(exe)],cwd=workdir,check=True)
    return exe


def native_run(exe: Path, x: np.ndarray, method: str, metric: str,
               criterion: str, threshold: float, depth: int):
    lines = [f'{len(x)} {x.shape[1]}',method,metric,criterion,
             f'{threshold:.17e} {depth}']
    lines += [' '.join(f'{v:.17e}' for v in row) for row in x]
    result = subprocess.run([str(exe)],input='\n'.join(lines)+'\n',text=True,
                            capture_output=True,check=True)
    output = result.stdout.strip().splitlines()
    status, size, corr_status, manual_size, corr = output[0].split()
    status, size, corr_status, manual_size = map(int,(status,size,corr_status,manual_size))
    if status != 0 or corr_status != 0 or size != len(x) or manual_size != len(x):
        raise AssertionError(f'Unexpected native status {output}, stderr={result.stderr}')
    labels = np.array([int(v) for v in output[1].split()],dtype=int)
    manual = np.array([int(v) for v in output[2].split()],dtype=int)
    return labels, manual, float(corr)


def same_partition(a: np.ndarray,b: np.ndarray) -> bool:
    return bool(np.array_equal(a[:,None]==a[None,:],b[:,None]==b[None,:]))


def run_one(exe: Path, rng: np.random.Generator, i: int):
    n = int(rng.integers(4,17))
    d = int(rng.integers(1,6))
    x = rng.normal(size=(n,d)) + rng.normal(size=d)*0.4
    method = METHODS[i%len(METHODS)]
    metric = METRICS[(i//len(METHODS))%len(METRICS)]
    criterion = CRITERIA[(i//(len(METHODS)*len(METRICS)))%len(CRITERIA)]
    depth = int(rng.integers(0,5))
    y = pdist(x,metric=metric)
    z = hierarchy.linkage(y,method=method)
    if criterion == 'maxclust':
        threshold = float(rng.integers(1,n+1))
    elif criterion == 'distance':
        height = float(rng.choice(z[:,2]))
        threshold = height*(1.0+float(rng.choice([-1,1]))*1.0e-7)
    else:
        r = hierarchy.inconsistent(z,d=depth)
        m = float(np.max(r[:,3]))
        threshold = float(rng.uniform(-0.3,m+0.4))
    flat = hierarchy.fclusterdata(x,threshold,criterion=criterion,
                                  metric=metric,depth=depth,method=method)
    direct = hierarchy.fcluster(z,threshold,criterion=criterion,depth=depth)
    if not same_partition(flat,direct):
        raise AssertionError(f'SciPy fcluster vs fclusterdata differ case {i}')
    native, manual, corr = native_run(exe,x,method,metric,criterion,threshold,depth)
    if not same_partition(native,flat):
        raise AssertionError(f'case {i}: fclusterdata {criterion} {method} {metric} n={n} d={d} '
                             f't={threshold} depth={depth}: native {native}; scipy {flat}')
    if not same_partition(manual,flat):
        raise AssertionError(f'case {i}: fcluster manual differs from SciPy')
    scipy_corr = float(hierarchy.cophenet(z,y)[0])
    if not np.isclose(corr,scipy_corr,rtol=3e-9,atol=3e-10,equal_nan=True):
        raise AssertionError(f'case {i}: cophenetic correlation {corr}, SciPy {scipy_corr}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('checked','optimized','both'),default='both')
    parser.add_argument('--cases',type=int,default=252)
    args = parser.parse_args()
    if args.cases < 1:
        parser.error('cases must be positive')
    modes = ('checked','optimized') if args.mode == 'both' else (args.mode,)
    for mode in modes:
        with tempfile.TemporaryDirectory(prefix='scifort-flatdata-') as scratch:
            exe = compile_native(Path(scratch),mode)
            rng = np.random.default_rng(20261103)
            for i in range(args.cases):
                run_one(exe,rng,i)
        print(f'{mode}: {args.cases} randomized fclusterdata and correlation cases PASS')

if __name__=='__main__':
    main()
