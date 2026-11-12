#!/usr/bin/env python3
"""Independent NumPy/SciPy numerical checks for the SciFort vq module.

Builds the standalone Fortran executable in a temporary directory; fpm is
not required. Uses an explicitly supplied centroid matrix so results do not
rely on NumPy and SciFort having identical random-number generators.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile
import warnings

import numpy as np
from scipy.cluster import vq as scipy_vq

ROOT = Path(__file__).resolve().parents[2]


def build(folder: Path, mode: str) -> Path:
    flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if mode == 'checked' else ['-O3']
    command = ['gfortran', '-std=f2018', *flags, '-J', str(folder), '-I', str(folder)]
    objects = []
    for name in ('kinds', 'random', 'cluster_vq'):
        obj = folder / f'{name}.o'
        subprocess.run([*command, '-c', str(ROOT/'src'/f'{name}.f90'),
                        '-o', str(obj)], cwd=folder, check=True, capture_output=True)
        objects.append(str(obj))
    exe = folder/'cluster_vq_driver'
    subprocess.run([*command, str(ROOT/'tools/cluster_reference/cluster_vq_driver.f90'),
                    *objects, '-o', str(exe)], cwd=folder, check=True, capture_output=True)
    return exe


def execute(exe: Path, operation: str, obs: np.ndarray, guess: np.ndarray, steps: int = 12):
    n, p = obs.shape
    k = len(guess)
    lines = [f'{operation} {n} {p} {k} {steps}']
    lines.extend(' '.join(f'{v:.17e}' for v in row) for row in obs)
    if operation != 'whiten':
        lines.extend(' '.join(f'{v:.17e}' for v in row) for row in guess)
    result = subprocess.run([str(exe)], input='\n'.join(lines)+'\n',
                            text=True, capture_output=True, check=True)
    parts = result.stdout.strip().splitlines()
    status = int(parts[0])
    if status != 0:
        raise AssertionError(f'Fortran {operation} status={status}, stderr={result.stderr}')
    if operation=='whiten':
        return np.array([[float(x) for x in row.split()] for row in parts[1:]])
    if operation=='vq':
        return np.fromstring(parts[1],sep=' ',dtype=int), np.fromstring(parts[2],sep=' ')
    if operation=='kmeans2':
        return (np.array([[float(x) for x in row.split()] for row in parts[1:k+1]]),
                np.fromstring(parts[k+1],sep=' ',dtype=int))
    if operation=='kmeans':
        nkeep=int(parts[1])
        return (np.array([[float(x) for x in row.split()] for row in parts[2:2+nkeep]]),
                float(parts[2+nkeep]))
    raise AssertionError(f'Unknown operation {operation}')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--cases',type=int,default=100)
    parser.add_argument('--seed',type=int,default=20261103)
    parser.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    args=parser.parse_args()
    modes=['checked','optimized'] if args.mode=='both' else [args.mode]
    for mode in modes:
        rng=np.random.default_rng(args.seed)
        with tempfile.TemporaryDirectory(prefix='scifort-vq-') as dirname:
            exe=build(Path(dirname),mode)
            for idx in range(args.cases):
                n=int(rng.integers(12,44))
                p=int(rng.integers(1,7))
                k=int(rng.integers(2,5))
                means=rng.uniform(-8,8,size=(k,p))
                means+=np.arange(k)[:,None]*7.0
                which=rng.integers(0,k,size=n)
                # All groups represented; well-separated means avoid unstable local optima.
                which[:k]=np.arange(k)
                obs=means[which]+rng.normal(0,0.45,size=(n,p))
                obs*=rng.uniform(0.1,5.0,size=p)
                initial=np.array([obs[which==j].mean(axis=0) for j in range(k)])
                initial+=rng.normal(size=(k,p))*0.01
                # vq and whiten, including a constant feature on selected cases.
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',RuntimeWarning)
                    expect_whiten=scipy_vq.whiten(obs)
                got_whiten=execute(exe,'whiten',obs,initial)
                np.testing.assert_allclose(got_whiten,expect_whiten,atol=3e-12,rtol=3e-12)
                expect_labels,expect_dist=scipy_vq.vq(obs,initial)
                got_labels,got_dist=execute(exe,'vq',obs,initial)
                np.testing.assert_array_equal(got_labels,expect_labels+1)
                np.testing.assert_allclose(got_dist,expect_dist,atol=2e-11,rtol=2e-11)
                expect_centers,expect_labels=scipy_vq.kmeans2(
                    obs,initial,iter=14,minit='matrix',missing='raise')
                got_centers,got_labels=execute(exe,'kmeans2',obs,initial,steps=14)
                np.testing.assert_allclose(got_centers,expect_centers,atol=5e-11,rtol=5e-11)
                np.testing.assert_array_equal(got_labels,expect_labels+1)
                expect_centers,expect_dist=scipy_vq.kmeans(obs,initial,thresh=1e-10)
                got_centers,got_dist=execute(exe,'kmeans',obs,initial)
                np.testing.assert_allclose(got_centers,expect_centers,atol=5e-11,rtol=5e-11)
                np.testing.assert_allclose(got_dist,expect_dist,atol=5e-11,rtol=5e-11)
        print(f'{mode}: {args.cases} SciPy whiten/vq/kmeans2/kmeans comparisons PASS', flush=True)


if __name__=='__main__':
    main()
