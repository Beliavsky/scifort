#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exploratory comparison of five bootstrap interval methods under VAR(1).

Uses the same complete 12-sample tapered Fourier windows for five methods,
fixed-length circular blocks and identical deterministic bootstrap draws.
The finite-window spectral estimator targets a band-smoothed quantity, so
coverage against the population zero-frequency coherence is exploratory.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from study_canonical_influence_var import (make_var_matrices,
                                         true_zero_frequency_coherence)


def run(driver: str, trials: int) -> None:
    rng = np.random.default_rng(20261030)
    nw, nperseg = 48, 12
    sigma = np.array([[1.,.5,0.],[.5,1.,0.],[0.,0.,900.]])
    scenarios = (
        ('persistent target',np.array([[.95,.035,0.],[0.,.94,0.],[0.,0.,.2]])),
        ('moderate target',np.array([[.60,.08,0.],[0.,.55,0.],[0.,0.,.96]])),
    )
    lines = [str(trials*len(scenarios))]
    for _,a in scenarios:
        for trial in range(trials):
            samples = make_var_matrices(rng,nw,nperseg,a,sigma)
            lines.append(f'{nw} 3 {1000+trial}')
            for matrix in samples:
                values = matrix.flatten(order='F')
                pairs = np.column_stack((values.real,values.imag)).ravel()
                lines.append(' '.join(f'{x:.17g}' for x in pairs))
    process = subprocess.run([driver],input='\n'.join(lines)+'\n',text=True,
                             capture_output=True,check=True)
    observations = [row.split() for row in process.stdout.splitlines()]
    assert len(observations) == trials*len(scenarios)
    methods = ['percentile','basic','bias-corrected','BCa','studentized']
    for i,(name,a) in enumerate(scenarios):
        truth = true_zero_frequency_coherence(a,sigma)
        rows = np.asarray([[float(x) for x in row[:10]]
                           for row in observations[i*trials:(i+1)*trials]])
        flags = np.asarray([[int(x) for x in row[-3:]] for row in
                            observations[i*trials:(i+1)*trials]],dtype=bool)
        print(f'VAR(1) {name}: true zero-frequency squared coherence={truth:.5f}; '
              f'{trials} datasets, 199 replicates, circular block length 4')
        for j,method in enumerate(methods):
            valid = flags[:,j-2] if j >= 2 else np.ones(trials,dtype=bool)
            if not valid.any():
                print(f'  {method}: 0 valid intervals')
                continue
            bounds = rows[valid,2*j:2*j+2]
            coverage = np.mean((bounds[:,0] <= truth)&(truth <= bounds[:,1]))
            width = np.mean(bounds[:,1]-bounds[:,0])
            print(f'  {method}: coverage={coverage:.3f}, '
                  f'mean width={width:.4f}, valid={np.sum(valid)}/{trials}')
    print('Exploratory only: no empirical-coverage claim or method ranking is implied.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    parser.add_argument('--trials',type=int,default=16)
    args = parser.parse_args()
    run(args.driver,args.trials)
