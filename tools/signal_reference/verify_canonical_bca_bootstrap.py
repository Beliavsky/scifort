#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent NumPy/SciPy block-aware BCa canonical bootstrap reference.

Regenerates all draws with an independent Park-Miller RNG, calculates all
canonical roots with NumPy, estimates block-jackknife accelerations with
NumPy and evaluates BCa endpoints using SciPy standard-normal functions.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.stats import norm

from verify_canonical_jackknife import canonical, integrated, regular
from verify_canonical_moving_bootstrap import cases, resampling_indices, BootstrapCase
from verify_canonical_bootstrap_intervals import input_payload


def acceleration(pseudovalues: np.ndarray, lengths: np.ndarray) -> float:
    if not np.isfinite(pseudovalues).all():
        return float('nan')
    weights = lengths / lengths.sum()
    delta = weights*(pseudovalues - np.dot(weights,pseudovalues))
    scale = np.max(np.abs(delta))
    if not np.isfinite(scale) or scale <= 64*np.finfo(float).eps:
        return float('nan')
    delta /= scale
    return float(np.sum(delta**3)/(6*np.sum(delta**2)**1.5))


def interval(values: np.ndarray, estimate: float, a: float) -> tuple[np.ndarray,bool]:
    invalid = np.array([np.nan,np.nan])
    if not np.isfinite(a):
        return invalid,False
    tie_tol = 64*np.finfo(float).eps*max(1.,abs(estimate))
    p = (np.count_nonzero(values < estimate-tie_tol) +
         .5*np.count_nonzero(np.abs(values-estimate) <= tie_tol))/len(values)
    if not 0. < p < 1.:
        return invalid,False
    z0 = norm.ppf(p)
    d = z0 + norm.ppf([.025,.975])
    denominators = 1.-a*d
    if np.any(denominators <= 0.):
        return invalid,False
    adjusted = norm.cdf(z0+d/denominators)
    if (not np.isfinite(adjusted).all() or adjusted[0] <= 0 or
            adjusted[1] >= 1 or adjusted[0] >= adjusted[1]):
        return invalid,False
    return np.quantile(values,adjusted,method='linear'),True


def reference(item: BootstrapCase):
    c = item.case
    x = integrated(c)
    nf,nw = x.shape[:2]
    nm = min(len(c.a),len(c.b))
    lengths = np.array([min(item.block_length,nw-k)
                        for k in range(0,nw,item.block_length)])
    indices = resampling_indices(nw,item.block_length,item.n_resamples,
                                 item.seed,item.circular)
    group_result = []
    modes_result = []
    for f in range(nf):
        point_modes,point_group = canonical(x[f].mean(axis=0),c.a,c.b,c.c,
                                           c.ridge,c.shrink)
        bootstrap_modes = []
        bootstrap_groups = []
        for idx in indices:
            roots,coh = canonical(x[f,idx].mean(axis=0),c.a,c.b,c.c,
                                  c.ridge,c.shrink)
            bootstrap_modes.append(roots)
            bootstrap_groups.append(coh)
        bootstrap_modes = np.asarray(bootstrap_modes)
        bootstrap_groups = np.asarray(bootstrap_groups)
        pseudo_groups = np.full(len(lengths),np.nan)
        pseudo_modes = np.full((len(lengths),nm),np.nan)
        if len(lengths) >= 4:
            for k,start in enumerate(range(0,nw,item.block_length)):
                mask = np.ones(nw,dtype=bool)
                mask[start:start+lengths[k]] = False
                leave_modes,leave_group = canonical(x[f,mask].mean(axis=0),
                    c.a,c.b,c.c,c.ridge,c.shrink)
                if (min(point_group,leave_group) <= c.tol or
                    max(point_group,leave_group) >= 1.-c.tol):
                    continue
                remaining = nw-lengths[k]
                pseudo_groups[k] = (nw*point_group-remaining*leave_group)/lengths[k]
                for m in range(nm):
                    if regular(point_modes,m,c.tol) and regular(leave_modes,m,c.tol):
                        pseudo_modes[k,m] = (nw*point_modes[m]-
                            remaining*leave_modes[m])/lengths[k]
        group_a = acceleration(pseudo_groups,lengths)
        group_interval, group_valid = interval(bootstrap_groups,point_group,group_a)
        group_stats = np.array([point_group,bootstrap_groups.std(ddof=1),
            bootstrap_groups.mean()-point_group,*group_interval,group_a])
        group_result.append((group_stats,group_valid))
        row = []
        for m in range(nm):
            mode_a = acceleration(pseudo_modes[:,m],lengths)
            valid_mode = (regular(point_modes,m,c.tol) and all(
                regular(roots,m,c.tol) for roots in bootstrap_modes))
            mode_interval, interval_ok = interval(bootstrap_modes[:,m],
                                                  point_modes[m],mode_a)
            if not valid_mode:
                mode_interval[:] = np.nan
            mode_stats = np.array([point_modes[m],bootstrap_modes[:,m].std(ddof=1),
                bootstrap_modes[:,m].mean()-point_modes[m],*mode_interval,mode_a])
            row.append((mode_stats,valid_mode and interval_ok))
        modes_result.append(row)
    return group_result,modes_result


def main(driver: str, max_cases: int) -> None:
    items = cases()[:max_cases]
    run = subprocess.run([driver],input=input_payload([(item,'bca') for item in items]),
                         capture_output=True,text=True,check=True)
    rows = iter(run.stdout.splitlines())
    max_error = 0.0
    comparisons = 0
    valid_group = 0
    invalid_group = 0
    for j,item in enumerate(items):
        status,nf,nm,nw,bl,nr,nb,wrap = map(int,next(rows).split())
        assert status in (0,2),(j,status)
        assert (nf,nm,nw,bl,nr,nb,wrap) == (
            1 if item.case.bands is not None else len(item.case.frequencies),
            min(len(item.case.a),len(item.case.b)),item.case.samples.shape[1],
            item.block_length,item.n_resamples,
            (item.case.samples.shape[1]-1)//item.block_length+1,int(item.circular))
        expected_group,expected_modes = reference(item)
        for f in range(nf):
            fields = next(rows).split()
            assert int(fields[0]) == f+1
            actual = np.array([float(value) for value in fields[1:-1]])
            expected,valid = expected_group[f]
            assert bool(int(fields[-1])) == valid,(j,f,'group validity',valid)
            assert np.array_equal(np.isnan(actual),np.isnan(expected)), (j,f,'group NaNs',actual,expected)
            if valid:
                valid_group += 1
            else:
                invalid_group += 1
                assert np.isnan(actual[-3:-1]).all(),(j,f,'invalid group bounds')
            diffs = np.abs(actual[np.isfinite(expected)]-expected[np.isfinite(expected)])
            if diffs.size:
                error = float(diffs.max())
                assert error < 4e-11,(j,f,'group',error,actual,expected)
                max_error = max(max_error,error)
                comparisons += len(diffs)
            for m in range(nm):
                fields = next(rows).split()
                assert (int(fields[0]),int(fields[1])) == (f+1,m+1)
                actual = np.array([float(value) for value in fields[2:-1]])
                expected,valid = expected_modes[f][m]
                assert bool(int(fields[-1])) == valid,(j,f,m,'mode validity')
                # When the mode is not isolated, Fortran intentionally skips
                # its moments entirely while the acceleration remains available.
                if not regular(np.array([expected_modes[f][k][0][0]
                    for k in range(nm)]),m,item.case.tol):
                    continue
                assert np.array_equal(np.isnan(actual),np.isnan(expected)),(j,f,m,'mode NaNs')
                diffs = np.abs(actual[np.isfinite(expected)]-expected[np.isfinite(expected)])
                if diffs.size:
                    error = float(diffs.max())
                    assert error < 4e-11,(j,f,m,'mode',error,actual,expected)
                    max_error = max(max_error,error)
                    comparisons += len(diffs)
    print(f'PASS {len(items)} independently regenerated BCa scenarios, '
          f'{comparisons} scalar quantities, {valid_group} valid and '
          f'{invalid_group} invalid group intervals; max error {max_error:.3g}',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    parser.add_argument('--max-cases',type=int,default=30)
    args = parser.parse_args()
    main(args.driver,args.max_cases)
