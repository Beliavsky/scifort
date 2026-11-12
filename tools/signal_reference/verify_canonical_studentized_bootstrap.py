#!/usr/bin/env python3
"""Independent NumPy/SciPy reference for block-jackknife studentized bootstrap-t.

Recreates Park-Miller bootstrap starts in Python, calculates conditional
canonical coherences with NumPy/SciPy SVD, and estimates delete-group
pseudovalue variances independently for original and resampled block series.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np

from verify_canonical_jackknife import canonical, integrated, regular
from verify_canonical_moving_bootstrap import BootstrapCase, cases, payload, resampling_indices


def block_standard_errors(windows: np.ndarray, groups: tuple[np.ndarray, np.ndarray, np.ndarray],
                          ridge: float, shrink: float, tol: float, block_length: int):
    """Weighted contiguous delete-block pseudovalues, not Fortran output.

    In particular, the final block can have fewer windows. All delete-group
    estimates are computed afresh using independent NumPy matrix operations.
    """
    a,b,c = groups
    n = len(windows)
    full_modes,full_group = canonical(windows.mean(axis=0),a,b,c,ridge,shrink)
    block_ranges = [np.arange(j, min(j+block_length,n)) for j in range(0,n,block_length)]
    m = len(block_ranges)
    se_modes = np.full(len(full_modes),np.nan)
    se_group = float('nan')
    if m < 4:
        return se_modes,se_group
    weights = np.array([len(rows)/n for rows in block_ranges])
    weights_square = weights**2
    group_pseudo = np.empty(m)
    mode_pseudo = np.full((m,len(full_modes)),np.nan)
    for i,drop in enumerate(block_ranges):
        included = np.ones(n,bool)
        included[drop] = False
        leave_modes,leave_group = canonical(windows[included].mean(axis=0),a,b,c,ridge,shrink)
        if (tol < full_group < 1-tol) and (tol < leave_group < 1-tol):
            group_pseudo[i] = (n*full_group-(n-len(drop))*leave_group)/len(drop)
        else:
            group_pseudo[i] = np.nan
        for k in range(len(full_modes)):
            if regular(full_modes,k,tol) and regular(leave_modes,k,tol):
                mode_pseudo[i,k] = (n*full_modes[k]-(n-len(drop))*leave_modes[k])/len(drop)
    denom = 1.0 - weights_square.sum()
    if np.isfinite(group_pseudo).all():
        mean = weights@group_pseudo
        variance = (weights_square@(group_pseudo-mean)**2)/denom
        if variance > (64*np.finfo(float).eps)**2:
            se_group = float(np.sqrt(variance))
    for k in range(len(full_modes)):
        if np.isfinite(mode_pseudo[:,k]).all():
            mean = weights@mode_pseudo[:,k]
            variance = (weights_square@(mode_pseudo[:,k]-mean)**2)/denom
            if variance > (64*np.finfo(float).eps)**2:
                se_modes[k] = float(np.sqrt(variance))
    return se_modes,se_group


def reference(item: BootstrapCase):
    c = item.case
    x = integrated(c)
    nf,nw = x.shape[:2]
    nm = min(len(c.a),len(c.b))
    draws = resampling_indices(nw,item.block_length,item.n_resamples,item.seed,item.circular)
    group_results = []
    mode_results = []
    for f in range(nf):
        block_args = ((c.a,c.b,c.c),c.ridge,c.shrink,c.tol,item.block_length)
        full_modes,full_group = canonical(x[f].mean(axis=0),c.a,c.b,c.c,c.ridge,c.shrink)
        se_modes0, se_group0 = block_standard_errors(x[f],*block_args)
        groups = np.empty(item.n_resamples)
        modes = np.empty((item.n_resamples,nm))
        group_t = np.full(item.n_resamples,np.nan)
        mode_t = np.full((item.n_resamples,nm),np.nan)
        for r,idx in enumerate(draws):
            sample = x[f,idx]
            modes[r],groups[r] = canonical(sample.mean(axis=0),c.a,c.b,c.c,c.ridge,c.shrink)
            se_mode_star,se_group_star = block_standard_errors(sample,*block_args)
            if np.isfinite(se_group_star) and se_group_star > 0:
                group_t[r] = (groups[r]-full_group)/se_group_star
            for k in range(nm):
                if np.isfinite(se_mode_star[k]) and se_mode_star[k] > 0:
                    mode_t[r,k] = (modes[r,k]-full_modes[k])/se_mode_star[k]
        if np.isfinite(group_t).all() and np.isfinite(se_group0) and se_group0 > 0:
            q = np.quantile(group_t,[.025,.975],method='linear')
            limits = np.clip(full_group - se_group0*q[::-1],0,1)
            valid_group = limits[0] <= limits[1]
        else:
            limits = [np.nan,np.nan]
            valid_group = False
        group_results.append((np.array([full_group,groups.std(ddof=1),
            groups.mean()-full_group,se_group0,*limits]),valid_group))
        current = []
        for k in range(nm):
            isolated = regular(full_modes,k,c.tol) and all(
                regular(modes[r],k,c.tol) for r in range(item.n_resamples))
            if isolated and np.isfinite(mode_t[:,k]).all() and np.isfinite(se_modes0[k]):
                q = np.quantile(mode_t[:,k],[.025,.975],method='linear')
                lims = np.clip(full_modes[k]-se_modes0[k]*q[::-1],0,1)
                valid_mode = lims[0] <= lims[1]
            else:
                lims = [np.nan,np.nan]
                valid_mode = False
            # Existing API withholds all moments for nonisolated modes.
            if isolated:
                values = [full_modes[k],modes[:,k].std(ddof=1),
                          modes[:,k].mean()-full_modes[k],se_modes0[k],*lims]
            else:
                values = [full_modes[k],np.nan,np.nan,se_modes0[k],*lims]
            current.append((np.asarray(values),valid_mode))
        mode_results.append(current)
    return group_results,mode_results


def main(driver: str, max_cases: int) -> None:
    regular_cases = [item for item in cases() if
        (len(item.case.samples[0])-1)//item.block_length+1 >= 4]
    insufficient_cases = [item for item in cases() if
        (len(item.case.samples[0])-1)//item.block_length+1 < 4]
    items = regular_cases[:max_cases] + insufficient_cases[:min(4,max_cases)]
    assert len(items) >= 8
    result = subprocess.run([driver],input=payload(items),text=True,capture_output=True,check=True)
    lines = iter(result.stdout.splitlines())
    maximum_error = 0.0
    checked = 0
    valid = 0
    invalid = 0
    for index,item in enumerate(items):
        status,nf,nm,nw,bl,nr,nb = map(int,next(lines).split())
        assert status in (0,2),(index,status)
        assert (nf,nm,nw,bl,nr,nb) == (
            1 if item.case.bands is not None else len(item.case.frequencies),
            min(len(item.case.a),len(item.case.b)),item.case.samples.shape[1],
            item.block_length,item.n_resamples,(nw-1)//bl+1)
        rg,rm = reference(item)
        for f in range(nf):
            g = next(lines).split()
            assert int(g[0]) == f+1
            actual = np.array([float(v) for v in g[1:-1]])
            expected,gvalid = rg[f]
            assert (int(g[-1]) == 1) == gvalid,(index,f,'group-valid')
            shared = np.isfinite(expected)
            assert np.array_equal(shared,np.isfinite(actual)),(index,f,'group finite')
            if shared.any():
                diff = float(np.max(np.abs(actual[shared]-expected[shared])))
                assert diff < 3e-11,(index,f,'group',diff,actual,expected)
                maximum_error = max(maximum_error,diff)
                checked += int(shared.sum())
            valid += int(gvalid)
            invalid += int(not gvalid)
            for k in range(nm):
                m = next(lines).split()
                assert (int(m[0]),int(m[1])) == (f+1,k+1)
                actual = np.array([float(v) for v in m[2:-1]])
                expected,mvalid = rm[f][k]
                assert (int(m[-1]) == 1) == mvalid,(index,f,k,'mode-valid')
                shared = np.isfinite(expected)
                assert np.array_equal(shared,np.isfinite(actual)),(index,f,k,'mode finite')
                if shared.any():
                    diff = float(np.max(np.abs(actual[shared]-expected[shared])))
                    assert diff < 3e-11,(index,f,k,'mode',diff,actual,expected)
                    maximum_error = max(maximum_error,diff)
                    checked += int(shared.sum())
                valid += int(mvalid)
                invalid += int(not mvalid)
    assert valid > 10 and invalid > 1 and checked > 100,(valid,invalid,checked)
    print(f'PASS {len(items)} independent NumPy bootstrap-t scenarios, {checked} finite '
          f'field comparisons, {valid} valid intervals, {invalid} invalid intervals, '
          f'max error {maximum_error:.3g}',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    parser.add_argument('--max-cases',type=int,default=20)
    args = parser.parse_args()
    main(args.driver,args.max_cases)
