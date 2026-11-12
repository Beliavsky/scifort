#!/usr/bin/env python3
"""Validate canonical moving-bootstrap basic and BC intervals against NumPy/SciPy.

Uses independent Python Park-Miller resampling, independent NumPy canonical
matrix SVD/solves, SciPy standard-normal quantiles and NumPy linear quantiles.
Neither the Fortran interval code nor Fortran-generated draws are imported.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np
from scipy.stats import norm

from verify_canonical_jackknife import canonical, integrated, regular
from verify_canonical_moving_bootstrap import (BootstrapCase, cases,
    payload, resampling_indices)


def input_payload(items: list[tuple[BootstrapCase,str]]) -> str:
    # Reuse the exact previous input format, adding method to each case's header.
    chunks = [str(len(items))]
    for item, method in items:
        lines = payload([item]).splitlines()[1:]
        lines[0] += ' '+method
        chunks.extend(lines)
    return '\n'.join(chunks)+'\n'


def interval(values: np.ndarray, point: float, method: str):
    sorted_values = np.sort(values)
    if method == 'percentile':
        return np.quantile(sorted_values,[0.025,0.975],method='linear'), True
    if method == 'basic':
        extremes = np.quantile(sorted_values,[0.025,0.975],method='linear')
        return np.clip(2*point-extremes[::-1],0,1), True
    tie_tol = 64*np.finfo(float).eps*max(1,abs(point))
    count_less = np.sum(values < point-tie_tol)
    count_equal = np.sum(np.abs(values-point) <= tie_tol)
    p0 = (count_less + 0.5*count_equal)/len(values)
    if not 0 < p0 < 1:
        return np.array([np.nan,np.nan]), False
    z0 = norm.ppf(p0)
    alphas = norm.cdf(2*z0+norm.ppf([.025,.975]))
    return np.quantile(sorted_values,alphas,method='linear'), True


def reference(item: BootstrapCase, method: str):
    case = item.case
    x = integrated(case)
    nf, nw = x.shape[:2]
    nm = min(len(case.a),len(case.b))
    indices = resampling_indices(nw,item.block_length,item.n_resamples,
                                 item.seed,item.circular)
    groups = np.zeros((nf,item.n_resamples))
    modes = np.zeros((nf,item.n_resamples,nm))
    gresults = []
    mresults = []
    for f in range(nf):
        root0,coh0 = canonical(x[f].mean(axis=0),case.a,case.b,case.c,
                              case.ridge,case.shrink)
        for r,idx in enumerate(indices):
            modes[f,r],groups[f,r] = canonical(x[f,idx].mean(axis=0),
                case.a,case.b,case.c,case.ridge,case.shrink)
        ci,ok = interval(groups[f],coh0,method)
        gresults.append((np.array([coh0,groups[f].std(ddof=1),
                                  groups[f].mean()-coh0,*ci]),ok))
        this_modes = []
        for k in range(nm):
            valid = (regular(root0,k,case.tol) and
                all(regular(modes[f,r],k,case.tol) for r in range(item.n_resamples)))
            ci,ok = interval(modes[f,:,k],root0[k],method)
            this_modes.append((np.array([root0[k],modes[f,:,k].std(ddof=1),
                modes[f,:,k].mean()-root0[k],*ci]),valid and ok))
        mresults.append(this_modes)
    return gresults,mresults


def main(driver: str, max_cases: int = 40) -> None:
    # Shuffle the method across cases while including all methods and both block types.
    samples = cases()[:max_cases]
    items = [(case,method) for case in samples
             for method in ('percentile','basic','bc')]
    output = subprocess.run([driver],input=input_payload(items),text=True,
                            capture_output=True,check=True)
    lines = iter(output.stdout.splitlines())
    max_error = 0.0
    compared = 0
    invalid_bc = 0
    for index,(item,method) in enumerate(items):
        header = list(map(int,next(lines).split()))
        status,nf,nm,nw,bl,nr,nb,wrap = header
        assert status in (0,2),(index,method,'status',status)
        assert (nf,nm,nw,bl,nr,nb,wrap) == (
            1 if item.case.bands is not None else len(item.case.frequencies),
            min(len(item.case.a),len(item.case.b)),item.case.samples.shape[1],
            item.block_length,item.n_resamples,
            (item.case.samples.shape[1]-1)//item.block_length+1,int(item.circular))
        rg,rm = reference(item,method)
        for f in range(nf):
            fields = next(lines).split()
            assert int(fields[0]) == f+1
            actual_valid = bool(int(fields[-1]))
            expected,valid = rg[f]
            assert actual_valid == valid,(index,f,method,'group validity')
            values = np.array([float(v) for v in fields[1:-1]])
            if not valid:
                assert np.isnan(values[-2:]).all(),(index,f,method,'invalid bounds')
                invalid_bc += method == 'bc'
            else:
                diff = np.max(np.abs(values-expected))
                assert diff < 4e-12,(index,f,method,diff,'group')
                max_error = max(max_error,float(diff))
                compared += 1
            for k in range(nm):
                fields = next(lines).split()
                assert (int(fields[0]),int(fields[1])) == (f+1,k+1)
                actual_valid = bool(int(fields[-1]))
                expected,valid = rm[f][k]
                assert actual_valid == valid,(index,f,k,method,'root validity')
                values = np.array([float(v) for v in fields[2:-1]])
                if not valid:
                    assert np.isnan(values[-2:]).all(),(index,f,k,method,'mode bounds')
                else:
                    diff = np.max(np.abs(values-expected))
                    assert diff < 4e-12,(index,f,k,method,diff,'root')
                    max_error = max(max_error,float(diff))
                    compared += 1
    assert compared > 2*max_cases
    print(f'PASS {len(items)} independent NumPy/SciPy interval cases, '
          f'{compared} valid scalar comparisons, {invalid_bc} unidentifiable '
          f'BC aggregate intervals, max error {max_error:.3g}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('driver')
    parser.add_argument('--max-cases',type=int,default=40)
    args = parser.parse_args()
    main(args.driver,args.max_cases)
