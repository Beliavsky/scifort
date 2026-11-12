#!/usr/bin/env python3
"""Compare SciFort's three new spatial distance metrics with SciPy 1.17+.

Builds the Fortran driver twice (bounds checked and -O3); no BLAS or FPM needed.
    python tools/spatial_reference/compare_covariance_distance_scipy.py
    python tools/spatial_reference/compare_covariance_distance_scipy.py --seed 84633

For Mahalanobis/seuclidean, default pdist and cdist statistics differ: pdist
uses XA and cdist uses the union of XA and XB. Explicit VI/V are compared as
well. SciPy's compiled jensenshannon cdist/pdist do not accept `base`; thus
base=2 is validated using scipy.spatial.distance.jensenshannon pointwise.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.spatial.distance import cdist, jensenshannon, pdist

ROOT = Path(__file__).resolve().parents[2]


def cases(seed: int) -> list[tuple]:
    rng = np.random.default_rng(seed)
    output = []
    for metric in ('seuclidean', 'mahalanobis', 'jensenshannon'):
        for ndim in (1, 2, 3, 5, 8):
            for sample in range(10):
                n = int(rng.integers(ndim+3, ndim+10))
                m = int(rng.integers(2, 8))
                if metric == 'jensenshannon':
                    xa = rng.gamma(1.1, 1.0, (n,ndim))
                    xb = rng.gamma(1.1, 1.0, (m,ndim))
                    if sample%3 == 0:
                        xa *= rng.uniform(0.01, 100, (n,1))
                        xb *= rng.uniform(0.01, 100, (m,1))
                    if sample%3 == 1:
                        xa = np.where(rng.random((n,ndim)) < 0.25, 0, xa)
                        xb = np.where(rng.random((m,ndim)) < 0.25, 0, xb)
                        xa[:,0] += 1e-8
                        xb[:,0] += 1e-8
                    if sample%5 == 0:
                        xb[0] = xa[0]
                    bases = (0.0,2.0)
                    options = (0,)
                else:
                    a = rng.normal(size=(ndim,ndim))
                    if ndim>1:
                        a += ndim*np.eye(ndim)
                    xa = rng.normal(size=(n,ndim)) @ a
                    xb = rng.normal(size=(m,ndim)) @ a
                    if sample%3 == 0:
                        xa *= 0.01
                        xb *= 0.01
                    elif sample%3 == 1:
                        xa *= 100
                        xb *= 100
                    if sample%5 == 0:
                        xb[0] = xa[0]
                    options = (0,1)
                    bases = (0.0,)
                for option in options:
                    variance, inv_cov = None, None
                    if option == 1 and metric == 'seuclidean':
                        variance = rng.uniform(0.3, 5.0, ndim)
                    elif option == 1 and metric == 'mahalanobis':
                        a = rng.normal(size=(ndim,ndim))
                        inv_cov = a.T@a + np.eye(ndim)*0.4
                    for base in bases:
                        output.append((metric,xa,xb,option,variance,inv_cov,base))
    return output


def run_fortran(executable: Path, inputs: list[tuple]) -> list[str]:
    lines = [str(len(inputs))]
    for metric,xa,xb,option,variance,inv_cov,base in inputs:
        lines.append(f'{len(xa)} {len(xb)} {xa.shape[1]} {metric} {option} {base:.17g}')
        if variance is not None:
            lines.append(' '.join(f'{v:.17e}' for v in variance))
        if inv_cov is not None:
            for row in inv_cov:
                lines.append(' '.join(f'{v:.17e}' for v in row))
        for row in np.vstack((xa,xb)):
            lines.append(' '.join(f'{v:.17e}' for v in row))
    proc = subprocess.run([str(executable)], input='\n'.join(lines)+'\n', text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return proc.stdout.splitlines()


def verify(lines: list[str], inputs: list[tuple]) -> tuple[int,float]:
    cursor = 0
    count = 0
    max_err = 0.0
    for case,(metric,xa,xb,option,variance,inv_cov,base) in enumerate(inputs):
        status = lines[cursor].split()
        cursor += 1
        if status != ['0','0']:
            raise AssertionError(f'{case} ({metric},{option}) invalid Fortran status {status}')
        n,m = len(xa),len(xb)
        nc, npairs = n*m, n*(n-1)//2
        got_c = np.array([float(s) for s in lines[cursor:cursor+nc]]).reshape(n,m)
        cursor += nc
        got_p = np.array([float(s) for s in lines[cursor:cursor+npairs]])
        cursor += npairs
        kwargs = {}
        if variance is not None: kwargs['V'] = variance
        if inv_cov is not None: kwargs['VI'] = inv_cov
        if metric != 'jensenshannon':
            exp_c = cdist(xa, xb,metric=metric,**kwargs)
            exp_p = pdist(xa,metric=metric,**kwargs)
        else:
            logbase = None if base == 0.0 else base
            exp_c = np.array([[jensenshannon(row_a,row_b,base=logbase)
                               for row_b in xb] for row_a in xa])
            exp_p = np.array([jensenshannon(xa[i],xa[j],base=logbase)
                               for i in range(n) for j in range(i+1,n)])
        # SciPy's C kernel can return spurious NaNs for JS with a single
        # feature (the mathematical distance is identically zero).
        if metric=='jensenshannon' and xa.shape[1]==1:
            exp_c = np.zeros_like(exp_c)
            exp_p = np.zeros_like(exp_p)
        for what,actual,expected in (('cdist',got_c,exp_c),('pdist',got_p,exp_p)):
            if not np.all(np.isfinite(expected)):
                raise AssertionError(f'nonfinite SciPy reference for case {case}')
            if not np.all(np.isfinite(actual)):
                raise AssertionError(f'nonfinite Fortran output for case {case}')
            if not np.allclose(actual,expected,rtol=1e-10,atol=1e-11):
                error = np.abs(actual-expected)
                at = np.unravel_index(np.argmax(error),error.shape)
                raise AssertionError(f'{case} ({metric},{option},{base}) {what} mismatch at {at}: '
                                     f'Fortran {actual[at]}, SciPy {expected[at]}, '
                                     f'abs error {error[at]}')
            err = np.abs(actual-expected)/np.maximum(1.0,np.abs(expected))
            max_err = max(max_err, float(err.max(initial=0.0)))
            count += actual.size
    if cursor != len(lines):
        raise AssertionError(f'Extra Fortran output lines: {len(lines)-cursor}')
    return count,max_err


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=20261102)
    parser.add_argument('--compiler',default=os.environ.get('FC','gfortran'))
    args = parser.parse_args()
    scenarios = cases(args.seed)
    for label,flags in (('checked',['-O0','-g','-fcheck=all','-fbacktrace']),('optimized',['-O3'])):
        with tempfile.TemporaryDirectory(prefix='scifort-spatial-covariance-') as tmp:
            path = Path(tmp)
            exe = path/'covariance_distance_driver'
            sources = [ROOT/'src'/f'{name}.f90' for name in
                       ('kinds','interpolate_kdtree','spatial_distance','spatial_hausdorff','spatial_slerp','spatial_procrustes','spatial_voronoi','spatial_convex_hull','spatial_halfspace','spatial')]
            sources.append(ROOT/'tools/spatial_reference/covariance_distance_driver.f90')
            subprocess.run([args.compiler,'-std=f2018','-ffree-line-length-none',*flags,
                            '-J',str(path),'-o',str(exe),*map(str,sources)], cwd=path, check=True)
            counts,err = verify(run_fortran(exe,scenarios),scenarios)
            print(f'{label}: {len(scenarios)} randomized cases, {counts} compared values, '
                  f'max scale-relative error {err:.3g}: PASS',flush=True)


if __name__ == '__main__':
    main()
