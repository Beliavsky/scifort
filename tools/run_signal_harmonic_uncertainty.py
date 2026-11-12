#!/usr/bin/env python3
"""Check SciFort's harmonic coefficient and curve SEs against NumPy WLS.

Run from any directory: python tools/run_signal_harmonic_uncertainty.py
Requires NumPy and SciPy for the previously established signal reference cases.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile

import numpy as np

from run_signal_harmonic_fit import make_cases, input_text

ROOT = Path(__file__).resolve().parents[1]


def expected_uncertainty(case):
    """Independent full-design WLS using NumPy, not the Fortran QR factors."""
    time, observed, freqs, weights0, new_times, nterms, center, floating, weighted = case
    weights = weights0.copy() if weighted else np.ones_like(time)
    weights /= np.max(weights)
    weights /= np.sum(weights)
    root_weights = np.sqrt(weights)
    shift = float(np.mean(observed)) if center else 0.0
    centered = observed - shift
    p = 2*nterms + int(floating)
    dof = int(np.count_nonzero(weights)) - p
    out = []
    for omega in freqs:
        harmonic = np.column_stack([
            part
            for h in range(1, nterms+1)
            for part in (np.cos(h*omega*time), np.sin(h*omega*time))
        ])
        new_harmonic = np.column_stack([
            part
            for h in range(1, nterms+1)
            for part in (np.cos(h*omega*new_times), np.sin(h*omega*new_times))
        ])
        design = np.column_stack((np.ones(len(time)),harmonic)) if floating else harmonic
        new_design = (np.column_stack((np.ones(len(new_times)),new_harmonic))
                      if floating else new_harmonic)
        weighted_design = root_weights[:,None]*design
        weighted_target = root_weights*centered
        coefficients = np.linalg.lstsq(weighted_design,weighted_target,rcond=None)[0]
        residual = weighted_target-weighted_design@coefficients
        residual_variance = float(residual@residual/dof)
        # The reference forms the full intercept+harmonics Gram matrix.
        # The Fortran calculation instead derives its covariance directly
        # from R and reconstructs intercept/harmonic correlations.
        covariance_fitted = residual_variance*np.linalg.inv(weighted_design.T@weighted_design)
        covariance = np.full((1+2*nterms,1+2*nterms),np.nan)
        if floating:
            covariance[:,:] = covariance_fitted
        else:
            covariance[1:,1:] = covariance_fitted
        std = np.sqrt(np.diag(covariance))
        prediction_var = np.einsum('ij,jk,ik->i',new_design,covariance_fitted,new_design)
        pred_std = np.sqrt(np.maximum(prediction_var,0.0))
        out.append((residual_variance,std,covariance,pred_std))
    return dof,out


def verify(exe: Path, cases, expected):
    completed = subprocess.run([str(exe)],input=input_text(cases),text=True,
                               check=True,capture_output=True)
    stream = iter(completed.stdout.splitlines())
    worst = 0.0
    values = 0
    for j,(case,(dof,results)) in enumerate(zip(cases,expected)):
        _,_,freqs,_,new_times,nterms,_,_,_ = case
        status,reported_dof = map(int,next(stream).split())
        assert (status,reported_dof) == (0,dof), (j,status,reported_dof,dof)
        p = 1+2*nterms
        for k,(variance,se,cov,pred_std) in enumerate(results):
            assert int(next(stream)) == 1, (j,k)
            got_variance = float(next(stream))
            got_se = np.array([float(next(stream)) for _ in range(p)])
            got_cov = np.array([float(next(stream)) for _ in range(p*p)]).reshape((p,p),order='F')
            assert int(next(stream)) == 0, (j,k)
            got_pred_std = np.array([float(next(stream)) for _ in new_times])
            for name,actual,ref in [('noise scale',got_variance,variance),
                                    ('coefficient SE',got_se,se),
                                    ('covariance',got_cov,cov),
                                    ('curve SE',got_pred_std,pred_std)]:
                np.testing.assert_allclose(actual,ref,rtol=5e-9,atol=5e-11,
                                           equal_nan=True,err_msg=f'case {j}, freq {k}, {name}')
                finite = np.isfinite(ref)
                distance = np.abs(np.asarray(actual)-ref)
                worst = max(worst,float(np.max(distance[finite])) if np.any(finite) else 0.0)
                values += int(np.count_nonzero(finite))
    assert next(stream,None) is None
    return values,worst


def build_and_run(directory: Path, optimized: bool, cases, expected):
    directory.mkdir()
    flags = ['-O3','-march=native'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    common = ['gfortran','-std=f2018',*flags,'-Wall','-Wextra',
              '-Werror=implicit-interface','-J',str(directory),'-I',str(directory)]
    sources = [ROOT/'src/kinds.f90',ROOT/'src/signal_irregular.f90']
    for title,tail in [('new native',ROOT/'test/test_signal_harmonic_uncertainty.f90'),
                       ('original native',ROOT/'test/test_signal_harmonic_fit.f90'),
                       ('multiharmonic native',ROOT/'test/test_signal_multiharmonic.f90'),
                       ('NumPy reference',ROOT/'tools/reference_signal_harmonic_uncertainty.f90')]:
        exe = directory/title.replace(' ','_')
        subprocess.run([*common,*(str(f) for f in sources),str(tail),'-o',str(exe)],
                       cwd=directory,check=True,capture_output=True,text=True)
        if title == 'NumPy reference':
            values,worst = verify(exe,cases,expected)
            print(f'  {len(cases)} cases, {values} finite values, max error {worst:.3e}')
        else:
            result = subprocess.run([str(exe)],check=True,capture_output=True,text=True)
            print(' ',result.stdout.strip())


def main():
    cases = make_cases()
    reference = [expected_uncertainty(case) for case in cases]
    with tempfile.TemporaryDirectory(prefix='scifort_harmonic_uncertainty_') as tmp:
        for optimized in (False,True):
            print('optimized -O3:' if optimized else 'checked -fcheck=all:')
            build_and_run(Path(tmp)/('opt' if optimized else 'checked'),optimized,cases,reference)
    print('All harmonic covariance and prediction-SE comparisons passed.')


if __name__ == '__main__':
    main()
