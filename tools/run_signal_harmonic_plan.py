#!/usr/bin/env python3
"""Validate prepared multiharmonic fitting against independent NumPy QR.

Usage: python tools/run_signal_harmonic_plan.py
Requires GNU Fortran, NumPy, and SciPy. Builds in temporary directories,
with both checked and optimized Fortran builds; no generated project files.
"""
from __future__ import annotations

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.signal import lombscargle

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261029)


def make_cases():
    cases = []
    for i in range(128):
        nh = 1 + i % 4
        n = int(RNG.integers(40, 92))
        m = int(RNG.integers(2, 8))
        nch = int(RNG.integers(2, 6))
        x = np.sort(RNG.uniform(-5.0, 20.0, size=n))
        f = RNG.uniform(0.32, 3.8, size=m)
        if i % 9 == 0:
            f[0] = 0.0  # Design rank deficiency confined to one bin.
        w = RNG.uniform(0.1, 2.8, size=n)
        weighted = (i % 5 != 0)
        if i % 7 == 0:
            w[:n // 8] = 0.0
        if i % 13 == 0:
            w *= 1e5
        floating = (i % 2 == 0)
        center = (i % 3 == 0)
        covariance = (i % 4 != 0)
        y = RNG.normal(0.0, 0.18, (n, nch))
        y += RNG.normal(0, 1.5, nch)
        for h in range(1, nh+1):
            y += np.cos(h*f[-1]*x[:,None])*RNG.normal(0, 0.8, nch)
            y += np.sin(h*f[-1]*x[:,None])*RNG.normal(0, 0.8, nch)
        cases.append((x,f,w,y,nh,weighted,floating,center,covariance))
    return cases


def input_text(cases):
    lines = [str(len(cases))]
    for x,f,w,y,nh,weighted,floating,center,covariance in cases:
        lines.append(f'{len(x)} {len(f)} {y.shape[1]} {nh} '
                     f'{int(weighted)} {int(floating)} {int(center)} {int(covariance)}')
        for a in [x,f,w,*[y[:,j] for j in range(y.shape[1])]]:
            lines.append(' '.join(f'{v:.17e}' for v in a))
    return '\n'.join(lines)+'\n'


def expected_fit(x, omega, w0, y, nh, weighted, floating, center, cov_requested):
    w = w0.copy() if weighted else np.ones(len(x))
    w /= np.sum(w)
    dof = int(np.count_nonzero(w>0))-2*nh-int(floating)
    shift = float(y.mean()) if center else 0.0
    baseline = float(w @ (y-shift)) if floating else 0.0
    target = y-shift-baseline
    yy = float(w @ (target*target))
    basis = np.column_stack([v for h in range(1,nh+1) for v in
                             (np.cos(h*omega*x),np.sin(h*omega*x))])
    means = w @ basis if floating else np.zeros(2*nh)
    design = (basis-means)*np.sqrt(w)[:,None]
    if np.linalg.matrix_rank(design) < 2*nh:
        return dof,None
    q,r = np.linalg.qr(design,mode='reduced')
    coeff = np.linalg.solve(r,q.T@(np.sqrt(w)*target))
    intercept = shift+baseline-float(means @ coeff)
    model = intercept+basis@coeff
    rss = float(w @ ((y-model)**2))
    improvement = float(np.linalg.norm(q.T@(np.sqrt(w)*target))**2)
    improvement = min(yy,max(0,improvement))
    values = [intercept,*coeff,0.5*len(x)*improvement,improvement/yy,np.sqrt(rss)]
    if nh == 1 and omega != 0:
        kw = dict(precenter=center, floating_mean=floating)
        if weighted:
            kw['weights'] = w0
        sp = float(np.asarray(lombscargle(x,y,np.array([omega]),**kw)))
        np.testing.assert_allclose(values[1+2*nh],sp,rtol=2e-8,atol=2e-8)
    if cov_requested:
        variance = rss/dof
        invr = np.linalg.solve(r,np.eye(2*nh))
        harmonic_cov = variance*invr@invr.T
        covariance = np.full((1+2*nh,1+2*nh),np.nan)
        covariance[1:,1:] = harmonic_cov
        if floating:
            covariance[0,0] = variance+means@harmonic_cov@means
            covariance[0,1:] = -(harmonic_cov@means)
            covariance[1:,0] = covariance[0,1:]
        standard_error = np.sqrt(np.maximum(0,np.diag(covariance)))
        values.extend([variance,*standard_error,*covariance.ravel(order='F')])
    return dof, values


def verify(exe, cases):
    raw = subprocess.run([str(exe)], input=input_text(cases),
                         text=True, capture_output=True, check=True).stdout
    stream = iter(raw.splitlines())
    checked = 0
    worst = 0.0
    badbins = 0
    for icase, case in enumerate(cases):
        x,f,w,y,nh,weighted,floating,center,covariance = case
        expect_plan_status = 2 if np.any(f == 0) else 0
        status,channels = map(int,next(stream).split())
        assert (status,channels) == (expect_plan_status,y.shape[1]), (icase,status)
        for j in range(channels):
            channel_status,reported_dof = map(int,next(stream).split())
            assert channel_status == expect_plan_status, (icase,j,channel_status)
            for k,omega in enumerate(f):
                dof,values = expected_fit(x,omega,w,y[:,j],nh,weighted,
                                          floating,center,covariance)
                assert dof == reported_dof, (icase,j,k,dof,reported_dof)
                good = int(next(stream))
                if values is None:
                    assert good == 0, (icase,j,k)
                    badbins += 1
                    continue
                assert good == 1, (icase,j,k)
                got = np.array([float(next(stream)) for _ in values])
                ref = np.array(values)
                np.testing.assert_allclose(got,ref,rtol=8e-9,atol=6e-10,
                                           equal_nan=True,
                                           err_msg=f'case {icase} channel {j} freq {k}')
                finite = np.isfinite(ref)
                if np.any(finite):
                    worst = max(worst,float(np.max(np.abs(got[finite]-ref[finite]))))
                    checked += int(np.count_nonzero(finite))
    assert next(stream,None) is None, 'extra output'
    return checked,worst,badbins


def run_build(directory, optimized, cases):
    directory.mkdir()
    flags = ['-O3','-march=native'] if optimized else ['-O0','-g','-fcheck=all','-fbacktrace']
    common = ['gfortran','-std=f2018',*flags,'-Werror=implicit-interface',
              '-J',str(directory),'-I',str(directory)]
    src = [ROOT/'src/kinds.f90', ROOT/'src/signal_irregular.f90']
    objects = []
    for source in src:
        obj = directory/(source.stem+'.o')
        subprocess.run([*common,'-c',str(source),'-o',str(obj)],
                       cwd=directory,check=True,capture_output=True,text=True)
        objects.append(obj)
    for title,main in [('prepared native',ROOT/'test/test_signal_harmonic_plan.f90'),
                       ('earlier fit native',ROOT/'test/test_signal_harmonic_fit.f90'),
                       ('earlier uncertainty native',ROOT/'test/test_signal_harmonic_uncertainty.f90'),
                       ('NumPy reference',ROOT/'tools/reference_signal_harmonic_plan.f90')]:
        exe = directory/title.replace(' ','_')
        subprocess.run([*common,str(main),*(str(obj) for obj in objects),'-o',str(exe)],
                       cwd=directory,check=True,capture_output=True,text=True)
        if title == 'NumPy reference':
            values,worst,singular = verify(exe,cases)
            print(f'  {len(cases)} cases, {values} finite values, '
                  f'{singular} singular bins, max abs difference {worst:.3e}')
        else:
            result = subprocess.run([str(exe)],capture_output=True,text=True,check=True)
            print(' ',result.stdout.strip())


def main():
    cases = make_cases()
    with tempfile.TemporaryDirectory(prefix='scifort_harmonic_plan_') as temp:
        for mode in [False,True]:
            print('optimized -O3:' if mode else 'checked -fcheck=all:')
            run_build(Path(temp)/('optimized' if mode else 'checked'),mode,cases)
    print('All prepared harmonic fit comparisons passed.')


if __name__ == '__main__':
    main()
