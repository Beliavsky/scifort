#!/usr/bin/env python3
"""Randomized validation of periodic make_interp_spline vs SciPy/NumPy.

Builds SciFort with gfortran, with no FITPACK, BLAS or LAPACK dependency.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import make_interp_spline

ROOT = Path(__file__).resolve().parents[2]


def fmt(a):
    return ' '.join(f'{float(v):.17g}' for v in np.ravel(a, order='F'))


def cases(rng, count):
    for idx in range(count):
        k = 1 + idx % 15
        n = max(k+3, 5) + idx % 9
        spacings = rng.uniform(0.35, 1.0, n-1)
        x = np.r_[0, np.cumsum(spacings)]
        x = -1.5 + 7.0*x/x[-1]
        # Use a smooth periodic signal and optionally small random noise.
        phase = 2*np.pi*(x-x[0])/(x[-1]-x[0])
        channels = 1 + idx % 4
        y = np.empty((n,channels))
        for j in range(channels):
            y[:,j] = (np.sin((1+j%3)*phase+0.4*j) +
                      0.15*np.cos(2*phase-j*0.23) +
                      0.02*rng.normal(size=n))
        y[-1,:] = y[0,:]
        extrap = idx % 5 != 0
        q = np.r_[x, x[0] - 1.5*(x[-1]-x[0]), x[-1]+1.3*(x[-1]-x[0]),
                  rng.uniform(x[0]-12,x[-1]+12,size=20)]
        period = x[-1]-x[0]
        integ = np.array([[x[0],x[-1]], [x[-1],x[0]],
                          [x[0]-1.3*period,x[-1]+2.1*period],
                          [x[1],x[-2]],[x[0],x[0]],
                          [x[-2]+0.4*period,x[-2]+1.8*period]])
        if not extrap:
            integ[-1] = [x[1],x[-2]]
        yield dict(idx=idx,n=n,k=k,x=x,y=y,extrap=extrap,q=q,integ=integ)


def check(got, expected, label, idx, scale_tol):
    if got.shape != expected.shape:
        raise AssertionError(f'{idx}: {label} shapes {got.shape} != {expected.shape}')
    err = np.abs(got-expected)
    scale = np.maximum(1, np.abs(expected))
    normalized = np.where(np.isnan(got) & np.isnan(expected),0,err/scale)
    peak = float(np.nanmax(normalized, initial=0))
    if not np.all(np.isfinite(normalized)) or peak > scale_tol:
        imax = np.unravel_index(int(np.nanargmax(normalized)),normalized.shape)
        raise AssertionError(f'case {idx}: {label} max error {peak:.4g}, location {imax}, '
                             f'Fortran {got[imax]} vs SciPy {expected[imax]}')
    return peak


def validate(flags, testcases):
    with tempfile.TemporaryDirectory() as folder:
        tmp = Path(folder)
        cmd = ['gfortran','-std=f2018','-ffree-line-length-none', *flags,
               '-J'+str(tmp),'-I'+str(tmp),
               str(ROOT/'src/kinds.f90'), str(ROOT/'src/interpolate_bspline.f90'),
               str(ROOT/'tools/interpolate_reference/periodic_bspline_driver.f90'),
               '-o',str(tmp/'compare_periodic')]
        subprocess.run(cmd,check=True,capture_output=True,text=True)
        input_lines = [str(len(testcases))]
        references = []
        for case in testcases:
            n,k,x,y,extrap,q,integ = (case[key] for key in ('n','k','x','y','extrap','q','integ'))
            input_lines.append(f'{n} {y.shape[1]} {k} {len(q)} {len(integ)} {int(extrap)}')
            input_lines.extend((fmt(x), fmt(y), fmt(q), fmt(integ[:,0]), fmt(integ[:,1])))
            spline = make_interp_spline(x,y,k=k,bc_type='periodic')
            spline.extrapolate = 'periodic' if extrap and k >= 2 else extrap
            references.append(spline)
        run = subprocess.run([str(tmp/'compare_periodic')],input='\n'.join(input_lines)+'\n',
                             capture_output=True,text=True)
        if run.returncode:
            raise RuntimeError(f'Fortran exited {run.returncode}:\n{run.stderr[-3000:]}\n{run.stdout[-800:]}')
        lines = iter(run.stdout.splitlines())
        worst = dict(knots=0.,coeff=0.,derivs=0.,integral=0.)
        for caseno, (c,spline) in enumerate(zip(testcases,references),1):
            label = next(lines).split()
            if label != ['CASE', str(caseno),'0']:
                raise AssertionError(f'case {caseno}: status {label}')
            size = int(next(lines).split()[1])
            t = np.fromstring(next(lines),sep=' ')
            assert size == len(t)
            worst['knots'] = max(worst['knots'], check(t,spline.t,'knots',caseno,8e-14))
            clen,channels = map(int,next(lines).split()[1:])
            coeff = np.array([np.fromstring(next(lines),sep=' ') for _ in range(channels)]).T
            assert coeff.shape == (clen,channels)
            worst['coeff'] = max(worst['coeff'],check(coeff,spline.c,'coefficients',caseno,1e-6))
            for order in range(min(3,c['k'])+1):
                mark = next(lines).split()
                assert mark == ['V',str(order),str(len(c['q']))],mark
                v = np.array([np.fromstring(next(lines),sep=' ') for _ in range(channels)]).T
                ref = spline(c['q'],nu=order)
                worst['derivs'] = max(worst['derivs'],check(v,ref,f'deriv={order}',caseno,2e-6))
            mark = next(lines).split()
            assert mark == ['I',str(len(c['integ']))],mark
            got = np.array([np.fromstring(next(lines),sep=' ') for _ in range(channels)]).T
            ref = np.array([spline.integrate(a,b) for a,b in c['integ']])
            worst['integral'] = max(worst['integral'],check(got,ref,'integral',caseno,3e-6))
        print(f"{len(testcases)}/{len(testcases)} periodic spline cases PASS; "
              f"maximum normalized differences {worst}",flush=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cases',type=int,default=90)
    ap.add_argument('--seed',type=int,default=20261101)
    ap.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    args=ap.parse_args()
    rng = np.random.default_rng(args.seed)
    testcases=list(cases(rng,args.cases))
    if args.mode in ('checked','both'):
        print('checked:',end=' ',flush=True)
        validate(['-O0','-g','-fcheck=all','-fbacktrace'],testcases)
    if args.mode in ('optimized','both'):
        print('optimized:',end=' ',flush=True)
        validate(['-O3'],testcases)


if __name__ == '__main__': main()
