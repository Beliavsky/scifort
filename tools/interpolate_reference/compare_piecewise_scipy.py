#!/usr/bin/env python3
"""Independently compare real PPoly/BPoly and conversions against SciPy 1.17.

Examples:
    python tools/interpolate_reference/compare_piecewise_scipy.py --seed 20261031
    python tools/interpolate_reference/compare_piecewise_scipy.py --seed 8675309

Both the checked and optimized GNU Fortran builds are tested by default.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import PPoly, BPoly

ROOT = Path(__file__).resolve().parents[2]


def generate(seed: int, count: int):
    rng = np.random.default_rng(seed)
    result = []
    for i in range(count):
        degree = i % 9
        intervals = 1 + (i % 7)
        channels = 1 + (i % 4)
        spacing = rng.uniform(0.22, 1.8, size=intervals)
        x = np.r_[rng.uniform(-2.5, 2.5), 0.0 + spacing.cumsum()]
        x[1:] += x[0]
        if i % 4 == 1:
            x = x[::-1].copy()
        cp = rng.normal(0, 0.35, size=(degree+1, intervals, channels))
        cb = rng.normal(0, 0.5, size=(degree+1, intervals, channels))
        if i % 9 == 0:
            cp *= 0.0
            cb *= 0.0
        if i % 8 == 0:
            cp[-1] += 2.0
            cb[:] = 0.25
        q = np.r_[x, np.nextafter(x[1], x[-1]), np.nextafter(x[1], x[0]),
                  rng.uniform(min(x)-0.25,max(x)+0.25, size=8),
                  min(x)-0.35, max(x)+0.35]
        # Repeated integration beyond the original range tests periodic cycles
        pairs = [(float(a), float(b)) for a,b in rng.uniform(min(x)-0.35,max(x)+0.35,(5,2))]
        pairs += [(float(x[0]), float(x[-1])),(float(x[-1]),float(x[0])),
                  (float(min(x)-0.35),float(max(x)+0.35))]
        periodic = i % 5 == 0
        extrapolate = i % 7 != 0
        result.append((x,cp,cb,q,pairs,periodic,extrapolate))
    return result


def serialize(cases):
    chunks = [str(len(cases))]
    def fmt(values):
        return ' '.join(f'{float(x):.17g}' for x in np.ravel(values,order='F'))
    for x,cp,cb,q,pairs,periodic,extrapolate in cases:
        k,m,nchan=cp.shape
        chunks.append(f'{k-1} {m} {nchan} {len(q)} {len(pairs)} {int(periodic)} {int(extrapolate)}')
        chunks.append(fmt(x))
        chunks.append(fmt(cp))
        chunks.append(fmt(cb))
        chunks.append(fmt(q))
        for _ in range(4):
            chunks.extend(f'{a:.17g} {b:.17g}' for a,b in pairs)
    return '\n'.join(chunks)+'\n'


def reference_objects(case):
    x,cp,cb,_,_,periodic,extrapolate=case
    flag = 'periodic' if periodic else extrapolate
    p=PPoly(cp,x,extrapolate=flag)
    b=BPoly(cb,x,extrapolate=flag)
    return (p,b,PPoly.from_bernstein_basis(b),BPoly.from_power_basis(p))


def compare(actual,expected,context,rtol=1.5e-8,atol=1.0e-9):
    a=np.asarray(actual,dtype=float)
    e=np.asarray(expected,dtype=float)
    if a.shape != e.shape:
        raise AssertionError(f'{context}: shape {a.shape} != {e.shape}')
    equal=np.isclose(a,e,rtol=rtol,atol=atol,equal_nan=True)
    if not np.all(equal):
        idx=tuple(np.argwhere(~equal)[0])
        raise AssertionError(f'{context} at {idx}: Fortran {a[idx]!r}; SciPy {e[idx]!r}; '
                             f'error={abs(a[idx]-e[idx]):.4g}')
    return int(np.isfinite(e).sum()),float(np.nanmax(np.abs(a-e))) if np.any(np.isfinite(e)) else 0.0


def run(flags, seed, count):
    cases=generate(seed,count)
    with tempfile.TemporaryDirectory(prefix='scifort_piecewise_') as temp:
        wd=Path(temp)
        sources=['src/kinds.f90','src/interpolate_piecewise.f90','src/interpolate_bspline.f90',
                 'src/interpolate_regular_grid.f90','src/interpolate_barycentric.f90', 'src/interpolate_global.f90', 'src/interpolate_kdtree.f90', 'src/interpolate_rbf.f90', 'src/interpolate_legacy.f90', 'src/interpolate_scattered.f90', 'src/interpolate_clough_tocher.f90', 'src/interpolate_fitpack_bivariate.f90', 'src/interpolate_fitpack_compat.f90', 'src/interpolate.f90',
                 'tools/interpolate_reference/piecewise_driver.f90']
        cmd=['gfortran','-std=f2018','-ffree-line-length-none',*flags,'-J'+str(wd),'-I'+str(wd),
             *(str(ROOT/s) for s in sources),'-o',str(wd/'piecewise')]
        cp=subprocess.run(cmd,cwd=wd,capture_output=True,text=True)
        if cp.returncode:
            raise RuntimeError(f'compile failed: {cp.stderr[-3000:]}')
        cp=subprocess.run([str(wd/'piecewise')],input=serialize(cases),capture_output=True,text=True,cwd=wd)
        if cp.returncode:
            raise RuntimeError(f'Fortran driver failed: {cp.stderr[-3000:]}\n{cp.stdout[-1500:]}')
        lines=iter(cp.stdout.splitlines())
        checks=0
        maxerr=0.0
        for caseid,case in enumerate(cases,1):
            x,pc,bc,q,pairs,_,_=case
            header=next(lines).split()
            if header!=['CASE',str(caseid)]:
                raise AssertionError(f'case {caseid}: header: {header}')
            objects=reference_objects(case)
            for j,obj in enumerate(objects):
                for nu in range(5):
                    native = np.array([np.fromstring(next(lines),sep=' ') for _ in q])
                    if isinstance(obj,BPoly) and nu > pc.shape[0]:
                        # SciPy 1.17's BPoly evaluator raises above
                        # degree+1; the zero polynomial is the defined
                        # mathematical derivative, preserving the same
                        # bounds/periodic semantics via SciPy PPoly.
                        expected=PPoly(np.zeros((1,len(x)-1,pc.shape[-1])), x,
                                       extrapolate=obj.extrapolate)(q)
                    else:
                        expected = obj(q,nu=nu)
                    nc,err=compare(native,expected,f'case={caseid} object={j} derivative={nu}',
                                   atol=2e-6 if pc.shape[0]>=7 else 1e-9)
                    checks+=nc
                    maxerr=max(err,maxerr)
                native = np.array([np.fromstring(next(lines),sep=' ') for _ in pairs])
                if x[0] > x[-1] and case[5]:
                    # SciPy 1.17 does not consistently support periodic
                    # descending-axis integration. Use its antiderivative
                    # with an independent periodic-domain decomposition.
                    anti = obj.antiderivative()
                    lo,hi=float(min(x)),float(max(x))
                    period=hi-lo
                    def primitive_at(t):
                        return anti(lo + ((t-lo) % period), extrapolate=True)
                    total=anti(hi,extrapolate=True)-anti(lo,extrapolate=True)
                    expected = np.array([primitive_at(b)-primitive_at(a) +
                                         (np.floor((b-lo)/period)-np.floor((a-lo)/period))*total
                                         for a,b in pairs])
                else:
                    expected = np.array([obj.integrate(a,b) for a,b in pairs])
                nc,err=compare(native,expected,f'case={caseid} object={j} definite integration',
                               rtol=2e-6 if pc.shape[0]>5 else 1.5e-8,
                               atol=2e-6 if pc.shape[0]>5 else 1e-9)
                checks+=nc
                maxerr=max(err,maxerr)
                native_d = np.empty((len(q),pc.shape[-1]))
                native_a = np.empty_like(native_d)
                for v in range(len(q)):
                    native_d[v]=np.fromstring(next(lines),sep=' ')
                    native_a[v]=np.fromstring(next(lines),sep=' ')
                objd=obj.derivative()
                obja=obj.antiderivative()
                nc,err=compare(native_d,objd(q),f'case={caseid} object={j} derivative object',
                               atol=2e-6 if pc.shape[0]>=7 else 1e-9)
                checks+=nc
                maxerr=max(err,maxerr)
                nc,err=compare(native_a,obja(q,extrapolate=True),f'case={caseid} object={j} primitive object',
                               atol=2e-6 if pc.shape[0]>=7 else 1e-9)
                checks+=nc
                maxerr=max(err,maxerr)
        if next(lines,None) is not None:
            raise AssertionError('unexpected trailing output from driver')
        print(f'{flags}: {count} cases, {checks} numerical entries, max abs error {maxerr:.3g}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=20261031)
    parser.add_argument('--cases',type=int,default=180)
    parser.add_argument('--optimized-only',action='store_true')
    args=parser.parse_args()
    builds=[['-O3']] if args.optimized_only else [['-O0','-g','-fcheck=all','-fbacktrace'],['-O3']]
    for flags in builds:
        run(flags,args.seed,args.cases)
