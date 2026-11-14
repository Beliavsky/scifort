#!/usr/bin/env python3
"""Compile and compare SciFort B-splines against SciPy across randomized cases.

From repository root: python tools/interpolate_reference/compare_bspline_scipy.py
Requires a Fortran 2018 compiler, NumPy, and SciPy; no SciPy runtime dependency.
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.interpolate import BSpline, make_interp_spline

ROOT = Path(__file__).resolve().parents[2]
PEAK_SCALED = 0.0


def cases(seed):
    rng = np.random.default_rng(seed)
    for k in range(0, 8):
        for case in range(70):
            n = int(rng.integers(max(2, k+1), max(2,k+1)+15))
            x = np.r_[rng.uniform(-1, 1), np.cumsum(10.**rng.uniform(-1.,.5,n-1))]
            x[1:] += x[0]
            channels = 1+case % 4
            y = rng.normal(size=(n,channels))
            # Known low-order polynomials, plus varying ordinates.
            if case % 17 == 0:
                y = np.arange(channels)[None,:] + .8 * x[:,None]
            bc = ('not-a-knot','natural','clamped')[case%3] if k == 3 else 'not-a-knot'
            ext = case % 4 != 0
            explicit_knots = case%11 == 0 and k >= 2 and bc == 'not-a-knot'
            t = None
            if explicit_knots:
                t = make_interp_spline(x,y,k=k).t
            target = make_interp_spline(x,y,k=k,bc_type=bc,t=t)
            target.extrapolate = ext
            # two or three repeated knot insertions, avoiding boundary repeats
            knot_in = x[0]+(x[-1]-x[0])*float(rng.uniform(.1,.9))
            insert_m = int(case%3) if k > 0 else 0
            q = np.r_[x[0], x[-1],x[1:-1],np.linspace(x[0]-.3*(x[-1]-x[0]),x[-1]+.3*(x[-1]-x[0]),26)]
            lo, hi = x[0],x[-1]
            ints = np.array([[lo,hi],[lo-.3*(hi-lo),hi+.3*(hi-lo)],
                             [hi,lo],[lo-.5,lo-.5],[lo+.1,lo+.4]])
            if k == 0:
                ints=np.array([[lo+.1*(hi-lo),lo+.8*(hi-lo)],
                               [lo+.3*(hi-lo),lo+.5*(hi-lo)],
                               [lo+.8*(hi-lo),lo+.1*(hi-lo)],
                               [lo+.1*(hi-lo),lo+.1*(hi-lo)]])
            yield dict(mode=1,x=x,y=y,k=k,bc=bc,ext=ext,periodic=False,t=t,
                       q=q,ints=ints,knot_in=knot_in,insert_m=insert_m,target=target)

    for k in range(0, 8):
        for case in range(70):
            nc = int(rng.integers(k+1,k+17))
            lo = float(rng.uniform(-2.,1.)); hi = lo+float(rng.uniform(1.,10.))
            interior = np.sort(rng.uniform(lo+.1,hi-.1,nc-k-1))
            if len(interior) >= 3 and case % 4 == 0:
                interior[1] = interior[0]
            if len(interior) >= 6 and case % 9 == 0:
                interior[2] = interior[1]
            t = np.r_[np.repeat(lo,k+1),interior,np.repeat(hi,k+1)]
            if len(t)!=nc+k+1: raise AssertionError('bad knot length')
            channels=1+case%3
            c = rng.normal(size=(nc,channels))
            ext = case%4 != 0
            periodic = case%9 == 0
            target = BSpline(t,c,k,extrapolate='periodic' if periodic else ext)
            knot_in=float(rng.uniform(lo+.25*(hi-lo),hi-.25*(hi-lo)))
            insert_m=int(case%3) if k>0 else 0
            q = np.r_[lo,hi,interior,np.linspace(lo-.3*(hi-lo),hi+.3*(hi-lo),27)]
            ints=np.array([[lo,hi],[lo-.3*(hi-lo),hi+.3*(hi-lo)],
                           [hi,lo],[lo-.5,lo-.5],[lo+0.1*(hi-lo),lo+.9*(hi-lo)],
                           [lo-2.4*(hi-lo), hi+1.7*(hi-lo)]])
            if k == 0:
                ints=np.array([[lo+.1*(hi-lo),lo+.8*(hi-lo)],
                               [lo+.3*(hi-lo),lo+.5*(hi-lo)],
                               [lo+.8*(hi-lo),lo+.1*(hi-lo)],
                               [lo+.1*(hi-lo),lo+.1*(hi-lo)]])
            yield dict(mode=2,x=np.zeros(nc),y=c,k=k,bc='not-a-knot',
                       ext=ext,periodic=periodic,t=t,q=q,ints=ints,
                       knot_in=knot_in,insert_m=insert_m,target=target)

    # General (not necessarily clamped) knot vectors, with the base interval
    # inside the full knot sequence. These exercise extrapolation polynomials
    # and the span selection at repeated boundary-adjacent knots.
    for k in range(1, 8):
        for case in range(35):
            nc=int(rng.integers(k+1,k+11))
            z=np.cumsum(10.**rng.uniform(-1.,0.5,nc+k))
            t=np.r_[z[0]-.7,z]
            if case%4 == 0 and nc>k+2:
                t[k+2]=t[k+1]  # repeat a knot adjacent to the base interval
            if case%5 == 0 and nc>k+3:
                t[k+3]=t[k+2]
            lo=t[k];hi=t[nc]
            if hi<=lo:continue
            channels=2
            c=rng.normal(size=(nc,channels))
            ext=case%3!=0
            periodic=case%7==0
            target=BSpline(t,c,k,extrapolate='periodic' if periodic else ext)
            knot_in=float(rng.uniform(lo+.2*(hi-lo),hi-.2*(hi-lo)))
            # SciPy/FITPACK rewrites exterior knots when inserting into an
            # unclamped representation. Verify shape invariance separately;
            # no byte-for-byte insertion comparison for these cases.
            insert_m=-1
            q=np.r_[lo,hi,t[k+1:nc],np.linspace(lo-.2*(hi-lo),hi+.2*(hi-lo),21)]
            ints=np.array([[lo+.2*(hi-lo),hi-.1*(hi-lo)],[lo-.2*(hi-lo),hi+.2*(hi-lo)],
                           [hi,lo],[lo-.5,lo-.5]])
            yield dict(mode=2,x=np.zeros(nc),y=c,k=k,bc='not-a-knot',
                       ext=ext,periodic=periodic,t=t,q=q,ints=ints,
                       knot_in=knot_in,insert_m=insert_m,target=target)

    # Extra high-degree checks (degrees 8, 9, 11, 13, 15) on well-scaled
    # coordinates. Avoid extreme extrapolation of ill-conditioned polynomials.
    for k in (8,9,11,13,15):
        for case in range(5):
            n=k+1+case
            x=np.linspace(-1.0,1.0,n)
            y=np.column_stack([np.sin(x*2.0),np.cos(x*1.3)])
            target=make_interp_spline(x,y,k=k)
            q=np.r_[x,np.linspace(-1.,1.,29)]
            ints=np.array([[-.8,.8],[-.5,.3],[.7,-.7],[.2,.2]])
            yield dict(mode=1,x=x,y=y,k=k,bc='not-a-knot',ext=True,
                       periodic=False,t=None,q=q,ints=ints,
                       knot_in=0.15,insert_m=0,target=target)


def arr_string(v):
    return ' '.join(f'{z:.17g}' for z in np.asarray(v).flatten(order='F'))


def input_text(allcases):
    lines = [str(len(allcases))]
    for item in allcases:
        x,y,k,q,ints=item['x'],item['y'],item['k'],item['q'],item['ints']
        t=item['t']
        lines.append(' '.join(map(str,[item['mode'],len(x),y.shape[1],k,
                      ['not-a-knot','natural','clamped'].index(item['bc']),
                      int(item['ext']),int(item['periodic']),0 if t is None else len(t),
                      len(q),len(ints),item['insert_m']])))
        lines.extend([arr_string(x),arr_string(y)])
        if t is not None: lines.append(arr_string(t))
        lines.extend([arr_string(q),arr_string(ints[:,0]),arr_string(ints[:,1]),
                      arr_string([item['knot_in']])])
    return '\n'.join(lines)+'\n'


def one_row(lines, n):
    raw = np.fromstring(next(lines),sep=' ')
    if len(raw) != n:
        raise AssertionError(f'Fortran output length: {len(raw)} expected {n}')
    return raw


def check_array(got,target,key,flags,case_num,atol=4e-9,rtol=2e-9):
    global PEAK_SCALED
    if got.shape != target.shape:
        raise AssertionError(f'{case_num} {key} shape {got.shape}!={target.shape}')
    if not np.allclose(got,target,rtol=rtol,atol=atol,equal_nan=True):
        mask=np.abs(got-target)
        i=np.unravel_index(np.nanargmax(mask),mask.shape)
        raise AssertionError(f'case {case_num} {key} {flags}: {i} got {got[i]} expected {target[i]} '
                             f'error {mask[i]}')
    PEAK_SCALED=max(PEAK_SCALED,float(np.nanmax(np.abs(got-target)/(1.+np.abs(target)),initial=0)))
    return float(np.nanmax(np.abs(got-target),initial=0))


def test_build(flags, seed):
    global PEAK_SCALED
    PEAK_SCALED=0.0
    allcases=list(cases(seed))
    with tempfile.TemporaryDirectory(prefix='scifort_bspline_') as td:
        tmp=Path(td)
        sources=[ROOT/'src/kinds.f90',ROOT/'src/interpolate_bspline.f90',
                 ROOT/'src/interpolate_regular_grid.f90',ROOT/'src/interpolate_piecewise.f90', ROOT/'src/interpolate_smoothing.f90', ROOT/'src/interpolate_barycentric.f90', ROOT/'src/interpolate_global.f90', ROOT/'src/interpolate_kdtree.f90', ROOT/'src/interpolate_rbf.f90', ROOT/'src/interpolate_legacy.f90', ROOT/'src/interpolate_scattered.f90', ROOT/'src/interpolate_clough_tocher.f90', ROOT/'src/interpolate_fitpack_bivariate.f90', ROOT/'src/interpolate_curfit_core.f90', ROOT/'src/interpolate_fitpack_compat.f90', ROOT/'src/interpolate.f90',
                 ROOT/'tools/interpolate_reference/bspline_driver.f90']
        cmd=[os.environ.get('FC','gfortran'),'-std=f2018',*flags,
             '-J'+str(tmp),'-I'+str(tmp),*map(str,sources),'-o',str(tmp/'test_bspline')]
        p=subprocess.run(cmd,cwd=tmp,capture_output=True,text=True)
        if p.returncode:
            raise RuntimeError('Compiler failed:\n'+p.stderr[-5000:])
        p=subprocess.run([str(tmp/'test_bspline')],input=input_text(allcases),
                         capture_output=True,text=True,cwd=tmp)
        if p.returncode:
            raise RuntimeError('Fortran program failed:\n'+p.stderr[-6000:]+'\n'+p.stdout[-2500:])
        lines=iter(p.stdout.splitlines())
        maxerr=0.
        num_values=0
        for case_num,item in enumerate(allcases,1):
            s=item['target'];k=item['k'];q=item['q'];m=item['y'].shape[1]
            label=next(lines).split()
            if label!=['CASE',str(case_num),'0']:
                raise AssertionError(f'case {case_num} {item} status {label}')
            nknots=int(next(lines).split()[1]); t=one_row(lines,nknots)
            if not np.array_equal(t,s.t):
                raise AssertionError(f'{case_num} knot mismatch\n{t}\n{s.t}')
            co_label=next(lines).split();nc=int(co_label[1])
            co=np.array([one_row(lines,nc) for _ in range(m)]).T
            maxerr=max(maxerr,check_array(co,s.c[:nc], 'coefficients',flags,case_num,atol=5e-8,rtol=6e-9))
            assert next(lines).split()==['V',str(len(q)),str(m)]
            for order in [0,1,2]:
                got=np.array([one_row(lines,len(q)) for _ in range(m)]).T
                ref=s(q,nu=order)
                maxerr=max(maxerr,check_array(got,ref,f'values order {order}',flags,case_num,atol=5e-8,rtol=2e-8))
                num_values+=ref.size
            assert next(lines).split()==['I',str(len(item['ints'])),str(m)]
            got=np.array([one_row(lines,len(item['ints'])) for _ in range(m)]).T
            ref=np.array([s.integrate(a,b) for a,b in item['ints']])
            maxerr=max(maxerr,check_array(got,ref,'integrals',flags,case_num,atol=1e-7,rtol=3e-8))
            num_values+=ref.size
            label=next(lines).split()
            if label[0]!='INS': raise AssertionError(label)
            if item['insert_m']>=0:
                if label!=['INS','0']:
                    raise AssertionError(f'case {case_num} insert failure {label}')
                target=(s if item['insert_m']==0 else s.insert_knot(item['knot_in'],item['insert_m']))
                nknots=int(next(lines).split()[1]); t=one_row(lines,nknots)
                maxerr=max(maxerr,check_array(t,target.t,'inserted knot locations',flags,case_num,atol=0,rtol=0))
                nc=int(next(lines).split()[1]); co=np.array([one_row(lines,nc) for _ in range(m)]).T
                maxerr=max(maxerr,check_array(co,target.c[:nc],'inserted coefficients',flags,case_num,atol=5e-8,rtol=1e-8))
                got=np.array([one_row(lines,len(q)) for _ in range(m)]).T
                maxerr=max(maxerr,check_array(got,target(q),'inserted values',flags,case_num,atol=5e-8,rtol=2e-8))
                num_values+=got.size
            label=next(lines).split()
            if label != ['DM','0']: raise AssertionError(f'{case_num} matrix status {label}')
            mat=np.array([one_row(lines,len(s.c)) for _ in range(len(q))])
            ref=BSpline.design_matrix(q,s.t,k,extrapolate=True).toarray()
            maxerr=max(maxerr,check_array(mat,ref,'design matrix',flags,case_num,atol=1e-8,rtol=2e-8))
            num_values+=ref.size
            label=next(lines).split()
            if k>=1:
                try:
                    ref=s.derivative()(q)
                except ValueError:
                    if label!=['DER','1']:
                        raise AssertionError(f'{case_num} repeated-knot derivative status {label}')
                else:
                    if label!=['DER','0']: raise AssertionError(f'{case_num} derivative_spline {label}')
                    got=one_row(lines,len(q))
                    maxerr=max(maxerr,check_array(got,ref[:,0],'derivative object',flags,case_num,atol=5e-8,rtol=2e-8))
            elif label!=['DER','-1']:
                raise AssertionError(label)
            label=next(lines).split()
            if label!=['ANT','0']: raise AssertionError(f'{case_num} antiderivative {label}')
            got=one_row(lines,len(q));ref=s.antiderivative()(q)
            # SciPy 1.17's degree-zero antiderivative has a discontinuity at
            # the final knot; SciFort retains the mathematically continuous integral.
            if k == 0:
                base_hi = s.t[len(s.c)]
                mask = (q >= s.t[0]) & (q < base_hi)
                got, ref = got[mask], ref[mask]
            maxerr=max(maxerr,check_array(got,ref[:,0],'antiderivative object',flags,case_num,atol=5e-8,rtol=2e-8))
            num_values+=len(got)
        trailing=list(lines)
        if trailing:
            raise AssertionError(f'leftover lines {trailing[:4]}')
        print(f'{len(allcases)} SciPy cases, {num_values:,} matrix/value entries, '
              f'max absolute error {maxerr:.3g}, max |error|/(1+|ref|) {PEAK_SCALED:.3g}, '
              f'flags={flags}, seed={seed}',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=20261031)
    args=parser.parse_args()
    for flags in (['-O0','-g','-fcheck=all','-fbacktrace'],['-O3']):
        test_build(flags,args.seed)


if __name__=='__main__':
    main()
