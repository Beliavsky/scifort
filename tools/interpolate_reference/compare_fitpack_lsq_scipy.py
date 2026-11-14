#!/usr/bin/env python3
"""SciPy LSQUnivariateSpline / InterpolatedUnivariateSpline numerical oracle.

Builds real GNU Fortran binaries (checked and -O3); checks coefficients, knots,
weighted residuals indirectly through coefficients, values, derivatives, integrals.
Includes multichannel and nonuniform/weighted sample sets. No plots.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import LSQUnivariateSpline, InterpolatedUnivariateSpline

ROOT=Path(__file__).resolve().parents[2]
SOURCES=('kinds','interpolate_piecewise','interpolate_bspline','interpolate_fitpack_bivariate','interpolate_curfit_core', 'interpolate_fitpack_compat','interpolate_smoothing',
         'interpolate_regular_grid','interpolate_barycentric','interpolate_global',
         'interpolate_kdtree','interpolate_rbf','interpolate_legacy',
         'interpolate_scattered','interpolate_clough_tocher','interpolate')
DRIVER=ROOT/'tools/interpolate_reference/fitpack_lsq_driver.f90'

def build(work:Path, flags:list[str])->Path:
    common=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    for name in SOURCES:
        cmd=common+['-I'+str(work),'-J'+str(work),'-c',str(ROOT/'src'/f'{name}.f90'),
                    '-o',str(work/f'{name}.o')]
        run=subprocess.run(cmd,cwd=work,text=True,capture_output=True)
        if run.returncode:raise RuntimeError(' '.join(cmd)+'\n'+run.stderr[-4000:])
    exe=work/'fitpack_driver'
    run=subprocess.run(common+['-I'+str(work),str(DRIVER),
                       *(str(work/f'{name}.o') for name in SOURCES),'-o',str(exe)],
                       cwd=work,text=True,capture_output=True)
    if run.returncode:raise RuntimeError(run.stderr[-4000:])
    return exe

def case(exe:Path,rng:np.random.Generator,i:int):
    k=1+i%5
    interpolating=(i%6==0)
    n=k+8+(i*7)%55
    x=np.cumsum(rng.uniform(.12,1.5,size=n))
    x=(x-x[0])/max(1.,(x[-1]-x[0]))*7.0-2.0
    nc=1+(i//3)%4
    w=np.exp(rng.uniform(-1.,1.,n))
    y=np.column_stack([
        np.sin((j+1)*0.7*x)+0.15*np.cos(2.5*x)+0.04*x**2+
        0.12*rng.normal(size=n) for j in range(nc)])
    if i%7==0:y[:,0]=1.2+0.7*x
    nt=0 if interpolating else (i//5)%min(9,max(1,(n-k-1)//4))+1
    if interpolating:
        interior=np.array([],dtype=float)
    else:
        # interior knots are interior observation sites, separated so every
        # basis column has enough measurements for full-rank LSQ fitting.
        ids=np.linspace(k+1,n-k-2,nt,dtype=int)
        interior=x[ids]
        if len(np.unique(interior))!=nt:
            raise AssertionError('bad knot construction')
    bbox=np.array([x[0],x[-1]])
    if i%7==2 and not interpolating:
        bbox[0]-=.2
        bbox[1]+=.3
    q=np.r_[np.linspace(x[0]-.25,x[-1]+.25,15), x,
            (x[1:]+x[:-1])*.5]
    lines=[f'{n} {nc} {nt} {k} {len(q)} {int(interpolating)}',
           ' '.join(f'{z:.17g}' for z in bbox),
           ' '.join(f'{z:.17g}' for z in x)]
    if nt: lines.append(' '.join(f'{z:.17g}' for z in interior))
    lines.append(' '.join(f'{z:.17g}' for z in w))
    lines.extend(' '.join(f'{z:.17g}' for z in row) for row in y)
    lines.append(' '.join(f'{z:.17g}' for z in q))
    run=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',
                       text=True,capture_output=True)
    if run.returncode: raise RuntimeError(run.stderr[-1000:])
    result=run.stdout.splitlines()
    status=int(result[0])
    if status:raise ValueError(f'case={i} got status {status}')
    nkn,ncoeff,ncols=map(int,result[1].split())
    assert ncols==nc
    knots=np.fromstring(result[2],sep=' ')
    assert len(knots)==nkn
    coeff=np.vstack([np.fromstring(line,sep=' ') for line in result[3:3+ncoeff]])
    base=3+ncoeff
    vals=np.vstack([np.fromstring(line,sep=' ') for line in result[base:base+len(q)*nc]])
    vals=vals.reshape(len(q),nc,4)
    integ=np.fromstring(result[base+len(q)*nc],sep=' ')
    expvals=np.empty_like(vals)
    expcoeff=np.empty_like(coeff)
    expintegral=np.empty_like(integ)
    for j in range(nc):
        if interpolating:
            ref=InterpolatedUnivariateSpline(x,y[:,j],k=k)
        else:
            ref=LSQUnivariateSpline(x,y[:,j],interior,w=w,k=k,bbox=tuple(bbox))
        # FITPACK pads the b-spline coefficient array beyond ncoeff; ignore.
        ref_knots=ref._eval_args[0]
        if len(ref_knots)!=len(knots):raise AssertionError('knot count mismatch')
        np.testing.assert_allclose(knots,ref_knots,rtol=1e-14,atol=1e-14)
        expcoeff[:,j]=ref.get_coeffs()
        for der in range(4):
            expvals[:,j,der] = ref(q,nu=der) if der<=k else 0.0
        expintegral[j]=ref.integral(x[0],x[-1])
    normalized=lambda a,b:np.max(np.abs(a-b)/np.maximum(1.,np.abs(b)))
    errors=np.array([normalized(coeff,expcoeff),normalized(vals,expvals),
                     normalized(integ,expintegral)])
    if np.any(errors>3e-8):
        raise AssertionError(f'case {i}: k={k}, interp={interpolating}, nt={nt}, '
                             f'maxscaled coef/eval/int = {errors}')
    return errors,coeff.size+vals.size+integ.size

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--cases',type=int,default=140)
    parser.add_argument('--seed',type=int,default=20261101)
    parser.add_argument('--mode',choices=['both','checked','optimized'],default='both')
    args=parser.parse_args()
    modes=[]
    if args.mode!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.mode!='checked':modes.append(('optimized',['-O3']))
    for mode,flags in modes:
        with tempfile.TemporaryDirectory() as d:
            exe=build(Path(d),flags)
            rng=np.random.default_rng(args.seed)
            err=np.zeros(3)
            total=0
            for i in range(args.cases):
                try:e,count=case(exe,rng,i)
                except Exception as ex:raise RuntimeError(f'{mode} seed={args.seed} case={i}') from ex
                err=np.maximum(err,e)
                total+=count
            print(f'{mode}: {args.cases} cases, {total} results compared; '
                  f'maxscaled(coef/eval/integral)={err}; PASS',flush=True)

if __name__=='__main__':main()
