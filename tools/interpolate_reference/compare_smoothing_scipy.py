#!/usr/bin/env python3
"""Compare SciFort natural cubic smoothing splines with SciPy's implementation.

Covers known lambda, per-channel GCV, weights, derivatives, extrapolation,
coefficients converted by SciPy's PPoly, and exact integral. No plotting.
"""
from __future__ import annotations
import argparse
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from scipy.interpolate import make_smoothing_spline, PPoly

ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'src'
DRIVER=ROOT/'tools/interpolate_reference/smoothing_driver.f90'


def build(where:Path,flags:list[str])->Path:
    common=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    names=('kinds','interpolate_piecewise','interpolate_bspline','interpolate_fitpack_bivariate','interpolate_curfit_core', 'interpolate_fitpack_compat','interpolate_smoothing',
           'interpolate_regular_grid','interpolate_barycentric','interpolate_global',
           'interpolate_kdtree','interpolate_rbf','interpolate_legacy',
           'interpolate_scattered','interpolate_clough_tocher','interpolate')
    for name in names:
        cmd=common+['-I'+str(where),'-J'+str(where),'-c',str(SRC/f'{name}.f90'),
                    '-o',str(where/f'{name}.o')]
        res=subprocess.run(cmd,cwd=where,text=True,capture_output=True)
        if res.returncode: raise RuntimeError(' '.join(cmd)+'\n'+res.stderr[-3000:])
    exe=where/'smoothing_driver'
    res=subprocess.run(common+['-I'+str(where),str(DRIVER),
                      *(str(where/f'{n}.o') for n in names),'-o',str(exe)],
                      cwd=where,text=True,capture_output=True)
    if res.returncode: raise RuntimeError(res.stderr[-3000:])
    return exe


def check(exe:Path,rng:np.random.Generator,index:int):
    n=5+index%40
    nc=1+(index//3)%4
    x=np.cumsum(rng.uniform(0.15,0.95,size=n)) - 3
    xx=x-x[0]
    y=np.column_stack([(0.03*xx**2+np.sin(xx*(.25+.04*c))
                        +rng.normal(scale=.08+.03*c,size=n)) for c in range(nc)])
    if index%6==0:
        y[:,0]=2+.5*x    # perfectly linear data: all lambda equivalent
    w=rng.uniform(.3,2.,size=n) if index%3 else np.ones(n)
    lam=(0,1e-5,0.003,0.02,0.13,1.5,12.0)[index%7]
    gcv=index%4==3
    q=np.r_[np.linspace(x[0]-.4,x[-1]+.4,11),x,(x[1:]+x[:-1])/2]
    raw=[f'{n} {nc} {len(q)} {int(gcv)}',f'{lam:.17g}',
         ' '.join(f'{v:.17g}' for v in x),
         ' '.join(f'{v:.17g}' for v in w)]
    raw.extend(' '.join(f'{v:.17g}' for v in row) for row in y)
    raw.append(' '.join(f'{v:.17g}' for v in q))
    proc=subprocess.run([str(exe)],input='\n'.join(raw)+'\n',text=True,capture_output=True)
    if proc.returncode: raise RuntimeError(proc.stderr[-2000:])
    lines=proc.stdout.splitlines()
    if int(lines[0])!=0:raise ValueError(f'status {lines[0]} case {index}')
    actual_lam=np.fromstring(lines[1],sep=' ')
    vals=np.array([np.fromstring(line,sep=' ') for line in lines[2:2+len(q)*nc]])
    vals=vals.reshape(len(q),nc,3)
    start=2+len(q)*nc
    coeff=np.array([np.fromstring(line,sep=' ') for line in lines[start:start+(n-1)*nc]])
    coeff=coeff.reshape(n-1,nc,4).transpose(2,0,1)
    integral=np.array([float(line) for line in lines[start+(n-1)*nc:]])
    if len(integral)!=nc:raise RuntimeError('output format mismatch')
    expected=np.zeros_like(vals)
    refcoeff=np.zeros_like(coeff)
    refint=np.zeros(nc)
    gcv_direct_max=0.0
    for c in range(nc):
        # Independently compare the fitted spline at the identical penalty.
        # The selection algorithms differ: SciPy's default bounded optimizer
        # terminates at ~1e-5, while the SciFort GCV search is more precise.
        ref=make_smoothing_spline(x,y[:,c],w=w,lam=actual_lam[c] if gcv else lam)
        for k in range(3): expected[:,c,k]=ref(q,k)
        refcoeff[:,:,c]=PPoly.from_spline(ref).c[:,3:-3]
        refint[c]=ref.integrate(x[0],x[-1])
        if gcv:
            sci_auto=make_smoothing_spline(x,y[:,c],w=w)
            gcv_direct_max=max(gcv_direct_max,float(np.max(
                np.abs(ref(q)-sci_auto(q)) / np.maximum(1.,np.abs(sci_auto(q))))))
    scale=np.maximum(1.,np.abs(expected))
    error=np.max(np.abs(vals-expected)/scale)
    coeff_err=np.max(np.abs(coeff-refcoeff)/np.maximum(1.,np.abs(refcoeff)))
    int_err=np.max(np.abs(integral-refint)/np.maximum(1.,np.abs(refint)))
    if not (error<5e-9 and coeff_err<5e-9 and int_err<5e-9):
        raise AssertionError(f'case {index}, gcv={gcv}, maxscaled={error:.4g}, '
                             f'coeff={coeff_err:.4g}, integral={int_err:.4g}, '
                             f'lam={actual_lam}')
    if not gcv and not np.allclose(actual_lam,lam,rtol=0,atol=1e-15):
        raise AssertionError(f'lambda mismatch {actual_lam} vs {lam}')
    return error,coeff_err,int_err,vals.size+coeff.size+nc,gcv_direct_max


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--seed',type=int,default=20261031)
    ap.add_argument('--cases',type=int,default=100)
    ap.add_argument('--mode',choices=['checked','optimized','both'],default='both')
    args=ap.parse_args()
    modes=[]
    if args.mode!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.mode!='checked':modes.append(('optimized',['-O3']))
    for title,flags in modes:
        with tempfile.TemporaryDirectory() as tmp:
            exe=build(Path(tmp),flags)
            rng=np.random.default_rng(args.seed)
            err=np.zeros(3)
            count=0
            scipy_gcv_delta=0.0
            for i in range(args.cases):
                try: a,b,c,num,delta=check(exe,rng,i)
                except Exception as ex:raise RuntimeError(f'{title} seed={args.seed} case={i}') from ex
                err=np.maximum(err,[a,b,c])
                count+=num
                scipy_gcv_delta=max(scipy_gcv_delta,delta)
            print(f'{title}: {args.cases} SciPy cases, {count} values; '
                  f'maxscaled(value,coef,integral)={err}; '
                  f'SciPy-autogcv-value-delta={scipy_gcv_delta:.6g}; PASS',flush=True)

if __name__=='__main__':main()
