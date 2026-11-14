#!/usr/bin/env python3
"""Independent scipy.interpolate.splprep(s=0/task=-1) reference checks.

No perazz/FITPACK dependency is needed: only SciPy is used as an oracle.
The SciFort API uses points(nsample, ndim), whereas SciPy uses a sequence of
coordinate arrays. Tests cover automatic chord-length u, explicit u, weighted
fixed-knot least squares, scalar/vector derivatives, coefficients and knots.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import splprep, splev

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('kinds', 'interpolate_bspline', 'interpolate_curfit_core', 'interpolate_fitpack_compat')

def build(work: Path, flags: list[str]) -> Path:
    base=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    for module in SOURCES:
        cmd=[*base,'-I'+str(work),'-J'+str(work),'-c',str(ROOT/'src'/f'{module}.f90'),
             '-o',str(work/f'{module}.o')]
        p=subprocess.run(cmd,cwd=work,capture_output=True,text=True)
        if p.returncode: raise RuntimeError(f'{cmd}:\n{p.stderr[-2500:]}')
    exe=work/'fitpack_parametric_driver'
    cmd=[*base,'-I'+str(work),str(ROOT/'tools/interpolate_reference/fitpack_parametric_driver.f90'),
         *(str(work/f'{module}.o') for module in SOURCES),'-o',str(exe)]
    p=subprocess.run(cmd,cwd=work,capture_output=True,text=True)
    if p.returncode: raise RuntimeError(f'{cmd}:\n{p.stderr[-2500:]}')
    return exe

def trial(exe: Path, rng: np.random.Generator, i: int) -> tuple[float,float,int]:
    k=1+i%5
    n=17+(i*13)%37
    d=2+i%5
    raw=rng.normal(size=(n,d))
    x=np.cumsum(raw*.08 + .1*np.linspace(.8,1.2,d),axis=0)
    x[:,0]+=np.sin(np.linspace(0,3,n))
    x[:,1]+=np.cos(np.linspace(0,5,n))
    weights=np.exp(rng.uniform(-.7,.7,n))
    has_u=i%3!=0
    if has_u:
        dif=rng.uniform(.06,1.5,n-1)
        u=np.r_[0.,np.cumsum(dif)]
        u=(u/u[-1])*4.0-2.0
    else:
        dif=np.linalg.norm(np.diff(x,axis=0),axis=1)
        u=np.r_[0.,np.cumsum(dif)]
        u=u/u[-1]
    mode=-1 if (i%4==0) else 0
    nt=0
    full_t=np.empty(0)
    if mode==-1:
        nt=1+(i//4)%4
        interior=np.linspace(u[0],u[-1],nt+2)[1:-1]
        full_t=np.r_[np.repeat(u[0],k+1),interior,np.repeat(u[-1],k+1)]
    q=np.r_[u[0]-.1*(u[-1]-u[0]),np.linspace(u[0],u[-1],12),u[-1]+.1*(u[-1]-u[0])]
    parts=[f'{n} {d} {len(q)} {k} {mode} {len(full_t)} {int(has_u)}']
    parts.extend(' '.join(f'{v:.17g}' for v in row) for row in x)
    parts.append(' '.join(f'{v:.17g}' for v in weights))
    if has_u:parts.append(' '.join(f'{v:.17g}' for v in u))
    if mode==-1:parts.append(' '.join(f'{v:.17g}' for v in full_t))
    parts.append(' '.join(f'{v:.17g}' for v in q))
    proc=subprocess.run([str(exe)],cwd=exe.parent,input='\n'.join(parts)+'\n',text=True,capture_output=True)
    if proc.returncode:raise AssertionError(f'case {i}, Fortran crash: {proc.stderr[-1000:]}')
    lines=iter(proc.stdout.splitlines())
    status=int(next(lines))
    if status: raise AssertionError(f'case {i}, status={status}')
    nkn,ncoeff,degree=map(int,next(lines).split())
    actual_u=np.fromstring(next(lines),sep=' ')
    actual_t=np.fromstring(next(lines),sep=' ')
    actual_c=np.array([np.fromstring(next(lines),sep=' ') for _ in range(d)])
    kwargs=dict(k=k,task=mode,w=weights,s=0.0)
    if has_u:kwargs['u']=u
    if mode==-1:kwargs['t']=full_t
    reference,ref_u=splprep(x.T,**kwargs)
    assert nkn==len(reference[0]) and ncoeff==len(reference[1][0]) and degree==k
    worst=0.0
    worst_scaled=0.0
    def check(a,b,tol=2.e-9):
        nonlocal worst, worst_scaled
        a=np.asarray(a);b=np.asarray(b)
        if a.shape!=b.shape:raise AssertionError((i,'shape',a.shape,b.shape))
        err=np.abs(a-b)
        worst=max(worst,float(np.max(err)) if err.size else 0.0)
        if err.size: worst_scaled=max(worst_scaled,float(np.max(err/(1.0+np.abs(b)))))
        np.testing.assert_allclose(a,b,atol=tol,rtol=2.e-9,
                                   err_msg=f'case={i}, k={k}, mode={mode}')
    check(actual_u,ref_u,2.e-12)
    check(actual_t,reference[0],2.e-12)
    check(actual_c,reference[1],5.e-9)
    for der in range(min(k,3)+1):
        assert int(next(lines))==0
        actual=np.array([np.fromstring(next(lines),sep=' ') for _ in q])
        expected=np.asarray(splev(q,reference,der=der)).T
        check(actual,expected,5.e-8)
    assert int(next(lines))==0
    check(np.fromstring(next(lines),sep=' '),np.array(splev(q[0],reference,der=1)),5.e-8)
    return worst,worst_scaled,d*n

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=20261101)
    parser.add_argument('--cases',type=int,default=90)
    parser.add_argument('--only',choices=('checked','optimized','both'),default='both')
    args=parser.parse_args()
    modes=[]
    if args.only!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.only!='checked':modes.append(('optimized',['-O3']))
    for name,flags in modes:
        with tempfile.TemporaryDirectory() as temp:
            exe=build(Path(temp),flags)
            rng=np.random.default_rng(args.seed)
            worst=0.0
            worst_scaled=0.0
            total=0
            for i in range(args.cases):
                err,scaled,num=trial(exe,rng,i)
                worst=max(worst,err);worst_scaled=max(worst_scaled,scaled);total+=num
            print(f'{name}: {args.cases}/{args.cases} SciPy parametric cases PASS; '
                  f'{total} samples x channels; max abs error {worst:.3g}; '
                  f'max scaled error {worst_scaled:.3g}',flush=True)
if __name__=='__main__':main()
