#!/usr/bin/env python3
"""Compare dependency-free FITPACK TCK wrappers with SciPy/FITPACK.

This checker covers s=0 and task=-1 fixed knots. Positive-s native smoothing
is checked separately by validation/check_splrep_smoothing.py.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.interpolate import splrep, splev, splint, splder, splantider, spalde

ROOT=Path(__file__).resolve().parents[2]
SOURCES=('kinds','interpolate_bspline','interpolate_fitpack_bivariate','interpolate_curfit_core', 'interpolate_fitpack_compat')

def run(command:list[str],cwd:Path,source:str|None=None)->str:
    p=subprocess.run(command, cwd=cwd, input=source,text=True,capture_output=True)
    if p.returncode:
        raise RuntimeError(f"{' '.join(map(str,command))}\nstdout:\n{p.stdout[-1000:]}\nstderr:\n{p.stderr[-3000:]}")
    return p.stdout

def build(work:Path,flags:list[str])->Path:
    common=['gfortran','-std=f2018','-ffree-line-length-none',*flags]
    for module in SOURCES:
        run([*common, '-I'+str(work),'-J'+str(work),'-c',str(ROOT/'src'/f'{module}.f90'),
             '-o',str(work/f'{module}.o')],work)
    exe=work/'fitpack_tck_driver'
    run([*common,'-I'+str(work),str(ROOT/'tools/interpolate_reference/fitpack_tck_driver.f90'),
         *(str(work/f'{module}.o') for module in SOURCES),'-o',str(exe)],work)
    return exe

def case(exe:Path,rng:np.random.Generator,i:int)->tuple[int,float]:
    k=1+i%5
    n=14+(i*13)%65
    x=np.cumsum(rng.uniform(.08,1.0,size=n))
    x=(x-x[0])/(x[-1]-x[0])*6-2
    w=np.exp(rng.uniform(-0.8,0.8,n))
    y=np.sin((i%6+1)*x*.4)+.2*np.cos(x*1.7)+.05*rng.normal(size=n)
    if i%9==0:
        y=.7+1.2*x-.5*x**2
    task=-1 if i%2 else 0
    nt=min(6,(i//5)%7)
    if task==0:nt=0
    ids=np.linspace(k+1,n-k-2,nt,dtype=int) if nt else np.empty(0,dtype=int)
    t=x[ids] if nt else np.array([],dtype=float)
    bbox=[x[0],x[-1]]
    if task==-1 and i%7==0:
        bbox[0]-=.1; bbox[1]+=.2
    q=np.r_[x[0]-.3,np.linspace(x[0],x[-1],21),x[-1]+.3]
    values=[f'{n} {len(q)} {k} {task} {nt}',
            ' '.join(f'{z:.17g}' for z in bbox),
            *(' '.join(f'{z:.17g}' for z in a) for a in (x,y,w))]
    if nt: values.append(' '.join(f'{z:.17g}' for z in t))
    values.append(' '.join(f'{z:.17g}' for z in q))
    output=run([str(exe)],exe.parent,'\n'.join(values)+'\n').splitlines()
    pos=0
    def pop()->str:
        nonlocal pos
        result=output[pos];pos+=1
        return result
    status=int(pop())
    if status: raise AssertionError((i,'Fortran status',status))
    nkn,ncoef,degree=map(int,pop().split())
    knots=np.fromstring(pop(),sep=' ')
    coeff=np.fromstring(pop(),sep=' ')
    assert nkn==ncoef==len(knots)==len(coeff) and degree==k
    kwargs={'k':k,'w':w,'xb':bbox[0],'xe':bbox[1]}
    if task==-1: kwargs.update(task=-1,t=t)
    else: kwargs['s']=0
    ref=splrep(x,y,**kwargs)
    maxerr=0.
    def check(a,b,atol=2.e-8,rtol=1.e-9):
        nonlocal maxerr
        a=np.asarray(a); b=np.asarray(b)
        if a.shape!=b.shape:raise AssertionError((i,'shape',a.shape,b.shape))
        dif=np.max(np.abs(a-b)) if a.size else 0.
        maxerr=max(maxerr,float(dif))
        np.testing.assert_allclose(a,b,atol=atol,rtol=rtol,err_msg=f'case {i}')
    check(knots,ref[0],atol=2.e-12)
    check(coeff[:nkn-k-1],ref[1][:nkn-k-1],atol=5.e-9)
    for nu in range(min(k,3)+1):
        check(np.fromstring(pop(),sep=' '),splev(q,ref,der=nu))
    for ext in (1,2,3):
        status=int(pop())
        if ext==2:
            assert status==4, (i,'out of bounds status',status)
        else:
            assert status==0
            check(np.fromstring(pop(),sep=' '),splev(q,ref,ext=ext))
    check(float(pop()),splint(q[0],q[-1],ref))
    assert int(pop())==0
    check(np.fromstring(pop(),sep=' '),spalde(.5*(x[1]+x[2]),ref))
    assert int(pop())==0
    check(np.fromstring(pop(),sep=' '),splev(q,splder(ref)))
    assert int(pop())==0
    check(np.fromstring(pop(),sep=' '),splev(q,splantider(ref)))
    return n,maxerr

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',type=int,default=24680)
    parser.add_argument('--cases',type=int,default=90)
    parser.add_argument('--only',choices=['checked','optimized','both'],default='both')
    args=parser.parse_args()
    modes=[]
    if args.only!='optimized':modes.append(('checked',['-O0','-g','-fcheck=all','-fbacktrace']))
    if args.only!='checked':modes.append(('optimized',['-O3']))
    for label,flags in modes:
        with tempfile.TemporaryDirectory() as tmp:
            work=Path(tmp)
            exe=build(work,flags)
            rng=np.random.default_rng(args.seed)
            error=0.
            for i in range(args.cases):
                _,e=case(exe,rng,i)
                error=max(error,e)
            print(f'{label}: {args.cases}/{args.cases} SciPy TCK cases PASS; max absolute difference {error:.3g}',flush=True)

if __name__=='__main__':main()
