#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent nested-order-integral references for censored kstwo (mpmath)."""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import subprocess
import mpmath as mp
from scipy import stats, optimize
ROOT=Path(__file__).resolve().parents[1]

@lru_cache(None)
def volume(n,d):
    # Nested integrals over ordered uniforms. Global polynomial antiderivatives
    # provide an independent oracle, rather than count convolution/face flux.
    lo=[max(mp.mpf(0),mp.mpf(i)/n-d) for i in range(1,n+1)]
    hi=[min(mp.mpf(1),mp.mpf(i-1)/n+d) for i in range(1,n+1)]
    grid=sorted(set([mp.mpf(0),mp.mpf(1)]+lo+hi))
    polys=[[mp.mpf(1)] for _ in grid[:-1]]
    for lower,upper in zip(lo,hi):
        carry=mp.mpf(0);new=[]
        for a,b,poly in zip(grid,grid[1:],polys):
            if b<=lower: q=[mp.mpf(0)]
            elif a>=upper: q=[carry]
            else:
                q=[mp.mpf(0)]+[v/(j+1) for j,v in enumerate(poly)]
                q[0]=carry-mp.polyval(q[::-1],a)
                carry=mp.polyval(q[::-1],b)
            new.append(q)
        polys=new
    return mp.factorial(n)*mp.polyval(polys[-1][::-1],1)

def cdf(n,z):
    if z<=mp.mpf(1)/(2*n):return mp.mpf(0)
    if z>=1:return mp.mpf(1)
    if z<=mp.mpf(1)/n:return mp.factorial(n)*(2*z-mp.mpf(1)/n)**n
    if 1-z<=mp.mpf(1)/n:return 1-2*(1-z)**n
    return volume(n,z)

def logs(n,x,loc=0.,scale=1.):
    lower=loc+scale/(2*n);upper=loc+scale;knot=loc+scale/n
    x,loc,scale=map(mp.mpf,(x,loc,scale))
    if x<=lower:return -mp.inf,-mp.inf,mp.mpf(0)
    if x>=upper:return -mp.inf,mp.mpf(0),-mp.inf
    if x<knot:
        t=2*(x-mp.mpf(lower))/scale
        lp=mp.loggamma(n+1)+n*mp.log(t)
        ld=mp.log(2*n)+mp.loggamma(n+1)+(n-1)*mp.log(t)-mp.log(scale)
        return ld,lp,mp.log1p(-mp.exp(lp))
    if x==knot:
        lp=mp.loggamma(n+1)-n*mp.log(n)
        ld=mp.log(2*n)+mp.loggamma(n+1)-(n-1)*mp.log(n)-mp.log(scale)
        return ld,lp,mp.log1p(-mp.exp(lp))
    t=(mp.mpf(upper)-x)/scale
    if t<=mp.mpf(1)/n:
        lq=mp.log(2)+n*mp.log(t)
        return mp.log(2*n)+(n-1)*mp.log(t)-mp.log(scale),mp.log1p(-mp.exp(lq)),lq
    z=max((x-loc)/scale,mp.mpf(math.nextafter(1/n,math.inf)))
    probability=cdf(n,z)
    density=mp.diff(lambda y:cdf(n,y),z)
    return mp.log(density)-mp.log(scale),mp.log(probability),mp.log1p(-probability)

def interval(n,a,b,loc=0.,scale=1.):
    # High precision retains adjacent interval differences in polynomial pieces.
    lower=loc+scale/(2*n);upper=loc+scale;knot=loc+scale/n
    a=max(a,lower);b=min(b,upper)
    if b<=a:return -mp.inf
    if a<knot<b:
        return mp.log(mp.exp(interval(n,a,knot,loc,scale))+mp.exp(interval(n,knot,b,loc,scale)))
    if b<=knot:
        def power(x):return mp.factorial(n)*(2*(mp.mpf(x)-lower)/mp.mpf(scale))**n
        return mp.log(power(b)-power(a))
    if (mp.mpf(upper)-a)/scale<=mp.mpf(1)/n:
        return mp.log(2*((mp.mpf(upper)-a)/scale)**n-2*((mp.mpf(upper)-b)/scale)**n)
    za=(mp.mpf(a)-loc)/scale;zb=(mp.mpf(b)-loc)/scale
    # The stored affine knot defines the change of polynomial at binary64 level.
    if a==knot:
        za=mp.mpf(1)/n
        zb=za+(mp.mpf(b)-a)/scale
    return mp.log(cdf(n,zb)-cdf(n,za))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse',action='store_true');args=parser.parse_args()
    cache=ROOT/'build'/'kstwo-reference';cache.mkdir(exist_ok=True)
    cases=[];inputs=[]
    def add(n,x,part,loc=0.,scale=1.,end=None):
        expected=logs(n,x,loc,scale)[part] if part<3 else interval(n,x,end,loc,scale)
        cases.append(float(expected))
        inputs.append(f'107 {n} {loc:.17g} {scale:.17g} '+' '.join(str(int(part==j)) for j in range(4)))
        inputs.append(f'{x:.17g}'+(' '+f'{end:.17g}' if end is not None else ''))
    if args.reuse:
        cases=json.loads((cache/'cases.json').read_text());inputs=(cache/'input.txt').read_text().splitlines()
    else:
        with mp.workdps(100):
            for n in (1,2,3,4,10,32):
                for x in sorted(set([.5/n,math.nextafter(.5/n,math.inf),1/n,.475,.7,math.nextafter(1.,0.),1.])):
                    for part in range(3):add(n,x,part)
                a=max(.55,1.2/n);b=.85
                if a<b:add(n,a,3,end=b)
                if n>1:
                    a=max(.475,1.1/n)
                    add(n,a,3,end=math.nextafter(a,math.inf))
                print('MP fixtures ready n=',n,flush=True)
            for n in (3,10,32):
                for scale in (1e-200,1.,1e200,1e307):
                    for loc in (0.,-scale):
                        for part in range(3):add(n,loc+.475*scale,part,loc,scale)
                        add(n,loc+.475*scale,3,loc,scale,loc+.7*scale)
                for loc in (-100.,100.):
                    knot=loc+1/n
                    add(n,knot,3,loc,1.,math.nextafter(knot,math.inf))
                    add(n,math.nextafter(knot,-math.inf),3,loc,1.,math.nextafter(knot,math.inf))
        (cache/'cases.json').write_text(json.dumps(cases));(cache/'input.txt').write_text('\n'.join(inputs)+'\n')
    for name,opt in (('checked','-O0'),('optimized','-O3')):
        out=cache/name;out.mkdir(exist_ok=True)
        cmd=['gfortran',opt,'-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-J'+str(out),'-I'+str(out),'-Ibuild/gfortran_2806A042DAA8BBF7','src/censored.f90','validation/censored_reference_driver.f90','build/gfortran_3C159EB802C22FB9/scifort/libscifort.a','-o',str(out/'driver.exe')]
        print(subprocess.list2cmdline(cmd),flush=True);subprocess.run(cmd,cwd=ROOT,check=True)
        result=subprocess.run([str(out/'driver.exe')],input='\n'.join(inputs)+'\n0\n',text=True,capture_output=True,check=True)
        actual=list(map(float,result.stdout.split()));assert len(actual)==len(cases)
        errors=[0. if a==b else (abs(a/b-1) if 0<abs(b)<1e-10 else
                 abs(a-b)/max(1.,abs(b))) for a,b in zip(actual,cases)]
        bad=[(i,a,b,e) for i,(a,b,e) in enumerate(zip(actual,cases,errors)) if not math.isfinite(e) or e>3e-11]
        assert not bad,(name,bad)
        print(name,len(cases),'MP cases passed; maximum scaled log error',max(errors),flush=True)
    sample=stats.CensoredData(uncensored=[.3,.4,.6],left=[.25],right=[.8],interval=[[.3,.35],[.5,.55]])
    fit=optimize.minimize_scalar(lambda s:stats.kstwo._penalized_nnlf((4.,0.,s),sample),bounds=(.9,1.5),method='bounded',options={'xatol':1e-13})
    print('SciPy smooth local scale fit:',fit.x,'NNLF',fit.fun,flush=True)
    knot_objective=stats.kstwo._penalized_nnlf((4.,0.,1.2),sample)
    assert knot_objective<fit.fun
    print('SciPy first-knot scale 1.2 NNLF:',knot_objective,flush=True)
if __name__=='__main__':main()
