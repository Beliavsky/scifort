#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-precision generalized normal, half generalized normal, and R distribution checks."""
import argparse
from functools import lru_cache
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('gennorm', 'halfgennorm', 'rdist')


def softplus(q):
    return max(q, 0)+mp.log1p(mp.exp(-abs(q)))


def logadd(a, b):
    if a == -mp.inf: return b
    if b == -mp.inf: return a
    return max(a,b)+mp.log1p(mp.exp(-abs(a-b)))


def difference(a,b):
    if a == b: return -mp.inf
    return a+mp.log(-mp.expm1(b-a))


def half_logs(t,beta,s):
    a=1/beta
    if not t: return -mp.loggamma(1+a)-mp.log(s),-mp.inf,mp.mpf(0)
    ly=beta*mp.log(t/s)
    if ly > mp.log(mp.mpf('1e310')):
        return -mp.inf,mp.mpf(0),-mp.inf
    y=mp.exp(ly)
    density=-mp.loggamma(1+a)-mp.log(s)-y
    if y < a/10:
        term=mp.mpf(1)
        total=term
        n=0
        while True:
            n+=1
            term*=y/(a+n)
            total+=term
            if abs(term) < mp.eps*abs(total):break
        lp=a*ly-y-mp.loggamma(a+1)+mp.log(total)
        lq=mp.log1p(-mp.exp(lp))
    else:
        lp=mp.log(mp.gammainc(a,0,y))-mp.loggamma(a)
        lq=mp.log(mp.gammainc(a,y,mp.inf))-mp.loggamma(a)
    return density,lp,lq


def logs(x,family,shape,second,loc,scale):
    return cached_logs(x,family,shape,second,loc,scale,mp.mp.dps)


@lru_cache(maxsize=None)
def cached_logs(x,family,shape,second,loc,scale,precision):
    x,c,mu,s=map(mp.mpf,(x,shape,loc,scale))
    t=x-mu
    if family<2:
        if family==1 and t<0: return -mp.inf,-mp.inf,mp.mpf(0)
        density,lp,lq=half_logs(abs(t),c,s)
        if family==0:
            density-=mp.log(2)
            small=lq-mp.log(2)
            if t<0: lp,lq=small,mp.log1p(-mp.exp(small))
            else: lp,lq=mp.log1p(-mp.exp(small)),small
        return density,lp,lq
    lo=mp.mpf(float(loc-scale))
    width=mp.mpf(float(2*scale))
    a=mp.mpf(float(shape/2))
    u=(x-lo)/width
    if u<0: return -mp.inf,-mp.inf,mp.mpf(0)
    if u>1: return -mp.inf,mp.mpf(0),-mp.inf
    density=(a-1)*(mp.log(u)+mp.log1p(-u))-mp.log(mp.beta(a,a))-mp.log(width)
    if u==mp.mpf('.5'): return density,-mp.log(2),-mp.log(2)
    if a>=5000:
        # Direct density integral in concentration-scaled coordinates.
        z=2*u-1
        root=mp.sqrt(2*a)
        q=abs(z)*root
        normalization=mp.beta(mp.mpf('.5'),a)*root
        def integrand(t):
            if t>=root:return mp.mpf(0)
            return mp.exp((a-1)*mp.log1p(-t*t/(2*a)))/normalization
        tail=mp.quad(integrand,[q,q+1,q+4,max(q+10,20),root])
        small=mp.log(tail)
        other=mp.log1p(-tail)
        return (density,other,small) if z>0 else (density,small,other)
    lp=mp.log(mp.betainc(a,a,0,u,regularized=True))
    lq=mp.log(mp.betainc(a,a,u,1,regularized=True))
    return density,lp,lq


def interval(a,b,family,shape,second,loc,scale):
    if family < 2:
        lo,hi,c,mu,s=map(mp.mpf,(a,b,shape,loc,scale))
        lo-=mu
        hi-=mu
        if family==0 and hi<=0:lo,hi=-hi,-lo
        if family==0 and lo<0:
            return logadd(half_logs(-lo,c,s)[1],half_logs(hi,c,s)[1])-mp.log(2)
        lo=max(lo,mp.mpf(0))
        _,pa,qa=half_logs(lo,c,s)
        _,pb,qb=half_logs(hi,c,s)
        value=difference(pb,pa) if pb<=qa else difference(qa,qb)
        return value-(mp.log(2) if family==0 else 0)
    _,pa,qa=logs(a,family,shape,second,loc,scale)
    _,pb,qb=logs(b,family,shape,second,loc,scale)
    if pb <= qa: return difference(pb,pa)
    if qa <= -mp.log(2): return difference(qa,qb)
    return mp.log1p(-mp.exp(pa)-mp.exp(qb))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module-dir', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--cases', type=int, default=100)
    parser.add_argument('--seed', type=int, default=20261117)
    args = parser.parse_args()
    mp.mp.dps = 100
    rng = random.Random(args.seed)
    cases, inputs = [], []

    def add(family, shape, second=1.2, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed'):
        value = sum(logs(x, family, shape, second, loc, scale)[0] for x in exact)
        value += sum(logs(x, family, shape, second, loc, scale)[1] for x in left)
        value += sum(logs(x, family, shape, second, loc, scale)[2] for x in right)
        value += sum(interval(a, b, family, shape, second, loc, scale) for a, b in intervals)
        expected = float(value)
        if not math.isfinite(expected):
            raise AssertionError((FAMILIES[family], label, value))
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{66+family} {shape:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    for family in range(3):
        for _ in range(args.cases):
            shape=10**rng.uniform(-.5,1.1) if family<2 else 10**rng.uniform(-1,2)
            scale=10**rng.uniform(-100,100)
            loc=scale*rng.uniform(-3,3)
            point=lambda p:loc+scale*p
            low=-.85 if family!=1 else .01
            exact=[point(rng.uniform(low,.85)) for _ in range(3)]
            left=[point(rng.uniform(low,.2))]
            right=[point(rng.uniform(.3,.9))]
            a=point(rng.uniform(max(low,.01),.7))
            b=a+scale*10**rng.uniform(-14,-2)
            if b==a:b=math.nextafter(a,math.inf)
            add(family,shape,loc=loc,scale=scale,exact=exact,left=left,right=right,
                intervals=[(a,b),(point(.1),point(.8))])
    with mp.workdps(1100):
        tiny=math.nextafter(0.,math.inf)
        for family in range(3):
            shapes=(1e-300,.1,1.,2.,100.,1e100,1e300) if family<2 else (2*tiny,1e-300,.1,1.,2.,3.,100.)
            for shape in shapes:
                x=1. if family<2 else .2
                add(family,shape,exact=[x],left=[x],right=[x],label='shape stress')
                add(family,shape,intervals=[(x,math.nextafter(x,math.inf))],label='adjacent shape')
            for scale in (1e-300,1e300):
                x=.2*scale
                add(family,2.,scale=scale,exact=[x],left=[x],right=[x],label='scale stress')
                add(family,2.,scale=scale,intervals=[(x,math.nextafter(x,math.inf))],label='adjacent scale')
        for family in (0,1):
            largest=float.fromhex('0x1.fffffffffffffp+1023')
            add(family,largest,exact=[.1],left=[.1],right=[.1],label='largest finite power')
            add(family,largest,exact=[1.],left=[1.],right=[1.],label='largest power boundary')
            add(family,largest,intervals=[(.1,.2)],label='largest finite power interval')
            add(family,1e6,intervals=[(math.exp(math.log(.0003)/1e6),
                math.exp(math.log(.001)/1e6))],label='steep small-power interval')
            add(family,1.,scale=1e300,intervals=[(0.,tiny)],label='subnormal interval')
            add(family,1e300,loc=-1.,left=[1e-300],label='amplified scale point')
            add(family,1e300,loc=-1.,intervals=[(0.,1e-300)],label='amplified scale interval')
            add(family,1e300,exact=[0.],label='nearly uniform endpoint log density')
        add(0,2.,intervals=[(-1e-300,1e-300)],label='tiny center crossing')
        add(0,2.,intervals=[(-1e100,1e100)],label='nearly full support')
        for shape in (1e4,1e6):
            add(2,shape,exact=[0.],left=[0.],right=[0.],label='large symmetric shape center')
        add(2,2.,scale=21*tiny,exact=[0.],intervals=[(0.,tiny)],label='subnormal R width')
    with mp.workdps(90):
        for x in (-.01,-.001,.001,.01):
            add(2,1e6,exact=[x],left=[x],right=[x],label='maximum R shape tails')
            add(2,1e6,intervals=[(x,math.nextafter(x,math.inf))],label='maximum R adjacent')
    for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow']),
                        ('optimized', ['-O3', '-fcheck=bounds'])]:
        out = ROOT/'build'/'generalized_normal_reference'/mode
        out.mkdir(parents=True, exist_ok=True)
        (out/'input.txt').write_text('\n'.join(inputs)+'\n',encoding='ascii')
        cmd = [args.compiler, *flags, '-J'+str(out), '-I'+str(out),
               '-I'+str(args.module_dir), 'src/special/log_gamma.f90',
               'src/empirical_distribution.f90', 'src/censored.f90',
               'validation/censored_reference_driver.f90', str(args.library),
               '-o', str(out/'driver.exe')]
        print(subprocess.list2cmdline(cmd), flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)
        result = subprocess.run([str(out/'driver.exe')], input='\n'.join(inputs)+'\n',
                                text=True, capture_output=True)
        if result.returncode:
            completed=len(result.stdout.split())
            raise RuntimeError((completed,cases[completed:completed+1],result.stderr))
        actual = [float(s) for s in result.stdout.split()]
        assert len(actual) == len(cases), (len(actual), len(cases))
        worst = (0., None)
        for value, (family, label, expected) in zip(actual, cases):
            error = abs(value-expected)/max(1., abs(expected))
            if not math.isfinite(value) or error > 2e-11:
                raise AssertionError((mode, family, label, value, expected, error))
            if error > worst[0]:
                worst = (error, (family, label))
        print(f'{mode}: {len(cases)} comparisons passed; max scaled error={worst}', flush=True)
    sample=stats.CensoredData(uncensored=[.2,.4,.6],left=[.1],right=[.7],interval=[[.25,.35],[.5,.55]])
    for name in FAMILIES:
        distribution=getattr(stats,name)
        fit=optimize.minimize_scalar(lambda c:distribution._penalized_nnlf((c,0.,1.3),sample),
            bounds=(.3,5.),method='bounded',options={'xatol':1e-13})
        assert fit.success
        print(f'{name}: SciPy shape={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__=='__main__':main()
