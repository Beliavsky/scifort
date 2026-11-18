#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent ARGUS, Pearson III, and kappa4 censored references."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess
import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('argus', 'pearson3', 'kappa4')


def gamma_mass(y, width, shape):
    if y == 0:
        return mp.gammainc(shape,0,width)/mp.gamma(shape)
    if shape == 1:
        return mp.exp(-y)*(-mp.expm1(-width))
    if width*(1+abs(shape-1)/y) < mp.mpf('.001'):
        ref = (shape-1)*mp.log(y)-y-mp.loggamma(shape)
        return mp.exp(ref)*width*mp.quad(
            lambda t:mp.exp((shape-1)*mp.log1p(width*t/y)-width*t),
            [0,1],method='gauss-legendre',maxdegree=3)
    return mp.gammainc(shape,y,y+width)/mp.gamma(shape)


def normal_logs(z):
    p, q = mp.erfc(-z/mp.sqrt(2))/2, mp.erfc(z/mp.sqrt(2))/2
    return mp.log(p) if z <= 0 else mp.log1p(-q), mp.log(q) if z >= 0 else mp.log1p(-p)


def density(x, family, shapes, loc, scale):
    z = (mp.mpf(x)-loc)/scale
    a,b = map(mp.mpf,shapes)
    if family == 0:
        cap = a*a/2
        y = cap*(1-z*z)
        norm = mp.gammainc(mp.mpf('1.5'),0,cap)
        return 2*mp.log(a)+mp.log(z)+mp.log(y)/2-y-mp.log(norm)-mp.log(scale)
    if family == 1:
        if abs(a) < mp.mpf('1.6e-5'):
            return -z*z/2-mp.log(2*mp.pi)/2-mp.log(scale)
        beta=2/a
        alpha=beta*beta
        y=beta*(z+beta)
        return mp.log(abs(beta))+(alpha-1)*mp.log(y)-y-mp.loggamma(alpha)-mp.log(scale)
    lt = mp.log1p(-b*z)/b if b else -z
    lu = mp.log1p(-b*z) if b else 0
    if not a:
        return -mp.exp(lt)+lt-lu-mp.log(scale)
    lv = mp.log(-mp.expm1(mp.log(a)+lt)) if a>0 else mp.log1p(-a*mp.exp(lt))
    return (1/a-1)*lv+lt-lu-mp.log(scale)


def logs(x, family, shapes, loc, scale):
    z = (mp.mpf(x)-loc)/scale
    a,b = map(mp.mpf,shapes)
    if family == 0:
        cap = a*a/2
        y = cap*(1-z*z)
        norm = mp.gammainc(mp.mpf('1.5'),0,cap)/mp.gamma(mp.mpf('1.5'))
        p = gamma_mass(y,cap*z*z,mp.mpf('1.5'))/norm
        q = mp.gammainc(mp.mpf('1.5'),0,y)/mp.gamma(mp.mpf('1.5'))/norm
    elif family == 1:
        if abs(a) < mp.mpf('1.6e-5'):
            lp,lq = normal_logs(z)
            return density(x,family,shapes,loc,scale),lp,lq
        beta=2/a
        alpha=beta*beta
        y=beta*(z+beta)
        p=mp.gammainc(alpha,0,y)/mp.gamma(alpha)
        q=mp.gammainc(alpha,y,mp.inf)/mp.gamma(alpha)
        if beta < 0:
            p,q=q,p
    else:
        lt = mp.log1p(-b*z)/b if b else -z
        lv = mp.log(-mp.expm1(mp.log(a)+lt)) if a>0 else mp.log1p(-a*mp.exp(lt))
        lp = lv/a if a else -mp.exp(lt)
        p,q=mp.exp(lp),-mp.expm1(lp)
    lp = mp.log(p) if p < mp.mpf('.5') else mp.log1p(-q)
    lq = mp.log(q) if q < mp.mpf('.5') else mp.log1p(-p)
    return density(x,family,shapes,loc,scale),lp,lq


def interval(lo, hi, family, shapes, loc, scale):
    lo,hi=mp.mpf(lo),mp.mpf(hi)
    za,zb=(lo-loc)/scale,(hi-loc)/scale
    a,b=map(mp.mpf,shapes)
    if family == 0:
        cap=a*a/2
        norm=mp.gammainc(mp.mpf('1.5'),0,cap)/mp.gamma(mp.mpf('1.5'))
        return mp.log(gamma_mass(cap*(1-zb*zb),cap*(zb+za)*(hi-lo)/scale,mp.mpf('1.5'))/norm)
    if family == 1 and abs(a) >= mp.mpf('1.6e-5'):
        beta=2/a
        y1,y2=beta*(za+beta),beta*(zb+beta)
        return mp.log(gamma_mass(min(y1,y2),abs(beta)*(hi-lo)/scale,beta*beta))
    if (hi-lo)/scale < mp.mpf('1e-6'):
        ref=density(lo,family,shapes,loc,scale)
        return ref+mp.log(hi-lo)+mp.log(mp.quad(
            lambda t:mp.exp(density(lo+t*(hi-lo),family,shapes,loc,scale)-ref),
            [0,1],method='gauss-legendre',maxdegree=3))
    _,lp,lq=logs(lo,family,shapes,loc,scale)
    _,up,uq=logs(hi,family,shapes,loc,scale)
    return up+mp.log(-mp.expm1(lp-up)) if up<lq else lq+mp.log(-mp.expm1(uq-lq))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='gfortran')
    parser.add_argument('--module-dir',type=Path,required=True)
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--cases',type=int,default=100)
    args=parser.parse_args()
    rng=random.Random(20261117)
    cases,inputs=[],[]

    def add(family,shapes,loc=0.,scale=1.,exact=(),left=(),right=(),intervals=(),label='mixed'):
        loc,scale=mp.mpf(loc),mp.mpf(scale)
        value=sum(logs(x,family,shapes,loc,scale)[0] for x in exact)
        value+=sum(logs(x,family,shapes,loc,scale)[1] for x in left)
        value+=sum(logs(x,family,shapes,loc,scale)[2] for x in right)
        value+=sum(interval(a,b,family,shapes,loc,scale) for a,b in intervals)
        expected=float(value)
        assert math.isfinite(expected),(family,label,value)
        cases.append((FAMILIES[family],label,expected))
        inputs.append(f'{84+family} {shapes[0]:.17g} {float(loc):.17g} {float(scale):.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        if family==2:
            inputs.append(f'{shapes[1]:.17g}')
        for group in (exact,left,right,*intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    with mp.workdps(100):
        for family in range(3):
            for _ in range(args.cases):
                shapes=((10**rng.uniform(-2,1),0.) if family==0 else
                        (rng.choice((-1,1))*rng.uniform(.2,2),0.) if family==1 else
                        (rng.choice((-1,0,1))*rng.uniform(.1,.8),
                         rng.choice((-1,0,1))*rng.uniform(.05,.3)))
                scale=10**rng.uniform(-200,200)
                loc=scale*rng.uniform(-.3,.3)
                point=lambda x:loc+scale*x
                lo=point(rng.uniform(.2,.6))
                hi=lo+scale*10**rng.uniform(-15,-2)
                if hi==lo:
                    hi=math.nextafter(lo,math.inf)
                add(family,shapes,loc,scale,exact=[point(.4)],left=[point(.1)],
                    right=[point(.8)],intervals=[(lo,hi),(point(.2),point(.7))])
            print(f'{FAMILIES[family]}: mixed references ready',flush=True)
    with mp.workdps(300):
        for chi in (1e-6,.1,1.,10.,100.):
            for x in (1e-100,.5,math.nextafter(1.,0.)):
                add(0,(chi,0.),exact=[x],left=[x],right=[x],label='ARGUS shape/tail')
                hi=math.nextafter(x,math.inf)
                if hi<1:
                    add(0,(chi,0.),intervals=[(x,hi)],label='ARGUS adjacent')
        print('ARGUS stress references ready',flush=True)
    with mp.workdps(100):
        for skew in (-2000.,-2.,-.002,0.,1e-6,.002,2.,2000.):
            add(1,(skew,0.),exact=[0.],left=[0.],right=[0.],label='Pearson shape range')
            x=.0001
            add(1,(skew,0.),intervals=[(x,math.nextafter(x,math.inf))],label='Pearson adjacent')
        print('Pearson shape references ready',flush=True)
    with mp.workdps(100):
        for skew in (-2.,2.):
            x=math.copysign(1e100,skew)
            add(1,(skew,0.),exact=[x],left=[x],right=[x],label='Pearson exponential tail')
            add(1,(skew,0.),intervals=[(x,math.nextafter(x,math.inf))],label='Pearson steep adjacent')
        for h in (-10.,-1.,-1e-6,0.,1e-6,1.,10.):
            for k in (-10.,-1.,0.,1.,10.):
                if h>0 and k*math.log(h)>math.log(10):
                    continue
                hm,km,p=mp.mpf(h),mp.mpf(k),mp.mpf('.4')
                t=-mp.expm1(hm*mp.log(p))/hm if h else -mp.log(p)
                x=float(-mp.expm1(km*mp.log(t))/km if k else -mp.log(t)) if h>0 else .01
                add(2,(h,k),exact=[x],left=[x],right=[x],label='kappa shapes')
                add(2,(h,k),intervals=[(x,math.nextafter(x,math.inf))],label='kappa adjacent')
        for family,shapes in ((0,(2.,0.)),(1,(-1.,0.)),(2,(-.3,.2))):
            for scale in (1e-300,1e300):
                x=.4*scale
                add(family,shapes,scale=scale,exact=[x],left=[x],right=[x],label='extreme scale')
                add(family,shapes,scale=scale,intervals=[(x,math.nextafter(x,math.inf))],label='adjacent scale')
        for x in (1e-100,.4):
            add(2,(1.,1.),exact=[x],left=[x],right=[x],label='kappa uniform limit')
            add(2,(1.,1.),intervals=[(x,math.nextafter(x,math.inf))],label='uniform adjacent')
        for x in (-100.,100.,700.):
            add(2,(-1.,0.),exact=[x],left=[x],right=[x],label='kappa logistic tail')
            add(2,(-1.,0.),intervals=[(x,math.nextafter(x,math.inf))],label='logistic adjacent')
    cache = ROOT/'build'/'gamma_kappa_reference'
    cache.mkdir(parents=True, exist_ok=True)
    (cache/'cases.json').write_text(json.dumps(cases), encoding='ascii')
    (cache/'input.txt').write_text('\n'.join(inputs)+'\n', encoding='ascii')
    for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow']),
                        ('optimized', ['-O3', '-fcheck=bounds'])]:
        out = cache/mode
        out.mkdir(parents=True, exist_ok=True)
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
            n = len(result.stdout.split())
            raise RuntimeError((n, cases[n:n+1], result.stderr))
        actual = [float(s) for s in result.stdout.split()]
        assert len(actual) == len(cases)
        worst = 0.
        for value, (family, label, expected) in zip(actual, cases):
            error = abs(value-expected)/max(1., abs(expected))
            assert math.isfinite(value) and error <= 2e-11, \
                (mode, family, label, value, expected, error)
            worst = max(worst, error)
        print(f'{mode}: {len(cases)} comparisons passed; max scaled error={worst}', flush=True)
    sample=stats.CensoredData(uncensored=[.2,.4,.6],left=[.1],right=[.8],
                              interval=[[.3,.35],[.5,.55]])
    for name,shapes in zip(FAMILIES,((2.,),(1.,),(.3,.2))):
        fit=optimize.minimize_scalar(
            lambda v:getattr(stats,name)._penalized_nnlf((*shapes,0.,v),sample),
            bounds=(.3,3.),method='bounded',options={'xatol':1e-13})
        assert fit.success
        print(f'{name}: SciPy scale={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__=='__main__':
    main()
