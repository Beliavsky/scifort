#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent scaled adaptive references for censored studentized range, k<=100.

Run after validation/check_censored_final.py to build checked/optimized drivers.
Use --reuse to repeat comparisons without regenerating the reference fixtures.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import subprocess
import mpmath as mp
import numpy as np
from scipy import integrate, optimize, special, stats
ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT/'build/final-censored-reference'
LOGPHI = math.log(2*math.pi)/2
NODES, WEIGHTS = special.roots_legendre(16)


def logdiff(a,b):
    if b == -math.inf:
        return a
    if b >= a:
        return -math.inf
    return a+math.log(-math.expm1(b-a))


@lru_cache(maxsize=100000)
def normal_range(r,k,part):
    """Positive minimum-normal integrals, scaled before adaptive quadrature."""
    center = 0. if part == 1 else -r/2
    def kernel(u):
        z = center+u
        lp = -z*z/2-LOGPHI
        if part == 2:
            lq = special.log_ndtr(-z)
            ratio = min(0.,special.log_ndtr(-z-r)-lq)
            if ratio == 0.:
                lm = 0.
            elif ratio < -36:
                lm = math.log(k-1)+ratio
            else:
                ld = math.log1p(-math.exp(ratio)) if ratio < -math.log(2) else math.log(-math.expm1(ratio))
                lm = math.log(-math.expm1((k-1)*ld))
            return math.log(k)+lp+(k-1)*lq+lm
        if r <= 1.:
            # Independently integrate the Gaussian over the small interval.
            # Factoring phi(midpoint) avoids subtracting nearly equal CDFs.
            mid = z+r/2
            factor = np.dot(WEIGHTS,np.exp(-mid*r*NODES/2-r*r*NODES*NODES/8))/2
            ld = math.log(r)-mid*mid/2-LOGPHI+math.log(factor)
        elif z >= 0:
            ld = logdiff(special.log_ndtr(-z),special.log_ndtr(-z-r))
        else:
            ld = logdiff(special.log_ndtr(z+r),special.log_ndtr(z))
        if part == 0:
            return math.log(k*(k-1))+lp-(z+r)**2/2-LOGPHI+(k-2)*ld
        return math.log(k)+lp+(k-1)*ld
    peak = optimize.minimize_scalar(lambda u:-kernel(u),bounds=(-12.,12.),method='bounded').x
    shift = kernel(peak)
    value = integrate.quad(lambda u:math.exp(kernel(u)-shift),-12,12,
                           points=[peak],epsabs=2e-13,epsrel=2e-13,limit=300)[0]
    return shift+math.log(value)


def finite_range(q,k,df,part):
    if math.isinf(df):
        return normal_range(q,k,part)
    c = math.log(2)+df/2*math.log(df/2)-special.gammaln(df/2)
    # Direct scale integration differs from the runtime fixed log-scale rule.
    high = math.sqrt((df+k+400)/df)
    def kernel(s):
        if s <= 0:
            return -math.inf
        return c+(df-1)*math.log(s)-df*s*s/2+normal_range(q*s,k,part)+(math.log(s) if part == 0 else 0)
    peak = optimize.minimize_scalar(lambda s:-kernel(s),bounds=(1e-12,high),method='bounded').x
    shift = kernel(peak)
    value = integrate.quad(lambda s:math.exp(kernel(s)-shift),0,high,
                           points=[peak],epsabs=8e-12,epsrel=8e-12,limit=300)[0]
    return shift+math.log(value)


def distribution_logs(q,k,df):
    ld = finite_range(q,k,df,0)
    lc = finite_range(q,k,df,1)
    if lc < -math.log(2):
        return ld,lc,math.log1p(-math.exp(lc))
    ls = finite_range(q,k,df,2)
    return ld,math.log1p(-math.exp(ls)),ls


def lower_constant(k,df):
    k,df = map(mp.mpf,(k,df))
    m = k-1
    c = mp.log(k)/2-m*mp.log(2*mp.pi)/2
    if not mp.isinf(df):
        c += m/2*mp.log(2/df)+mp.loggamma((df+m)/2)-mp.loggamma(df/2)
    return c


def normal_range_mp(r,k,part):
    """Independent erf and mpmath quadrature for small normal-range values."""
    r,k = map(mp.mpf,(r,k))
    def interval(u):
        return (mp.erf((u+r/2)/mp.sqrt(2))-mp.erf((u-r/2)/mp.sqrt(2)))/2
    d0 = interval(0)
    if part == 0:
        shift = mp.log(k*(k-1))-mp.log(2*mp.pi)-r*r/4+(k-2)*mp.log(d0)
        kernel = lambda u:mp.exp(-u*u)*(interval(u)/d0)**(k-2)
    else:
        shift = mp.log(k)-mp.log(2*mp.pi)/2-r*r/8+(k-1)*mp.log(d0)
        kernel = lambda u:mp.exp(-u*u/2+u*r/2)*(interval(u)/d0)**(k-1)
    return shift+mp.log(mp.quad(kernel,[-12,-2,-1,-.5,0,.5,1,2,12]))


@lru_cache(maxsize=None)
def range_moment(k,df):
    kernel = lambda r:df*math.log(r)+normal_range(r,k,0)
    peak = optimize.minimize_scalar(lambda r:-kernel(r),bounds=(1e-5,40.),method='bounded').x
    shift = kernel(peak)
    integral = integrate.quad(lambda r:math.exp(kernel(r)-shift),0,40,
                              points=[peak],epsabs=3e-12,epsrel=3e-12,limit=300)[0]
    return shift+math.log(integral)


def generate():
    cases,rows = [],[]
    def add(k,df,q,part,ref,end=None,loc=0.,scale=1.,label='adaptive'):
        cases.append([label,float(ref),3e-10])
        rows.extend([f'109 {k:.17g} {loc:.17g} {scale:.17g} '+
                     ' '.join(str(int(part == j)) for j in range(4)),f'{df:.17g}',
                     f'{q:.17g}'+(' '+f'{end:.17g}' if end is not None else '')])
    with mp.workdps(80):
        for k in (math.nextafter(10.,math.inf),11.,20.,37.5,50.,75.,100.):
            for df in (1.,10.,100.,math.inf):
                c = lower_constant(k,df)
                q = 1e-100
                add(k,df,q,1,c+(k-1)*mp.log(q),label='MP lower CDF')
                add(k,df,q,0,c+mp.log(k-1)+(k-2)*mp.log(q),label='MP lower PDF')
        for k in (20.,37.5,100.):
            for q in (.1,1.):
                for part in (0,1):
                    ref = normal_range_mp(q,k,part)
                    assert abs(float(ref)-normal_range(q,k,part)) < 3e-11
                    add(k,math.inf,q,part,ref,label='MP normal range')
            for df in (1.,100.,math.inf):
                for q in (.1,1.,4.,8.):
                    refs = distribution_logs(q,k,df)
                    for part in range(3):
                        add(k,df,q,part,refs[part])
                q = 4.;end = math.nextafter(q,math.inf)
                add(k,df,q,3,finite_range(q,k,df,0)+math.log(end-q),end=end,
                    label='adaptive adjacent')
                lc1 = finite_range(3.,k,df,1);lc2 = finite_range(4.,k,df,1)
                add(k,df,3.,3,logdiff(lc2,lc1),end=4.,label='adaptive interval')
                if not math.isinf(df):
                    lm = range_moment(k,df)
                    c = math.log(2)+df/2*math.log(df/2)-special.gammaln(df/2)
                    for q in (1e100,1e300):
                        add(k,df,q,2,c+lm-math.log(df)-df*math.log(q),label='moment remote SF')
                        add(k,df,q,0,c+lm-(df+1)*math.log(q),label='moment remote PDF')
            for q in (12.,40.):
                for part in (0,2):
                    add(k,math.inf,q,part,normal_range(q,k,part),label='adaptive infinite-df remote')
            print('Large-k fixtures ready',k,flush=True)
        k,df=100.,100.
        for scale in (1e-200,1e200):
            loc=.3*scale;x=loc+4*scale
            z=float((mp.mpf(x)-loc)/scale)
            refs=distribution_logs(z,k,df)
            for part in range(3):
                ref=refs[part]-(math.log(scale) if part == 0 else 0)
                add(k,df,x,part,ref,loc=loc,scale=scale,label='adaptive affine')
        x=1e300;scale=1e-200
        c=math.log(2)+df/2*math.log(df/2)-special.gammaln(df/2)+range_moment(k,df)
        lz=mp.log(mp.mpf(x)/scale)
        ld=c-(df+1)*lz-mp.log(scale)
        add(k,df,x,0,ld,scale=scale,label='moment overflowed PDF')
        add(k,df,x,2,c-mp.log(df)-df*lz,scale=scale,label='moment overflowed SF')
        end=math.nextafter(x,math.inf)
        add(k,df,x,3,ld+mp.log(mp.mpf(end)-x),end=end,scale=scale,label='moment overflowed adjacent')
    return cases,rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse',action='store_true')
    parser.add_argument('--profile',choices=('both','checked','optimized'),default='both')
    args=parser.parse_args()
    if args.reuse:
        cases=json.loads((CACHE/'large-k-cases.json').read_text())
        rows=(CACHE/'large-k-input.txt').read_text().splitlines()
    else:
        cases,rows=generate()
        (CACHE/'large-k-cases.json').write_text(json.dumps(cases))
        (CACHE/'large-k-input.txt').write_text('\n'.join(rows)+'\n')
    profiles=('checked','optimized') if args.profile == 'both' else (args.profile,)
    for profile in profiles:
        result=subprocess.run([str(CACHE/profile/'driver.exe')],input='\n'.join(rows)+'\n0\n',
                              text=True,capture_output=True,check=True)
        (CACHE/f'large-k-{profile}-actual.txt').write_text(result.stdout)
        values=list(map(float,result.stdout.split()));assert len(values)==len(cases)
        worst=0.
        for i,(actual,(label,expected,tol)) in enumerate(zip(values,cases)):
            error=0. if actual == expected else abs(actual-expected)/max(1.,abs(expected))
            if 0 < abs(expected) < 1e-10:
                error=abs(actual/expected-1)
            assert math.isfinite(error) and error<=tol,(profile,i,label,actual,expected,error)
            worst=max(worst,error)
        print(profile,len(cases),'large-k references passed; maximum error',worst,flush=True)
    def objective(scale):
        rv=stats.studentized_range(100.,100.,scale=scale)
        return -(sum(rv.logpdf([3.,4.,5.]))+rv.logcdf(2.5)+rv.logsf(7.)+
                 math.log(rv.cdf(4.75)-rv.cdf(4.5)))
    fit=optimize.minimize_scalar(objective,bounds=(.5,2.),method='bounded',options={'xatol':1e-10})
    print('Independent SciPy k=100 fit: scale',fit.x,'NNLF',fit.fun,flush=True)


if __name__ == '__main__':
    main()
