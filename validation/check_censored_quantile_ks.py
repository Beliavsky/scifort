#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent mathematical checks for Tukey lambda, KSone and fixed-branch von Mises."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess
import mpmath as mp
from scipy import stats, optimize

ROOT = Path(__file__).resolve().parents[1]
NAMES = ('tukeylambda', 'ksone', 'vonmises')


def tl(x, shape, loc, scale):
    x, lam, loc, scale = map(mp.mpf, (x, shape, loc, scale))
    z = (x-loc)/scale
    if lam > 0:
        lower = mp.mpf(float(loc)-float(scale)/float(lam))
        upper = mp.mpf(float(loc)+float(scale)/float(lam))
        if x < lower or x > upper:
            return -mp.inf, (-mp.inf if x < lower else mp.mpf(0)), (mp.mpf(0) if x < lower else -mp.inf)
        if x == lower or x == upper:
            ld = -mp.log(scale) if lam > 1 else -mp.log(2*scale) if lam == 1 else -mp.inf
            return ld, (-mp.inf if x == lower else mp.mpf(0)), (mp.mpf(0) if x == lower else -mp.inf)
        gap = (x-lower)/scale if x < loc else (upper-x)/scale
        if abs(z) > 1/(2*lam):
            # Finite affine endpoints follow the documented binary64 convention.
            z = mp.sign(z)*(1/lam-gap)
    if z == 0:
        return (lam-2)*mp.log(2)-mp.log(scale), -mp.log(2), -mp.log(2)
    if lam == 0:
        small = -mp.log1p(mp.exp(abs(z)))
    else:
        def quantile(s):
            logq = -s
            logp = mp.log1p(-mp.exp(-s))
            return (mp.expm1(lam*logp)-mp.expm1(lam*logq))/lam
        lo, hi = mp.log(2), mp.mpf(1)
        while quantile(hi) < abs(z):
            hi *= 2
        for _ in range(180):
            mid = (lo+hi)/2
            if quantile(mid) < abs(z):
                lo = mid
            else:
                hi = mid
        small = -(lo+hi)/2
    large = mp.log1p(-mp.exp(small))
    ld = -mp.log(mp.exp((lam-1)*small)+mp.exp((lam-1)*large))-mp.log(scale)
    return ld, (small if z < 0 else large), (large if z < 0 else small)


def ko(x, shape, loc, scale, side=0):
    n = int(shape)
    x, loc, scale = map(mp.mpf, (x, loc, scale))
    upper = mp.mpf(float(loc)+float(scale))
    if x < loc or x > upper:
        return -mp.inf, (-mp.inf if x < loc else mp.mpf(0)), (mp.mpf(0) if x < loc else -mp.inf)
    if x == loc or x == upper:
        ld = -mp.log(scale) if n == 1 else -mp.inf
        return ld, (-mp.inf if x == loc else mp.mpf(0)), (mp.mpf(0) if x == loc else -mp.inf)
    z, u = (x-loc)/scale, (upper-x)/scale
    if n == 1:
        small = mp.log(min(z, u))
        large = mp.log1p(-mp.exp(small))
        return -mp.log(scale), (small if z < u else large), (large if z < u else small)
    if z <= mp.mpf(1)/n:
        lp = mp.log(z)+(n-1)*mp.log1p(z)
        lq = mp.log1p(-mp.exp(lp))
    elif u <= mp.mpf(1)/n:
        lq = n*mp.log(u)
        lp = mp.log1p(-mp.exp(lq))
    else:
        sf = (1-z)**n
        for j in range(1, int(mp.floor(n*(1-z)))+1):
            sf += mp.binomial(n,j)*z*(z+mp.mpf(j)/n)**(j-1)*(1-z-mp.mpf(j)/n)**(n-j)
        lq = mp.log(sf)
        lp = mp.log1p(-sf)
    # One-sided derivatives at the first knot, with the documented midpoint value.
    knot = mp.mpf(float(loc)+float(scale)/n)
    if abs(z-mp.mpf(1)/n) < mp.mpf('1e-13'):
        z = mp.mpf(1)/n+(x-knot)/scale
    if z < mp.mpf(1)/n:
        pdf = (1+z)**(n-2)*(1+n*z)
    elif u <= mp.mpf(1)/n and x != knot:
        pdf = n*u**(n-1)
    else:
        z = max(z, mp.mpf(1)/n)
        pdf = n*(1-z)**(n-1)
        for j in range(1, min(n-2, int(mp.floor(n*(1-z))))+1):
            a, b = z+mp.mpf(j)/n, 1-z-mp.mpf(j)/n
            if b <= 0:
                continue
            term = mp.binomial(n,j)*z*a**(j-1)*b**(n-j)
            pdf -= term*(1/z+(j-1)/a-(n-j)/b)
        if x == knot and side == 0:
            pdf += mp.mpf('.5')
    return mp.log(pdf)-mp.log(scale), lp, lq


def vm(x, shape, loc, scale):
    x, k, loc, scale = map(mp.mpf, (x, shape, loc, scale))
    lower = mp.mpf(float(loc)-math.pi*float(scale))
    upper = mp.mpf(float(loc)+math.pi*float(scale))
    if x < lower or x > upper:
        return -mp.inf, (-mp.inf if x < lower else mp.mpf(0)), (mp.mpf(0) if x < lower else -mp.inf)
    z = (x-loc)/scale
    norm = mp.log(2*mp.pi*mp.besseli(0,k))
    ld = k*mp.cos(z)-norm-mp.log(scale)
    if x == lower or x == upper:
        return ld, (-mp.inf if x == lower else mp.mpf(0)), (mp.mpf(0) if x == lower else -mp.inf)
    # Integrate the nearer tail from the exact standardized angular endpoint.
    width = mp.pi-abs(z)
    base = k*mp.cos(z)-norm
    integral = mp.quad(lambda t: mp.exp(k*(mp.cos(abs(z)+width*t)-mp.cos(z))),
                       [0, .5, 1], maxdegree=6)
    small = base+mp.log(width)+mp.log(integral)
    large = mp.log1p(-mp.exp(small))
    return ld, (small if z < 0 else large), (large if z < 0 else small)


ORACLES = (tl, ko, vm)


def interval(a, b, family, shape, loc, scale):
    a, b = mp.mpf(a), mp.mpf(b)
    if family == 0 and shape > 0:
        a = max(a, mp.mpf(float(loc)-float(scale)/shape))
        b = min(b, mp.mpf(float(loc)+float(scale)/shape))
    elif family == 1:
        a, b = max(a, mp.mpf(loc)), min(b, mp.mpf(float(loc)+float(scale)))
    elif family == 2:
        a = max(a, mp.mpf(float(loc)-math.pi*float(scale)))
        b = min(b, mp.mpf(float(loc)+math.pi*float(scale)))
    if b <= a:
        return -mp.inf
    oracle = ORACLES[family]
    if family == 0 and shape > 0:
        if a == mp.mpf(float(loc)-float(scale)/shape):
            return tl(b,shape,loc,scale)[1]
        if b == mp.mpf(float(loc)+float(scale)/shape):
            return tl(a,shape,loc,scale)[2]
    if family == 1:
        if a == mp.mpf(loc):
            return ko(b,shape,loc,scale)[1]
        if b == mp.mpf(float(loc)+float(scale)):
            return ko(a,shape,loc,scale)[2]
    if family == 0 and shape == 0:
        if a >= loc:
            large, small = tl(a,shape,loc,scale)[2], tl(b,shape,loc,scale)[2]
        else:
            large, small = tl(b,shape,loc,scale)[1], tl(a,shape,loc,scale)[1]
        return large+mp.log1p(-mp.exp(small-large))
    # Integrate positive relative density, splitting the KS first-knot jump.
    pieces = [a, b]
    if family == 1:
        knot = mp.mpf(float(loc)+float(scale)/shape)
        if a < knot < b:
            pieces.insert(1, knot)
    total = mp.mpf(0)
    for lo, hi in zip(pieces[:-1], pieces[1:]):
        width = hi-lo
        mid = (lo+hi)/2
        base = oracle(mid, shape, loc, scale)[0]
        if not mp.isfinite(base):
            return -mp.inf
        relative_width = width/mp.mpf(scale)
        if family == 0 and shape < 0:
            relative_width = width/(mp.mpf(scale)+abs(mid-mp.mpf(loc)))
        if relative_width < mp.mpf('1e-10'):
            # The density is smooth within each split piece. Midpoint error is O(width**2).
            mass = mp.exp(base)*width
        elif relative_width < mp.mpf('.001'):
            weights = (1,4,2,4,1)
            points = (mp.mpf('1e-20'),mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1-mp.mpf('1e-20'))
            total_weight = sum(w*mp.exp(oracle(lo+width*t,shape,loc,scale)[0]-base)
                               for t,w in zip(points,weights))/12
            mass = mp.exp(base)*width*total_weight
        else:
            mass = mp.exp(base)*width*mp.quad(
                lambda t: mp.exp(oracle(lo+width*t, shape, loc, scale)[0]-base),
                [0, .5, 1], maxdegree=4)
        total += mass
    return mp.log(total)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse', action='store_true')
    args = parser.parse_args()
    cache = ROOT/'build'/'quantile_ks-reference'
    cache.mkdir(exist_ok=True)
    cases, inputs = [], []
    def add(family, shape, x, part, loc=0., scale=1., end=None, label='reference'):
        expected = (ORACLES[family](x,shape,loc,scale)[part] if part < 3
                    else interval(x,end,family,shape,loc,scale))
        cases.append((NAMES[family],label,float(expected)))
        counts = [0]*4
        counts[part] = 1
        inputs.append(f'{102+family} {shape:.17g} {loc:.17g} {scale:.17g} '+
                      ' '.join(map(str,counts)))
        inputs.append(f'{x:.17g}'+(' '+f'{end:.17g}' if end is not None else ''))
    if args.reuse:
        cases = json.loads((cache/'cases.json').read_text())
        inputs = (cache/'input.txt').read_text().splitlines()
    else:
        with mp.workdps(100):
            for shape in (-5., -1., -1e-6, 0., 1e-6, .14, 1., 2., 5.):
                for x in (-.1, .01, .1):
                    for part in range(3):
                        add(0,shape,x,part,label='Tukey shape grid')
                add(0,shape,.1,3,end=math.nextafter(.1, math.inf),label='Tukey adjacent')
            for shape in (-5., -1., -1e-6, 0.):
                for x in (-1e300, -1e100, 1e100, 1e300):
                    add(0,shape,x,0,label='Tukey extreme density')
                    add(0,shape,x,1 if x < 0 else 2,label='Tukey extreme tail')
                    add(0,shape,x,3,end=math.nextafter(x,math.inf),label='Tukey extreme adjacent')
            for shape in (.14, 1., 2., 5.):
                endpoint = 1/shape
                x = math.nextafter(endpoint, 0.)
                add(0,shape,x,2,label='Tukey endpoint tail')
                add(0,shape,x,3,end=endpoint,label='Tukey endpoint interval')
            for n in (1., 2., 3., 10., 100.):
                for x in (1e-300, .01, .1, .3, .5, math.nextafter(1.,0.)):
                    for part in range(3):
                        add(1,n,x,part,label='KS sample/tail grid')
                knot = 1/n
                if n > 1:
                    for x,end in ((math.nextafter(knot,0.),knot),
                                  (knot,math.nextafter(knot,math.inf))):
                        add(1,n,x,3,end=end,label='KS first-knot interval')
                add(1,n,math.nextafter(1.,0.),3,end=1.,label='KS endpoint interval')
            for k in (0., 1., 100.):
                for x in (-3., -.1, .1, 3.):
                    for part in range(3):
                        add(2,k,x,part,label='angular tails')
                    add(2,k,x,3,end=math.nextafter(x,math.inf),label='angular adjacent')
            for family, shape in enumerate((-.3,10.,1.)):
                for scale in (1e-200,1.,1e200):
                    loc = -.3*scale
                    x = loc+.4*scale
                    for part in range(3):
                        add(family,shape,x,part,loc,scale,label='affine')
                    add(family,shape,x,3,loc,scale,math.nextafter(x,math.inf),'affine adjacent')
            for scale in (1e-200,1.,1e200):
                for sign in (-1,1):
                    loc = sign*1e12*scale
                    x = loc+.7*scale
                    add(0,-.3,x,3,loc,scale,math.nextafter(x,math.inf),'large-location interval')
            x = math.nextafter(0.,1.)
            for family, shape in enumerate((-.3,100.,1.)):
                add(family,shape,x,3,end=2*x,label='subnormal interval')
            for scale in (1e-200,1.,1e200):
                loc = -.3*scale
                endpoint = loc+scale/.14
                x = math.nextafter(endpoint,-math.inf)
                add(0,.14,x,2,loc,scale,label='positive affine endpoint')
                add(0,.14,x,3,loc,scale,endpoint,'positive affine interval')
            add(0,-5.,1e300,1,label='tiny complement')
            add(0,-1.,1e300,1,label='tiny complement')
            add(0,.14,math.nextafter(1/.14,0.),1,label='tiny complement')
            add(1,100.,math.nextafter(0.,1.),2,label='tiny complement')
            add(0,0.,400.,1,label='tiny complement')
            add(1,100.,1e-300,2,label='tiny complement')
            add(2,100.,-3.,2,label='tiny complement')
            print('Multiprecision fixtures ready:',len(cases),flush=True)
        (cache/'cases.json').write_text(json.dumps(cases),encoding='ascii')
        (cache/'input.txt').write_text('\n'.join(inputs)+'\n',encoding='ascii')
    for name,flags in (('checked',['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']),
                       ('optimized',['-O3','-fcheck=bounds'])):
        out = cache/name
        out.mkdir(exist_ok=True)
        cmd = ['gfortran',*flags,'-J'+str(out),'-I'+str(out),'-Ibuild/gfortran_2806A042DAA8BBF7',
               'src/math.f90','src/special/log_gamma.f90','src/empirical_distribution.f90','src/censored.f90',
               'validation/censored_reference_driver.f90',
               'build/gfortran_3C159EB802C22FB9/scifort/libscifort.a','-o',str(out/'driver.exe')]
        print(subprocess.list2cmdline(cmd),flush=True)
        subprocess.run(cmd,cwd=ROOT,check=True)
        r = subprocess.run([str(out/'driver.exe')],input='\n'.join(inputs)+'\n',
                           text=True,capture_output=True)
        if r.returncode:
            n = len(r.stdout.split())
            raise RuntimeError((n,cases[n:n+1],r.stderr))
        actual = list(map(float,r.stdout.split()))
        assert len(actual) == len(cases)
        worst, tiny, failures = 0., 0., []
        for v,case in zip(actual,cases):
            family,label,expected = case
            if not math.isfinite(expected):
                e = 0. if v == expected else math.inf
            elif label == 'tiny complement':
                e = abs(v/expected-1)
            else:
                e = abs(v-expected)/max(1.,abs(expected))
            if math.isnan(v) or e > 3e-11:
                failures.append((case,v,e))
            if label == 'tiny complement':
                tiny = max(tiny,e)
            else:
                worst = max(worst,e)
        assert not failures,(name,failures)
        print(name,len(cases),'MP comparisons passed; scaled log',worst,'tiny relative',tiny,flush=True)
        scipy_check(out/'driver.exe')
    sample = stats.CensoredData(uncensored=[.2,.4,.6],left=[.1],right=[.8],interval=[[.3,.35],[.5,.55]])
    for f,sh in ((stats.tukeylambda,.14),(stats.vonmises_line,1.)):
        opt = optimize.minimize_scalar(lambda scale:f._penalized_nnlf((sh,0.,scale),sample),
                                      bounds=(.3,3.),method='bounded',options={'xatol':1e-13})
        print(f.name,'SciPy scale',opt.x,'NNLF',opt.fun,flush=True)


def scipy_check(driver):
    rng = random.Random(2026111807)
    inputs, refs = [], []
    for family in range(3):
        for _ in range(25):
            shape = rng.uniform(-3,3) if family == 0 else float(rng.randint(1,100)) if family == 1 else rng.uniform(0,10)
            x = rng.uniform(-.2,.2) if family != 1 else rng.uniform(.02,.98)
            distribution = (stats.tukeylambda,stats.ksone,stats.vonmises_line)[family]
            for part,method in enumerate(('logpdf','logcdf','logsf')):
                ref = getattr(distribution,method)(x,shape)
                # Only stable ordinary probabilities are compared directly to SciPy.
                if not math.isfinite(ref) or (part != 0 and abs(ref) > 15):
                    continue
                counts = [0]*4
                counts[part] = 1
                inputs.extend([f'{102+family} {shape} 0 1 '+ ' '.join(map(str,counts)),str(x)])
                refs.append(float(ref))
    r = subprocess.run([str(driver)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    actual = list(map(float,r.stdout.split()))
    assert len(actual) == len(refs)
    for v,e in zip(actual,refs):
        assert abs(v-e) < 3e-9*max(1.,abs(e)),(v,e)
    print('Seeded ordinary SciPy comparisons passed:',len(refs),flush=True)


if __name__ == '__main__':
    main()
