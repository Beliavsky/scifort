#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent alpha, GEV and Crystal Ball censored likelihood checks."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess
import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('alpha', 'genextreme', 'crystalball')


def phi(x):
    return mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)


def cdf(x):
    return mp.erfc(-x/mp.sqrt(2))/2


def sf(x):
    return mp.erfc(x/mp.sqrt(2))/2


def mass(a, b):
    return sf(a)-sf(b) if a >= 0 else cdf(b)-cdf(a)


def cb_constants(beta, m):
    tail = m/beta/(m-1)*mp.exp(-beta*beta/2)
    return tail, tail+mp.sqrt(2*mp.pi)*cdf(beta)


def density(x, family, shapes, loc, scale):
    z = (mp.mpf(x)-loc)/scale
    a, b = map(mp.mpf, shapes)
    if family == 0:
        return -(a-1/z)**2/2-mp.log(2*mp.pi)/2-mp.log(cdf(a))-2*mp.log(z)-mp.log(scale)
    if family == 1:
        lt = mp.log1p(-a*z)/a if a else -z
        lu = mp.log1p(-a*z) if a else mp.mpf(0)
        return -mp.exp(lt)+lt-lu-mp.log(scale)
    _, den = cb_constants(a,b)
    val = -z*z/2 if z > -a else -a*a/2-b*mp.log1p(a*(-a-z)/b)
    return val-mp.log(den)-mp.log(scale)


def logs(x, family, shapes, loc, scale):
    z = (mp.mpf(x)-loc)/scale
    a, b = map(mp.mpf, shapes)
    if family == 0:
        p, q = cdf(a-1/z)/cdf(a), mass(a-1/z,a)/cdf(a)
    elif family == 1:
        t = mp.exp(mp.log1p(-a*z)/a if a else -z)
        p, q = mp.exp(-t), -mp.expm1(-t)
    else:
        tail, den = cb_constants(a,b)
        if z <= -a:
            power = (1-b)*mp.log1p(a*(-a-z)/b)
            p = tail*mp.exp(power)/den
            q = (mp.sqrt(2*mp.pi)*cdf(a)+tail*(-mp.expm1(power)))/den
        else:
            p = (tail+mp.sqrt(2*mp.pi)*mass(-a,z))/den
            q = mp.sqrt(2*mp.pi)*sf(z)/den
    lp = mp.log(p) if p < mp.mpf('.5') else mp.log1p(-q)
    lq = mp.log(q) if q < mp.mpf('.5') else mp.log1p(-p)
    return density(x,family,shapes,loc,scale), lp, lq


def interval(lo, hi, family, shapes, loc, scale):
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    za, zb = (lo-loc)/scale, (hi-loc)/scale
    a, b = map(mp.mpf,shapes)
    if family == 0:
        return mp.log(mass(a-1/za,a-1/zb)/cdf(a))
    if family == 1:
        ta = mp.exp(mp.log1p(-a*za)/a if a else -za)
        tb = mp.exp(mp.log1p(-a*zb)/a if a else -zb)
        return -tb+mp.log(-mp.expm1(-(ta-tb)))
    if zb <= -a:
        tail, den = cb_constants(a,b)
        da, db = b/a-a-za, b/a-a-zb
        return mp.log(tail/den)+(1-b)*mp.log(db/(b/a))+mp.log(-mp.expm1((1-b)*mp.log(da/db)))
    if za >= -a:
        _, den = cb_constants(a,b)
        return mp.log(mass(za,zb)*mp.sqrt(2*mp.pi)/den)
    tail, den = cb_constants(a,b)
    power = (1-b)*mp.log1p(a*(-a-za)/b)
    return mp.log((tail*(-mp.expm1(power))+mp.sqrt(2*mp.pi)*mass(-a,zb))/den)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--module-dir', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--cases', type=int, default=100)
    args = parser.parse_args()
    rng = random.Random(20261117)
    cases, inputs = [], []

    def add(family, shapes, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed'):
        loc, scale = mp.mpf(loc), mp.mpf(scale)
        value = sum(logs(x,family,shapes,loc,scale)[0] for x in exact)
        value += sum(logs(x,family,shapes,loc,scale)[1] for x in left)
        value += sum(logs(x,family,shapes,loc,scale)[2] for x in right)
        value += sum(interval(a,b,family,shapes,loc,scale) for a,b in intervals)
        expected = float(value)
        assert math.isfinite(expected), (family,label,value)
        cases.append((FAMILIES[family],label,expected))
        inputs.append(f'{81+family} {shapes[0]:.17g} {float(loc):.17g} {float(scale):.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        if family == 2:
            inputs.append(f'{shapes[1]:.17g}')
        for group in (exact,left,right,*intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    with mp.workdps(100):
        for family in range(3):
            for _ in range(args.cases):
                shapes = ((10**rng.uniform(-1,1),0.) if family == 0 else
                          (rng.choice((-1,0,1))*rng.uniform(.01,.3),0.) if family == 1 else
                          (10**rng.uniform(-1,1),1+10**rng.uniform(-1,2)))
                scale = 10**rng.uniform(-200,200)
                loc = scale*rng.uniform(-.3,.3)
                point = lambda x: loc+scale*x
                lo = point(rng.uniform(.3,2))
                hi = lo+scale*10**rng.uniform(-15,-2)
                if hi == lo:
                    hi = math.nextafter(lo,math.inf)
                intervals = [(lo,hi),(point(.4),point(1.7))]
                if family == 2:
                    intervals += [(point(-shapes[0]-2),point(-shapes[0]+.3))]
                add(family,shapes,loc,scale,exact=[point(.7)],left=[point(.2)],
                    right=[point(2.5)],intervals=intervals)
            print(f'{FAMILIES[family]}: mixed references ready',flush=True)
        for _ in range(100):
            shapes = (rng.uniform(.1,10),rng.uniform(1.1,100))
            scale = 10**rng.uniform(-200,200)
            loc = scale*rng.uniform(-.3,.3)
            transition = loc-scale*shapes[0]
            lo, hi = math.nextafter(transition,-math.inf), math.nextafter(transition,math.inf)
            add(2,shapes,loc,scale,intervals=[(lo,hi)],label='adjacent affine transition')
    with mp.workdps(1100):
        for family,shapes in ((0,(1.,0.)),(1,(0.,0.)),(2,(2.,3.))):
            for scale in (1e-300,1e300):
                x = .7*scale
                add(family,shapes,scale=scale,exact=[x],left=[x],right=[x],label='scale')
                add(family,shapes,scale=scale,intervals=[(x,math.nextafter(x,math.inf))],
                    label='adjacent scale')
        for shape in (1e-6,1.,100.):
            for x in (1e-100,1e100,1e300):
                add(0,(shape,0.),exact=[x],left=[x],right=[x],label='alpha tails')
                add(0,(shape,0.),intervals=[(x,math.nextafter(x,math.inf))],label='alpha adjacent')
        for shape in (-100.,-1.,-1e-6,0.,1e-6,1.,100.):
            x = .001
            add(1,(shape,0.),exact=[x],left=[x],right=[x],label='GEV shape range')
            add(1,(shape,0.),intervals=[(x,math.nextafter(x,math.inf))],label='GEV adjacent')
        for x in (-100.,100.,1000.):
            add(1,(0.,0.),exact=[x],left=[x],right=[x],label='Gumbel far tails')
            add(1,(0.,0.),intervals=[(x,math.nextafter(x,math.inf))],label='Gumbel adjacent')
        for beta in (1e-6,2.,100.):
            for m in (math.nextafter(1.,math.inf),3.,1e6):
                for x in (-beta-1.,40.):
                    add(2,(beta,m),exact=[x],left=[x],right=[x],label='Crystal Ball range')
                    add(2,(beta,m),intervals=[(x,math.nextafter(x,math.inf))],label='Crystal Ball adjacent')
        for x in (-1e100,-1e300):
            add(2,(2.,3.),exact=[x],left=[x],right=[x],label='power tail')
            add(2,(2.,3.),intervals=[(x,math.nextafter(x,math.inf))],label='power adjacent')
        for family,shapes,x in ((0,(1.,0.),1e300),(1,(-1.,0.),1e300),
                                (1,(-1e-6,0.),1e300),(2,(2.,3.),-1e300)):
            add(family,shapes,scale=1e-300,exact=[x],left=[x],right=[x],
                label='overflowed standardized coordinate')
            add(family,shapes,scale=1e-300,intervals=[(x,math.nextafter(x,math.inf))],
                label='overflowed standardized adjacent')
    cache = ROOT/'build'/'shape_reference'
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
    sample = stats.CensoredData(uncensored=[.5,.8,1.1],left=[.4],right=[1.7],
                                interval=[[.6,.7],[.9,1.]])
    for name,shapes in zip(FAMILIES,((1.7,),(.2,),(1.5,2.7))):
        fit = optimize.minimize_scalar(
            lambda v:getattr(stats,name)._penalized_nnlf((*shapes,0.,v),sample),
            bounds=(.3,3.),method='bounded',options={'xatol':1e-13})
        assert fit.success
        print(f'{name}: SciPy scale={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
