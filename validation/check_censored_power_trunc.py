#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent generalized-gamma and truncated-law censoring checks."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess
import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('gengamma', 'truncnorm', 'truncweibull_min')


def normal_mass(a, b):
    if a >= 0:
        return (mp.erfc(a/mp.sqrt(2))-mp.erfc(b/mp.sqrt(2)))/2
    return (mp.erfc(-b/mp.sqrt(2))-mp.erfc(-a/mp.sqrt(2)))/2


def density(x, family, shapes, loc, scale):
    z = (mp.mpf(x)-loc)/scale
    a, b, c = map(mp.mpf, shapes)
    if family == 0:
        return mp.log(abs(b))+(a*b-1)*mp.log(z)-z**b-mp.loggamma(a)-mp.log(scale)
    if family == 1:
        return -z*z/2-mp.log(2*mp.pi)/2-mp.log(normal_mass(a,b))-mp.log(scale)
    span = c**a-b**a
    shift = z**a-b**a
    return mp.log(a)+(a-1)*mp.log(z)-shift-mp.log(-mp.expm1(-span))-mp.log(scale)


def logs(x, family, shapes, loc, scale):
    z = (mp.mpf(x)-loc)/scale
    a, b, c = map(mp.mpf, shapes)
    if family == 0:
        y = z**b
        if a == 1:
            q = mp.exp(-y)
            p = -mp.expm1(-y)
        else:
            p = mp.gammainc(a,0,y)/mp.gamma(a)
            q = mp.gammainc(a,y,mp.inf)/mp.gamma(a)
        if b < 0:
            p, q = q, p
    elif family == 1:
        denominator = normal_mass(a,b)
        p, q = normal_mass(a,z)/denominator, normal_mass(z,b)/denominator
    else:
        denominator = -mp.expm1(-(c**a-b**a))
        p = -mp.expm1(-(z**a-b**a))/denominator
        q = mp.exp(-(z**a-b**a))*(-mp.expm1(-(c**a-z**a)))/denominator
    lp = mp.log(p) if p < mp.mpf('.5') else mp.log1p(-q)
    lq = mp.log(q) if q < mp.mpf('.5') else mp.log1p(-p)
    return density(x,family,shapes,loc,scale), lp, lq


def interval(lo, hi, family, shapes, loc, scale):
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    if family == 0 and shapes[0] == 1:
        power = mp.mpf(shapes[1])
        first, last = ((lo-loc)/scale)**power, ((hi-loc)/scale)**power
        if first > last:
            first, last = last, first
        return -first+mp.log(-mp.expm1(-(last-first)))
    if (hi-lo)/scale < mp.mpf('1e-5'):
        ref = density(lo,family,shapes,loc,scale)
        return ref+mp.log(hi-lo)+mp.log(mp.quad(
            lambda t: mp.exp(density(lo+t*(hi-lo),family,shapes,loc,scale)-ref),
            [0,1],method='gauss-legendre',maxdegree=3))
    _, lp, lq = logs(lo,family,shapes,loc,scale)
    _, up, uq = logs(hi,family,shapes,loc,scale)
    return up+mp.log(-mp.expm1(lp-up)) if up < lq else lq+mp.log(-mp.expm1(uq-lq))


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
        inputs.append(f'{78+family} {shapes[0]:.17g} {float(loc):.17g} {float(scale):.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        inputs.append(f'{shapes[1]:.17g}')
        if family == 2:
            inputs.append(f'{shapes[2]:.17g}')
        for group in (exact,left,right,*intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    with mp.workdps(100):
        for family in range(3):
            for _ in range(args.cases):
                if family == 0:
                    shapes = (10**rng.uniform(-1,1),rng.choice((-1,1))*rng.uniform(.3,3),0.)
                elif family == 1:
                    shapes = (-rng.uniform(.5,2),rng.uniform(3,5),0.)
                else:
                    shapes = (rng.uniform(.3,3),rng.uniform(.01,.1),rng.uniform(3,5))
                scale = 10**rng.uniform(-200,200)
                loc = scale*rng.uniform(-.3,.3)
                point = lambda x: loc+scale*x
                lo = point(rng.uniform(.3,2))
                hi = lo+scale*10**rng.uniform(-15,-2)
                if hi == lo:
                    hi = math.nextafter(lo,math.inf)
                add(family,shapes,loc,scale,exact=[point(.7)],left=[point(.2)],
                    right=[point(2.5)],intervals=[(lo,hi),(point(.4),point(1.7))])
            print(f'{FAMILIES[family]}: mixed references ready',flush=True)
    with mp.workdps(1100):
        for power in (-2.,-.5,.5,2.):
            for shape in (1e-300,.01,1.,100.):
                add(0,(shape,power,0.),exact=[.7],left=[.7],right=[.7],label='gamma shape/power')
                add(0,(shape,power,0.),intervals=[(.7,math.nextafter(.7,math.inf))],
                    label='gamma adjacent')
        for power in (-2.,2.):
            for x in (1e-100,1e100):
                add(0,(1.,power,0.),exact=[x],left=[x],right=[x],label='gamma extreme tail')
                add(0,(1.,power,0.),intervals=[(x,math.nextafter(x,math.inf))],
                    label='gamma extreme adjacent')
        for family,shapes in ((0,(1.7,-1.3,0.)),(1,(-1.,3.,0.)),(2,(1.3,.1,3.))):
            for scale in (1e-300,1e300):
                x = .7*scale
                add(family,shapes,scale=scale,exact=[x],left=[x],right=[x],label='scale')
                add(family,shapes,scale=scale,intervals=[(x,math.nextafter(x,math.inf))],
                    label='adjacent scale')
        for lo in (-100.,-81.,40.,80.):
            hi = lo+1.
            x = lo+.3
            add(1,(lo,hi,0.),exact=[x],left=[x],right=[x],label='normal far truncation')
        hi = math.nextafter(math.nextafter(1.,math.inf),math.inf)
        for family,shapes in ((1,(1.,hi,0.)),(2,(2.,1.,hi))):
            add(family,shapes,intervals=[(1.,math.nextafter(1.,math.inf))],label='two-ULP window')
        for power,lo,hi in ((.01,.1,3.),(2.,100.,101.),(2.,1e100,1.01e100)):
            x = lo+(hi-lo)*.3
            add(2,(power,lo,hi),exact=[x],left=[x],right=[x],label='Weibull hazard conditioning')
        for bound in (1e-100,1e-300):
            x = math.nextafter(bound,0.)
            add(2,(1e6,0.,bound),exact=[x],left=[x],right=[x],label='tiny Weibull normalization')
    cache = ROOT/'build'/'power_trunc_reference'
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
    for family,name in enumerate(FAMILIES):
        parameters = ((lambda v:(v,-1.3,0.,1.)),(lambda v:(-1.,3.,v,1.)),
                      (lambda v:(v,.1,3.,0.,1.)))[family]
        bounds = (-.1,.1) if family == 1 else (.3,5.)
        fit = optimize.minimize_scalar(lambda v:getattr(stats,name)._penalized_nnlf(parameters(v),sample),
                                       bounds=bounds,method='bounded',options={'xatol':1e-13})
        assert fit.success
        print(f'{name}: SciPy parameter={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
