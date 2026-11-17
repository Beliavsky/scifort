#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-precision folded normal, folded Cauchy, and generalized Pareto checks."""
import argparse
import importlib.util
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('reference', ROOT/'tools/generate_censored_reference.py')
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)
FAMILIES = ('foldnorm', 'foldcauchy', 'genpareto')


def softplus(q):
    return max(q, 0)+mp.log1p(mp.exp(-abs(q)))


def logadd(a, b):
    if a == -mp.inf: return b
    if b == -mp.inf: return a
    return max(a,b)+mp.log1p(mp.exp(-abs(a-b)))


def difference(a,b):
    if a == b: return -mp.inf
    return a+mp.log(-mp.expm1(b-a))


def normal_mass(a,b):
    pa,qa = reference.transformed_normal_logs(a)
    pb,qb = reference.transformed_normal_logs(b)
    if b <= 0: return difference(pb,pa)
    if a >= 0: return difference(qa,qb)
    return mp.log1p(-mp.exp(pa)-mp.exp(qb))


def cauchy_mass(a,b):
    if a >= 0 or b <= 0:
        return mp.log(mp.atan((b-a)/(1+a*b))/mp.pi)
    return mp.log((mp.atan(b)-mp.atan(a))/mp.pi)


def logs(x, family, shape, second, loc, scale):
    x,c,mu,s = map(mp.mpf,(x,shape,loc,scale))
    z = (x-mu)/s
    if z < 0: return -mp.inf,-mp.inf,mp.mpf(0)
    if family < 2:
        a,b = z-c,z+c
        if family == 0:
            density = logadd(-a*a/2,-b*b/2)-mp.log(2*mp.pi)/2-mp.log(s)
            lp = normal_mass(-b,a) if z else -mp.inf
            lq = logadd(reference.transformed_normal_logs(a)[1],
                        reference.transformed_normal_logs(b)[1])
        else:
            density = mp.log(1/(1+a*a)+1/(1+b*b))-mp.log(mp.pi*s)
            lp = cauchy_mass(-b,a) if z else -mp.inf
            lq = mp.log((mp.atan2(1,a)+mp.atan2(1,b))/mp.pi)
    else:
        if c < 0:
            width = mp.mpf(float(scale/(-shape)))
            w = (mu+width-x)/width
            if w < 0: return -mp.inf,mp.mpf(0),-mp.inf
            if not w:
                density = mp.inf if c < -1 else (-mp.inf if c > -1 else -mp.log(width))
                return density,mp.mpf(0),-mp.inf
            h = -mp.log(w)/(-c)
            density = -mp.log(-c*width)+(1/(-c)-1)*mp.log(w)
        else:
            h = mp.log1p(c*z)/c if c else z
            density = -h-mp.log1p(c*z)-mp.log(s)
        lp = mp.log(-mp.expm1(-h))
        lq = -h
    if lp <= lq: lq = mp.log1p(-mp.exp(lp))
    else: lp = mp.log1p(-mp.exp(lq))
    return density,lp,lq


def interval(a,b,family,shape,second,loc,scale):
    lo,hi,c,s,mu = map(mp.mpf,(a,b,shape,scale,loc))
    lo = max(lo,mu)
    if family < 2:
        za,zb = (lo-mu)/s,(hi-mu)/s
        mass = normal_mass if family == 0 else cauchy_mass
        return logadd(mass(za-c,zb-c),mass(za+c,zb+c))
    _,pa,qa=logs(float(lo),family,shape,second,loc,scale)
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
        inputs.append(f'{63+family} {shape:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    for family in range(3):
        for _ in range(args.cases):
            shape = rng.uniform(-2,3) if family == 2 else rng.uniform(0,10)
            scale = 10**rng.uniform(-100,100)
            width = scale/(-shape) if family == 2 and shape < 0 else scale
            loc = width*rng.uniform(-3,3)
            point = lambda p: loc+width*p
            exact = [point(rng.uniform(.01,.95)) for _ in range(3)]
            left = [point(rng.uniform(.01,.4))]
            right = [point(rng.uniform(.6,.99))]
            a = point(rng.uniform(.01,.9))
            b = a+width*10**rng.uniform(-14,-2)
            if b == a: b = math.nextafter(a,math.inf)
            add(family,shape,loc=loc,scale=scale,exact=exact,left=left,right=right,
                intervals=[(a,b),(point(.1),point(.8))])
    with mp.workdps(1100):
        tiny = math.nextafter(0.,math.inf)
        for family in range(3):
            shapes = (0.,tiny,1e-300,1.,100.) if family < 2 else (-1e300,-2.,-1.,-.5,-1e-300,0.,tiny,1e-300,1.,1e300)
            for shape in shapes:
                x = .5/(-shape) if family == 2 and shape < 0 else .5
                add(family,shape,exact=[x],left=[x],right=[x],label='shape stress')
                add(family,shape,intervals=[(x,math.nextafter(x,math.inf))],label='adjacent shape')
            for scale in (1e-300,1e300):
                x = .5*scale
                add(family,0.,scale=scale,exact=[x],left=[x],right=[x],label='scale stress')
                add(family,0.,scale=scale,intervals=[(x,math.nextafter(x,math.inf))],label='adjacent scale')
                add(family,0.,scale=scale,intervals=[(0.,tiny)],label='subnormal interval')
            for x in (tiny,1e-300,1e-15):
                add(family,1.,exact=[x],left=[x],right=[x],label='lower gaps')
                add(family,1.,intervals=[(x,math.nextafter(x,math.inf))],label='lower intervals')
        for family in (0,1):
            for shape in (1e100,1e150,1e300):
                add(family,shape,exact=[shape],left=[shape],right=[shape],label='large fold center')
            add(family,1.,intervals=[(tiny,1e100)],label='nearly full interval')
        add(0,0.,right=[1e150],label='extreme normal tail')
        add(1,1e300,exact=[tiny],left=[tiny],right=[tiny],label='tiny width huge offset')
        add(1,0.,scale=tiny,exact=[1.],right=[1.],intervals=[(1.,2.)],label='overflowed standardized distance')
        add(2,tiny,scale=tiny,exact=[tiny],left=[tiny],right=[tiny],label='tiny positive shape and scale')
        add(2,-tiny,scale=tiny,exact=[tiny],left=[tiny],right=[tiny],label='tiny negative shape and scale')
        add(2,-1.,scale=21*tiny,exact=[10*tiny],intervals=[(10*tiny,11*tiny)],label='subnormal bounded width')
    for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow']),
                        ('optimized', ['-O3', '-fcheck=bounds'])]:
        out = ROOT/'build'/'fold_pareto_reference'/mode
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
            raise RuntimeError(result.stderr)
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
    sample = stats.CensoredData(uncensored=[.2,.4,.6],left=[.1],right=[.7],
                               interval=[[.25,.35],[.5,.55]])
    for name in FAMILIES:
        distribution=getattr(stats,name)
        bounds=(-1.,3.) if name=='genpareto' else (0.,3.)
        fit=optimize.minimize_scalar(lambda c: distribution._penalized_nnlf((c,0.,1.3),sample),
            bounds=bounds,method='bounded',options={'xatol':1e-13})
        assert fit.success
        print(f'{name}: SciPy shape={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__=='__main__':
    main()
