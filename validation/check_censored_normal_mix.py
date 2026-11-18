#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-precision normal-mixture censoring checks."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('exponnorm', 'skewnorm', 'dpareto_lognorm')


def phi(z):
    return mp.exp(-z*z/2)/mp.sqrt(2*mp.pi)


def cdf(z):
    return mp.erfc(-z/mp.sqrt(2))/2


def emg_h(z, r):
    return mp.exp(r*r/2-r*z)*cdf(z-r)


def logdensity(x, family, shapes, loc, scale):
    x, loc, scale = map(mp.mpf, (x, loc, scale))
    z = (x-loc)/scale
    if family == 0:
        r = 1/mp.mpf(shapes[0])
        return mp.log(r*emg_h(z, r)/scale)
    if family == 1:
        return mp.log(2*phi(z)*cdf(mp.mpf(shapes[0])*z)/scale)
    u, sigma, a, b = map(mp.mpf, shapes)
    if z <= 0:
        return -mp.inf
    t = (mp.log(z)-u)/sigma
    return mp.log(a*b/(a+b)*(emg_h(t, a*sigma)+emg_h(-t, b*sigma))/(x-loc))


def skew_lower(z, a):
    if a == 0:
        return cdf(z)
    if a < 0:
        return 2*cdf(z)-skew_lower(z, -a)
    base = mp.log(2*phi(z)*cdf(a*z))
    rate = mp.sqrt(1+a*a)*(1+mp.sqrt(1+a*a)*abs(z))
    def relative(v):
        t = z-v/rate
        return mp.exp(mp.log(2*phi(t)*cdf(a*t))-base)
    return mp.exp(base)/rate*mp.quad(relative, [0, 1, mp.inf])


def oracle(x, family, shapes, loc, scale):
    x, loc, scale = map(mp.mpf, (x, loc, scale))
    s = tuple(map(mp.mpf, shapes))
    z = (x-loc)/scale
    if family == 0:
        r = 1/s[0]
        h = emg_h(z, r)
        density = r*h/scale
        p, q = cdf(z)-h, cdf(-z)+h
    elif family == 1:
        a = s[0]
        density = 2*phi(z)*cdf(a*z)/scale
        if z <= 0:
            p = skew_lower(z, a)
            q = 1-p
        else:
            q = skew_lower(-z, -a)
            p = 1-q
    else:
        u, sigma, a, b = s
        if z <= 0:
            return -mp.inf, -mp.inf, mp.mpf(0)
        t = (mp.log(z)-u)/sigma
        plus, minus = emg_h(t, a*sigma), emg_h(-t, b*sigma)
        # Signed CDF identities integrated from the density, unlike the positive mixture.
        p = cdf(t)+a/(a+b)*minus-b/(a+b)*plus
        q = cdf(-t)+b/(a+b)*plus-a/(a+b)*minus
        density = a*b/(a+b)*(plus+minus)/(x-loc)
    lp = mp.log(p) if p < mp.mpf('.5') else mp.log1p(-q)
    lq = mp.log(q) if q < mp.mpf('.5') else mp.log1p(-p)
    return mp.log(density), lp, lq


def interval(a, b, family, shapes, loc, scale):
    a, b = mp.mpf(a), mp.mpf(b)
    left_pdf = logdensity(a, family, shapes, loc, scale)
    right_pdf = logdensity(b, family, shapes, loc, scale)
    if abs(left_pdf-right_pdf) < 1:
        # Finite quadrature of relative density retains masses far below 1e-300.
        mid = (a+b)/2
        ref = logdensity(mid, family, shapes, loc, scale)
        def relative(t):
            x = a+(b-a)*t
            # Density-only evaluation avoids nesting the skew-CDF quadrature.
            z = (x-mp.mpf(loc))/mp.mpf(scale)
            if family == 0:
                r = 1/mp.mpf(shapes[0])
                f = r*emg_h(z, r)/scale
            elif family == 1:
                f = 2*phi(z)*cdf(mp.mpf(shapes[0])*z)/scale
            else:
                u, sigma, alpha, beta = map(mp.mpf, shapes)
                t0 = (mp.log(z)-u)/sigma
                f = alpha*beta/(alpha+beta)*(emg_h(t0, alpha*sigma)
                    +emg_h(-t0, beta*sigma))/(x-loc)
            return mp.exp(mp.log(f)-ref)
        return ref+mp.log(b-a)+mp.log(mp.quad(relative, [0, 1],
                                                    method='gauss-legendre', maxdegree=3))
    left, right = oracle(a, family, shapes, loc, scale), oracle(b, family, shapes, loc, scale)
    if right[1] <= left[2]:
        return mp.log(mp.exp(right[1])-mp.exp(left[1]))
    return mp.log(mp.exp(left[2])-mp.exp(right[2]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--module-dir', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--cases', type=int, default=20)
    parser.add_argument('--reuse', action='store_true')
    args = parser.parse_args()
    cache = ROOT/'build'/'normal_mix_reference'
    cache.mkdir(parents=True, exist_ok=True)
    cases, inputs = [], []
    rng = random.Random(2026111801)

    def add(family, shapes, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed'):
        value = sum(logdensity(x, family, shapes, loc, scale) for x in exact)
        value += sum(oracle(x, family, shapes, loc, scale)[1] for x in left)
        value += sum(oracle(x, family, shapes, loc, scale)[2] for x in right)
        value += sum(interval(a, b, family, shapes, loc, scale) for a, b in intervals)
        expected = float(value)
        assert math.isfinite(expected), (family, shapes, label, value)
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{90+family} {shapes[0]:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        inputs.extend(f'{s:.17g}' for s in shapes[1:])
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    if args.reuse:
        cases = json.loads((cache/'cases.json').read_text())
        inputs = (cache/'input.txt').read_text().splitlines()
    else:
        with mp.workdps(100):
            for family in range(3):
                for _ in range(args.cases):
                    shapes = ((10**rng.uniform(-2, 2),) if family == 0 else
                              (rng.uniform(-10, 10),) if family == 1 else
                              (rng.uniform(-1, 1), 10**rng.uniform(-1, .5),
                               10**rng.uniform(-1, 1), 10**rng.uniform(-1, 1)))
                    scale = 10**rng.uniform(-200, 200)
                    loc = scale*rng.uniform(-.3, .3)
                    point = lambda z: loc+scale*z
                    a = point(rng.uniform(.2, .6))
                    b = a+scale*10**rng.uniform(-15, -2)
                    if a == b:
                        b = math.nextafter(a, math.inf)
                    add(family, shapes, loc, scale, exact=[point(.4)], left=[point(.1)],
                        right=[point(.8)], intervals=[(a, b), (point(.2), point(.7))])
                print(FAMILIES[family], 'mixed references ready', flush=True)
            for k in (.01, .1, 1., 100.):
                for x in (-1000., -12., 0., 1000.):
                    add(0, (k,), exact=[x], left=[x], right=[x], label='EMG tails')
                    add(0, (k,), intervals=[(x, math.nextafter(x, math.inf))]
                        if x else [(.1, math.nextafter(.1, math.inf))], label='EMG adjacent')
            for a in (-10., -1., 0., 1., 10.):
                for x in (-1000., -12., -1., 0., 1., 12., 1000.):
                    add(1, (a,), exact=[x], left=[x], right=[x], label='skew tails')
                    if x:
                        add(1, (a,), intervals=[(x, math.nextafter(x, math.inf))],
                            label='skew adjacent')
                print('skew stress shape', a, 'ready', flush=True)
            for shapes in ((0., .01, .01, 100.), (0., 10., 100., .01),
                           (-100., 1., 2., 3.), (100., 1., 3., 2.),
                           (0., .01, .01, .01), (0., 10., 100., 100.)):
                for x in (1e-300, 1., 1e300):
                    add(2, shapes, exact=[x], left=[x], right=[x], label='Pareto tails')
                    add(2, shapes, intervals=[(x, math.nextafter(x, math.inf))], label='Pareto adjacent')
            for family, shapes in enumerate(((1.,), (2.,), (.1, .7, 2., 3.))):
                for scale in (1e-300, 1e300):
                    x = .4*scale
                    add(family, shapes, scale=scale, exact=[x], left=[x], right=[x],
                        intervals=[(x, math.nextafter(x, math.inf))], label='extreme scale')
            x = math.nextafter(0., 1.)
            for scale in (1., 1e300):
                add(2, (.1, .7, 2., 3.), scale=scale, exact=[x], left=[x], right=[x],
                    intervals=[(x, 2*x)], label='subnormal Pareto coordinate')
            add(0, (1.,), left=[700.], label='tiny complement')
            add(0, (1.,), right=[-35.], label='tiny complement')
            add(1, (2.,), left=[35.], label='tiny complement')
            add(1, (-2.,), right=[-35.], label='tiny complement')
            add(2, (.1, .7, 3., 2.), left=[1e100], label='tiny complement')
            add(2, (.1, .7, 2., 3.), right=[1e-100], label='tiny complement')
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
        failures, worst, worst_tiny = [], 0., 0.
        for value, (family, label, expected) in zip(actual, cases):
            error = (abs(value/expected-1) if label == 'tiny complement' else
                     abs(value-expected)/max(1., abs(expected)))
            if not math.isfinite(value) or error > 3e-11:
                failures.append((family, label, value, expected, error))
            if label == 'tiny complement':
                worst_tiny = max(worst_tiny, error)
            else:
                worst = max(worst, error)
        assert not failures, (mode, failures)
        print(f'{mode}: {len(cases)} comparisons passed; max scaled log error={worst}; '
              f'max relative tiny-log error={worst_tiny}', flush=True)
    sample = stats.CensoredData(uncensored=[.2, .4, .6], left=[.1], right=[.8],
                                interval=[[.3, .35], [.5, .55]])
    for name, shapes in zip(FAMILIES, ((1.7,), (1.,), (.1, .7, 2., 3.))):
        fit = optimize.minimize_scalar(
            lambda v: getattr(stats, name)._penalized_nnlf((*shapes, 0., v), sample),
            bounds=(.3, 3.), method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: SciPy scale={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
