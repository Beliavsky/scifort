#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-precision Student-t, F, and Jones-Faddy censoring checks."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('t', 'f', 'jf_skew_t')


def coordinates(x, family, first, second, loc, scale):
    x, c, d, mu, s = map(mp.mpf, (x, first, second, loc, scale))
    z = (x-mu)/s
    if family == 0:
        a = b = mp.mpf(float(first/2))
    elif family == 1:
        a, b = mp.mpf(float(first/2)), mp.mpf(float(second/2))
    else:
        a, b = c, d
    if family == 1:
        y = (c/d)*z
        u, v = y/(1+y), 1/(1+y)
    else:
        h = mp.sqrt(a+b+z*z)
        small = (a+b)/(2*h*(h+abs(z)))
        u, v = (1-small, small) if z >= 0 else (small, 1-small)
    return a, b, u, v, x, mu, s


def density(x, family, first, second, loc, scale):
    a, b, u, v, x, mu, s = coordinates(x, family, first, second, loc, scale)
    lu = mp.log1p(-v) if v < mp.mpf('.5') else mp.log(u)
    lv = mp.log1p(-u) if u < mp.mpf('.5') else mp.log(v)
    kernel = a*lu+b*lv-mp.log(mp.beta(a, b))
    if family == 1:
        return kernel-mp.log(x-mu)
    return kernel+(lu+lv)/2+2*mp.log(2)-mp.log(a+b)/2-mp.log(s)


def logs(x, family, first, second, loc, scale):
    a, b, u, v, _, _, _ = coordinates(x, family, first, second, loc, scale)
    log_density = density(x, family, first, second, loc, scale)
    if a == b and u == v:
        return log_density, -mp.log(2), -mp.log(2)
    if b == 1:
        lp = a*mp.log(u)
        return log_density, lp, mp.log(-mp.expm1(lp))
    if a == 1:
        lq = b*mp.log(v)
        return log_density, mp.log(-mp.expm1(lq)), lq
    if a == b and a >= 10000:
        logb = mp.log(mp.beta(a, b))
        mass = mp.quad(lambda w: mp.exp((a-1)*mp.log(w)+(b-1)*mp.log1p(-w)-logb),
                       [mp.mpf('.5'), u])
        return log_density, mp.log(mp.mpf('.5')+mass), mp.log(mp.mpf('.5')-mass)
    lp = mp.log(mp.betainc(a, b, 0, u, regularized=True))
    lq = mp.log(mp.betainc(b, a, 0, v, regularized=True))
    return log_density, lp, lq


def interval(lo, hi, family, first, second, loc, scale, quadrature=False):
    if quadrature:
        start, width = mp.mpf(lo), mp.mpf(hi)-lo
        integral = mp.quad(lambda t: mp.exp(density(start+width*t, family, first,
                                                   second, loc, scale)), [0, 1])
        return mp.log(width)+mp.log(integral)
    _, pa, qa = logs(lo, family, first, second, loc, scale)
    _, pb, qb = logs(hi, family, first, second, loc, scale)
    if pb <= qa:
        return pb+mp.log(-mp.expm1(pa-pb))
    if qa <= -mp.log(2):
        return qa+mp.log(-mp.expm1(qb-qa))
    return mp.log1p(-mp.exp(pa)-mp.exp(qb))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--module-dir', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--cases', type=int, default=100)
    args = parser.parse_args()
    rng = random.Random(20261117)
    cases, inputs = [], []

    def add(family, first, second=1.7, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed', quadrature=False):
        value = sum(logs(x, family, first, second, loc, scale)[0] for x in exact)
        value += sum(logs(x, family, first, second, loc, scale)[1] for x in left)
        value += sum(logs(x, family, first, second, loc, scale)[2] for x in right)
        value += sum(interval(a, b, family, first, second, loc, scale, quadrature)
                     for a, b in intervals)
        expected = float(value)
        assert math.isfinite(expected), (family, label, value)
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{72+family} {first:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        if family:
            inputs.append(f'{second:.17g}')
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    with mp.workdps(100):
        for family in range(3):
            for _ in range(args.cases):
                a, b = 10**rng.uniform(-1, 2), 10**rng.uniform(-1, 2)
                scale = 10**rng.uniform(-200, 200)
                loc = scale*rng.uniform(-3, 3)
                point = lambda p: loc+scale*p
                low = .01 if family == 1 else -3
                x = point(rng.uniform(low, 3))
                lo = point(rng.uniform(max(low, .01), 2))
                hi = lo+scale*10**rng.uniform(-15, -2)
                if hi == lo:
                    hi = math.nextafter(lo, math.inf)
                add(family, a, b, loc, scale, exact=[x], left=[point(.1)],
                    right=[point(4)], intervals=[(lo, hi), (point(.1), point(3))])
            print(f'{FAMILIES[family]}: generated mixed references', flush=True)
    with mp.workdps(1100):
        for family in range(3):
            for shape in (1e-300, .1, 1., 10., 100.):
                add(family, shape, shape, exact=[.7], left=[.7], right=[.7],
                    label='shape stress')
                add(family, shape, shape, intervals=[(.7, 2.)], label='wide shape interval')
                add(family, shape, shape, intervals=[(.7, math.nextafter(.7, math.inf))],
                    label='adjacent shape interval')
            for scale in (1e-300, 1e300):
                x = .2*scale
                add(family, 1.3, 1.7, scale=scale, exact=[x], left=[x], right=[x],
                    label='scale stress')
                add(family, 1.3, 1.7, scale=scale, intervals=[(x, math.nextafter(x, math.inf))],
                    label='adjacent scale')
            for x in (1e-300, 1e300):
                add(family, 1.3, 1.7, exact=[x], left=[x], right=[x], label='tail stress')
                add(family, 1.3, 1.7, intervals=[(x, math.nextafter(x, math.inf))],
                    label='adjacent tail')
            a = float.fromhex('0x1.0p-1022')*(2 if family < 2 else 1)
            add(family, a, a, exact=[.7], left=[.7], right=[.7], label='minimum beta shapes')
            add(family, a, a, intervals=[(.7, 2.)], label='minimum-shape interval')
            if family:
                for first, second in ((1e-300, 1.), (1., 1e-300), (1e-300, 1e-20)):
                    add(family, first, second, exact=[.7], left=[.7], right=[.7],
                        label='asymmetric shape stress')
                    add(family, first, second, intervals=[(.7, 2.)],
                        label='asymmetric shape interval')
            if family != 1:
                add(family, a, a, intervals=[(-2., .7)], label='tiny center crossing')
                add(family, 1., 1., scale=1e300, intervals=[(-1e-300, 1e-300)],
                    label='underflowed scaled width')
            print(f'{FAMILIES[family]}: generated stress references', flush=True)
    with mp.workdps(90):
        for family in range(3):
            a = 1e6 if family < 2 else 5e5
            center = 1. if family == 1 else 0.
            add(family, a, a, exact=[center], left=[center], right=[center],
                label='maximum shape center')
            add(family, a, a, intervals=[(center, math.nextafter(center, math.inf))],
                label='maximum shape adjacent', quadrature=True)
            point = 1.001 if family == 1 else .01
            add(family, a, a, exact=[point], left=[point], right=[point],
                label='maximum shape near center')
            add(family, a, a, intervals=[(point, math.nextafter(point, math.inf))],
                label='maximum shape near interval', quadrature=True)
        for family in (1, 2):
            first, second = (1999998., 2.) if family == 1 else (999999., 1.)
            point = .7 if family == 1 else 7e5
            add(family, first, second, exact=[point], left=[point], right=[point],
                label='maximum asymmetric shape')
            add(family, first, second, intervals=[(point, math.nextafter(point, math.inf))],
                label='maximum asymmetric interval', quadrature=True)
    cache = ROOT/'build'/'beta_t_reference'
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
        worst = (0., None)
        for value, (family, label, expected) in zip(actual, cases):
            error = abs(value-expected)/max(1., abs(expected))
            assert math.isfinite(value) and error <= 2e-11, \
                (mode, family, label, value, expected, error)
            if error > worst[0]:
                worst = (error, (family, label))
        print(f'{mode}: {len(cases)} comparisons passed; max scaled error={worst}', flush=True)
    sample = stats.CensoredData(uncensored=[.2, .4, .6], left=[.1], right=[.7],
                                interval=[[.25, .35], [.5, .55]])
    for family, name in enumerate(FAMILIES):
        distribution = getattr(stats, name)
        parameters = lambda a: (a, 0., 1.3) if not family else (a, 1.7, 0., 1.3)
        fit = optimize.minimize_scalar(lambda a: distribution._penalized_nnlf(parameters(a), sample),
                                       bounds=(.3, 5.), method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: SciPy shape={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
