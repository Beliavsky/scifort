#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-precision log-uniform, truncated Pareto, and generalized half-logistic checks."""
import argparse
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('loguniform', 'truncpareto', 'genhalflogistic')


def softplus(q):
    return max(q, 0)+mp.log1p(mp.exp(-abs(q)))


def logs(x, family, shape, second, loc, scale):
    """Direct arbitrary-precision mathematical density and complementary tails."""
    x, a, b, center, scale = map(mp.mpf, (x, shape, second, loc, scale))
    t = x-center
    if family < 2:
        lo = mp.mpf(float(shape*float(scale))) if family == 0 else scale
        hi = mp.mpf(float(second*float(scale)))
        power = mp.mpf(0) if family == 0 else a
        d = mp.log(hi/lo)
        if t < lo:
            return -mp.inf, -mp.inf, mp.mpf(0)
        if t > hi:
            return -mp.inf, mp.mpf(0), -mp.inf
        v, w = mp.log(t/lo), mp.log(hi/t)
        if not power:
            density = -mp.log(t)-mp.log(d)
            lp, lq = mp.log(v/d), mp.log(w/d)
        elif power > 0:
            normalizer = -mp.expm1(-power*d)
            density = mp.log(power)-mp.log(t)-power*v-mp.log(normalizer)
            lp = mp.log(-mp.expm1(-power*v))-mp.log(normalizer)
            lq = -power*v+mp.log(-mp.expm1(-power*w))-mp.log(normalizer)
        else:
            k = -power
            normalizer = -mp.expm1(-k*d)
            density = mp.log(k)-mp.log(t)-k*w-mp.log(normalizer)
            lp = -k*w+mp.log(-mp.expm1(-k*v))-mp.log(normalizer)
            lq = mp.log(-mp.expm1(-k*w))-mp.log(normalizer)
    else:
        width = mp.mpf(float(float(scale)/shape))
        w = (center+width-x)/width
        if t < 0:
            return -mp.inf, -mp.inf, mp.mpf(0)
        if w < 0:
            return -mp.inf, mp.mpf(0), -mp.inf
        if w == 0:
            density = -mp.inf if a < 1 else (mp.inf if a > 1 else mp.log(2/width))
            return density, mp.mpf(0), -mp.inf
        h = -mp.log(w)/a
        density = mp.log(2/a/width)-mp.log(w)-h-2*mp.log1p(mp.exp(-h))
        lp = mp.log(-mp.expm1(-h))-mp.log1p(mp.exp(-h))
        lq = mp.log(2)-h-mp.log1p(mp.exp(-h))
    if lp <= lq:
        lq = mp.log1p(-mp.exp(lp))
    else:
        lp = mp.log1p(-mp.exp(lq))
    return density, lp, lq


def interval(a, b, family, shape, second, loc, scale):
    _, pa, qa = logs(a, family, shape, second, loc, scale)
    _, pb, qb = logs(b, family, shape, second, loc, scale)
    if pb <= qa:
        return pb+mp.log(-mp.expm1(pa-pb))
    if qa <= -mp.log(2):
        return qa+mp.log(-mp.expm1(qb-qa))
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
        inputs.append(f'{60+family} {shape:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        if family < 2:
            inputs.append(f'{second:.17g}')
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    for family in range(3):
        for _ in range(args.cases):
            shape = 10**rng.uniform(-2, 0) if family == 0 else 10**rng.uniform(-2, 2)
            if family == 1 and rng.random() < .5:
                shape = -shape
            second = 10**rng.uniform(.1, 2)
            scale = 10**rng.uniform(-100, 100)
            lo = shape*scale if family == 0 else (scale if family == 1 else 0.)
            hi = second*scale if family < 2 else scale/shape
            width = hi-lo
            loc = width*rng.uniform(-3, 3)
            point = lambda p: loc+lo+width*p
            exact = [point(rng.uniform(.01, .95)) for _ in range(3)]
            left = [point(rng.uniform(.01, .4))]
            right = [point(rng.uniform(.6, .99))]
            a = point(rng.uniform(.01, .9))
            b = a+width*10**rng.uniform(-14, -2)
            if b == a:
                b = math.nextafter(a, math.inf)
            add(family, shape, second, loc, scale, exact, left, right,
                [(a, b), (point(.1), point(.8))])
    with mp.workdps(1100):
        for family in range(3):
            for shape in ((1e-300, 1e-100, .5) if family == 0 else
                          (-1e300, -1e-300, -.5, .5, 1e-300, 1e300) if family == 1 else
                          (1e-300, 1e-100, .5, 1., 2., 1e100, 1e300)):
                x = 1.5 if family < 2 else .5/shape
                second = 2.
                add(family, shape, second, exact=[x], left=[x], right=[x], label='shape stress')
                add(family, shape, second, intervals=[(x, math.nextafter(x, math.inf))],
                    label='adjacent shape stress')
            for scale in (1e-300, 1e300):
                shape = .5 if family == 0 else 1.
                x = scale if family == 0 else (1.5*scale if family == 1 else .5*scale)
                add(family, shape, 2., scale=scale, exact=[x], left=[x], right=[x],
                    label='scale stress')
                add(family, shape, 2., scale=scale,
                    intervals=[(x, math.nextafter(x, math.inf))], label='adjacent scale stress')
        tiny = math.nextafter(0., math.inf)
        add(0, 1., math.nextafter(1., math.inf), exact=[1.],
            label='adjacent cutoffs endpoint density')
        add(1, -1e300, math.nextafter(1., math.inf), exact=[math.nextafter(1., math.inf)],
            label='adjacent cutoffs negative power')
        add(2, tiny, scale=tiny, exact=[tiny], left=[tiny], right=[tiny],
            intervals=[(tiny, 2*tiny)], label='smallest shape and scale')
        for family in (0, 1):
            add(family, .5 if family == 0 else 1., 2., scale=21*tiny,
                exact=[22*tiny], intervals=[(22*tiny, 23*tiny)], label='rounded subnormal cutoffs')
        for x in (tiny, 1e-300, 1e-15, math.nextafter(1., -math.inf)):
            add(2, 1., exact=[x], left=[x], right=[x], label='endpoint gaps')
            add(2, 1., intervals=[(x, math.nextafter(x, math.inf))], label='endpoint intervals')
        for family in range(3):
            a = math.nextafter(1., math.inf) if family == 1 else (math.nextafter(.5, math.inf) if family == 0 else tiny)
            b = math.nextafter(2. if family < 2 else 1., -math.inf)
            add(family, .5 if family == 0 else 1., 2., intervals=[(a,b)], label='nearly full support')
    for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow']),
                        ('optimized', ['-O3', '-fcheck=bounds'])]:
        out = ROOT/'build'/'bounded_power_reference'/mode
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
    sample = stats.CensoredData(uncensored=[1.2, 1.4, 1.6], left=[1.1], right=[1.7],
                               interval=[[1.25, 1.35], [1.5, 1.55]])
    for name in FAMILIES:
        distribution = getattr(stats, name)
        params = lambda a: (a, 3., 0., 1.) if name != 'genhalflogistic' else (a, 0., 1.)
        bounds = (.1, 1.09) if name == 'loguniform' else ((.01, 5.) if name == 'truncpareto' else (.05, .58))
        fit = optimize.minimize_scalar(lambda a: distribution._penalized_nnlf(params(a), sample),
            bounds=bounds, method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: SciPy shape={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
