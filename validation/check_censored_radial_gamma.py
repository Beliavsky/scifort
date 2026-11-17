#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent mpmath references for chi2, chi, and Nakagami censoring."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('chi2', 'chi', 'nakagami')


def logs(x, family, shape, loc, scale):
    x, c, mu, s = map(mp.mpf, (x, shape, loc, scale))
    a = mp.mpf(float(shape/2)) if family < 2 else c
    power = 1 if family == 0 else 2
    rate = mp.mpf('.5') if family < 2 else c
    z = (x-mu)/s
    if z <= 0:
        return -mp.inf, -mp.inf, mp.mpf(0)
    y = rate*z**power
    density = mp.log(power)+a*mp.log(rate)-mp.loggamma(a) \
        +(power*a-1)*mp.log(z)-y-mp.log(s)
    if y <= mp.mpf('1.5'):
        # DLMF 8.7: alternating upper-tail series, retaining tiny shapes.
        base = a*mp.log(y)-mp.loggamma(1+a)
        term = -y
        total = term/(a+1)
        k = 1
        while True:
            k += 1
            term *= -y/k
            increment = term/(a+k)
            total += increment
            if abs(increment) < mp.eps*abs(total):
                break
        q = -mp.expm1(base)-mp.exp(base)*a*total
        lq = mp.log(q)
        lp = mp.log(mp.gammainc(a, 0, y))-mp.loggamma(a)
    else:
        lp = mp.log(mp.gammainc(a, 0, y))-mp.loggamma(a)
        lq = mp.log(mp.gammainc(a, y, mp.inf))-mp.loggamma(a)
    return density, lp, lq


def interval(lo, hi, family, shape, loc, scale):
    if lo <= loc:
        return logs(hi, family, shape, loc, scale)[1]
    _, pa, qa = logs(lo, family, shape, loc, scale)
    _, pb, qb = logs(hi, family, shape, loc, scale)
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

    def add(family, shape, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed'):
        value = sum(logs(x, family, shape, loc, scale)[0] for x in exact)
        value += sum(logs(x, family, shape, loc, scale)[1] for x in left)
        value += sum(logs(x, family, shape, loc, scale)[2] for x in right)
        value += sum(interval(a, b, family, shape, loc, scale) for a, b in intervals)
        expected = float(value)
        assert math.isfinite(expected), (family, label, value)
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{69+family} {shape:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    with mp.workdps(100):
        for family in range(3):
            for _ in range(args.cases):
                shape = 10**rng.uniform(-1, 2)
                scale = 10**rng.uniform(-200, 200)
                loc = scale*rng.uniform(-3, 3)
                point = lambda p: loc+scale*p
                x = point(rng.uniform(.01, 3))
                lo = point(rng.uniform(.01, 2))
                hi = lo+scale*10**rng.uniform(-15, -2)
                if hi == lo:
                    hi = math.nextafter(lo, math.inf)
                add(family, shape, loc, scale, exact=[x], left=[point(.1)],
                    right=[point(4)], intervals=[(lo, hi), (point(.1), point(3))])
            print(f'{FAMILIES[family]}: generated random references', flush=True)
    with mp.workdps(1100):
        for family in range(3):
            for shape in (1e-300, .1, 1., 2., 10., 100.):
                x = .7
                add(family, shape, exact=[x], left=[x], right=[x], label='shape stress')
                add(family, shape, intervals=[(x, math.nextafter(x, math.inf))],
                    label='adjacent shape')
            for scale in (1e-300, 1e300):
                x = .2*scale
                add(family, 1.3, scale=scale, exact=[x], left=[x], right=[x],
                    label='scale stress')
                add(family, 1.3, scale=scale, intervals=[(x, math.nextafter(x, math.inf))],
                    label='adjacent scale')
            tiny = math.nextafter(0., math.inf)
            c = 2*tiny if family < 2 else tiny
            add(family, c, exact=[.7], left=[.7], right=[.7], label='minimum shape')
            add(family, 1., intervals=[(0., tiny)], label='subnormal interval')
            add(family, 1., scale=1e300, intervals=[(0., tiny)],
                label='underflowed powered distance')
            add(family, .7, loc=-1e308, scale=1e308, exact=[1e308],
                left=[1e308], right=[1e308], label='overflowing physical distance')
            add(family, 1.3, intervals=[(.01, 100.)], label='nearly full support')
            print(f'{FAMILIES[family]}: generated references', flush=True)
    with mp.workdps(90):
        for family in range(3):
            shape = 2e6 if family < 2 else 1e6
            center = (2e6 if family == 0 else math.sqrt(2e6)) if family < 2 else 1.
            add(family, shape, exact=[center], left=[center], right=[center],
                label='maximum shape center')
            add(family, shape, intervals=[(center, math.nextafter(center, math.inf))],
                label='maximum shape adjacent')

    cache = ROOT/'build'/'radial_gamma_reference'
    cache.mkdir(parents=True, exist_ok=True)
    (cache/'cases.json').write_text(json.dumps(cases), encoding='ascii')
    (cache/'input.txt').write_text('\n'.join(inputs)+'\n', encoding='ascii')
    for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow']),
                        ('optimized', ['-O3', '-fcheck=bounds'])]:
        out = ROOT/'build'/'radial_gamma_reference'/mode
        out.mkdir(parents=True, exist_ok=True)
        (out/'input.txt').write_text('\n'.join(inputs)+'\n', encoding='ascii')
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
            completed = len(result.stdout.split())
            raise RuntimeError((completed, cases[completed:completed+1], result.stderr))
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
    for name in FAMILIES:
        distribution = getattr(stats, name)
        fit = optimize.minimize_scalar(
            lambda c: distribution._penalized_nnlf((c, 0., 1.3), sample),
            bounds=(.3, 5.), method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: SciPy shape={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
