#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent compact-support censoring checks from public density specifications.

Uses an integrated hypergeometric sine series, incomplete beta, and arbitrary
precision arithmetic. The rounded physical half-width matches the documented
native likelihood contract. No upstream implementation is reproduced.
"""
import argparse
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('cosine', 'anglit', 'semicircular')
HALF_WIDTHS = (math.pi, math.pi/4, 1.)


def small_tail(q, family):
    if q <= 0:
        return mp.mpf(0)
    if family == 0:
        return mp.pi**2*q**3/12*mp.hyper([1], [2, mp.mpf('2.5')], -mp.pi**2*q*q/4)
    if family == 1:
        return mp.sin(mp.pi*q/4)**2
    return mp.betainc(mp.mpf('1.5'), mp.mpf('1.5'), 0, q/2, regularized=True)


def tail(x, family, loc, radius, upper=False):
    # Form the endpoint gap before dividing; do not subtract a rounded unit coordinate.
    delta = mp.mpf(x)-loc
    q = (radius-abs(delta))/radius
    small = small_tail(q, family)
    near_upper = delta >= 0
    return small if upper == near_upper else 1-small


def log_tail(x, family, loc, radius, upper=False):
    delta = mp.mpf(x)-loc
    q = (radius-abs(delta))/radius
    small = small_tail(q, family)
    return mp.log(small) if upper == (delta >= 0) else mp.log1p(-small)


def logpdf(x, family, loc, radius):
    delta = mp.mpf(x)-loc
    q = (radius-abs(delta))/radius
    if q <= 0:
        return -mp.inf
    if family == 0:
        density = mp.sin(mp.pi*q/2)**2
    elif family == 1:
        density = mp.pi/4*mp.sin(mp.pi*q/2)
    else:
        density = 2/mp.pi*mp.sqrt(q*(2-q))
    return mp.log(density)-mp.log(radius)


def loginterval(a, b, family, loc, radius):
    width = mp.mpf(b)-a
    if width/radius < mp.mpf('1e-30') and abs(mp.mpf(a)-loc) < radius/2:
        # On this central interval, relative density variation is < 1e-30.
        # Integrate the positive density instead of subtracting unit-sized tails.
        return mp.log(width)+logpdf((mp.mpf(a)+b)/2, family, loc, radius)
    if a <= loc <= b:
        omitted = tail(a, family, loc, radius)+tail(b, family, loc, radius, True)
        return mp.log1p(-omitted)
    if b <= loc:
        mass = tail(b, family, loc, radius)-tail(a, family, loc, radius)
    else:
        mass = tail(a, family, loc, radius, True)-tail(b, family, loc, radius, True)
    return mp.log(mass)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module-dir', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--seed', type=int, default=20261116)
    parser.add_argument('--cases', type=int, default=100)
    args = parser.parse_args()
    mp.mp.dps = 100
    rng = random.Random(args.seed)
    inputs, cases = [], []

    def add(family, loc, scale, exact=(), left=(), right=(), intervals=(), label='mixed'):
        radius = mp.mpf(HALF_WIDTHS[family]*scale)
        center = mp.mpf(loc)
        value = sum(logpdf(x, family, center, radius) for x in exact)
        value += sum(log_tail(x, family, center, radius) for x in left)
        value += sum(log_tail(x, family, center, radius, True) for x in right)
        value += sum(loginterval(a, b, family, center, radius) for a, b in intervals)
        expected = float(value)
        if not math.isfinite(expected):
            raise AssertionError((FAMILIES[family], label, value))
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{54+family} 1 {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    for family in range(3):
        for _ in range(args.cases):
            scale = 10**rng.uniform(-100, 100)
            radius = HALF_WIDTHS[family]*scale
            loc = radius*rng.uniform(-3, 3)
            exact = [loc+radius*rng.uniform(-.95, .95) for _ in range(3)]
            left = [loc+radius*rng.uniform(-.99, .99)]
            right = [loc+radius*rng.uniform(-.99, .99)]
            a = loc+radius*rng.uniform(-.9, .9)
            b = a+radius*10**rng.uniform(-14, -2)
            if b == a:
                b = math.nextafter(a, math.inf)
            intervals = [(a, b), (loc-.3*radius, loc+.7*radius)]
            add(family, loc, scale, exact, left, right, intervals)
        # Approach both endpoints, and sample across each analytic-series switch.
        radius = HALF_WIDTHS[family]
        for q in (1e-1, 1e-4, 1e-8, 1e-12, 2**-50, .25, .5, .75):
            x = radius*(1-q)
            add(family, 0., 1., exact=[x], right=[x], label='endpoint density/tail')
            add(family, 0., 1., intervals=[(math.nextafter(x, -math.inf), x)],
                label='endpoint adjacent interval')
            add(family, 0., 1., intervals=[(-x, x)], label='nearly full interval')
        for scale in (1e-300, 1e300):
            radius = HALF_WIDTHS[family]*scale
            add(family, 0., scale, intervals=[(0., math.nextafter(0., math.inf))],
                label='subnormal interval')
            add(family, 0., scale, exact=[radius/2], left=[-radius/2], right=[radius/2],
                label='extreme scale')
        for halfspan in (1e-5, 1e-8, 1e-14):
            add(family, 0., 1., intervals=[(-halfspan, halfspan)],
                label='narrow interval crossing center')
        # Cancellation leaves a subnormal positive endpoint gap with unit radius.
        # 800 digits also resolve cosine's cubic mass when the gap is subnormal.
        with mp.workdps(800):
            scale = 1./HALF_WIDTHS[family]
            radius = HALF_WIDTHS[family]*scale
            add(family, radius, scale, exact=[math.nextafter(0., math.inf)],
                left=[math.nextafter(0., math.inf)], label='subnormal endpoint gap')
    for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow']),
                        ('optimized', ['-O3', '-fcheck=bounds'])]:
        out = ROOT/'build'/'compact_reference'/mode
        out.mkdir(parents=True, exist_ok=True)
        cmd = [args.compiler, *flags, '-J'+str(out), '-I'+str(out),
               '-I'+str(args.module_dir), 'src/special/log_gamma.f90',
               'src/empirical_distribution.f90', 'src/censored.f90',
               'validation/censored_reference_driver.f90', str(args.library),
               '-o', str(out/'driver.exe')]
        print(subprocess.list2cmdline(cmd), flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)
        run = subprocess.run([str(out/'driver.exe')], input='\n'.join(inputs)+'\n',
                             text=True, capture_output=True, check=True)
        actual = [float(s) for s in run.stdout.split()]
        assert len(actual) == len(cases), (len(actual), len(cases))
        worst = (0., None)
        for value, (family, label, expected) in zip(actual, cases):
            error = abs(value-expected)/max(1., abs(expected))
            if not math.isfinite(value) or error > 2e-11:
                raise AssertionError((mode, family, label, value, expected, error))
            if error > worst[0]:
                worst = (error, (family, label))
        print(f'{mode}: {len(cases)} comparisons passed; max scaled error={worst}', flush=True)
    # Independent SciPy optimization with fixed center; record scale benchmarks.
    sample = stats.CensoredData(uncensored=[-.2, .1, .35], left=[-.4], right=[.45],
                               interval=[[-.3, -.1], [.2, .25]])
    for name in FAMILIES:
        distribution = getattr(stats, name)
        fit = optimize.minimize_scalar(lambda scale: distribution._penalized_nnlf(
            (0., scale), sample), bounds=(.7, 2.), method='bounded',
            options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: fixed-center SciPy scale={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
