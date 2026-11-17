#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-precision truncated-exponential, Johnson SB, and kappa3 checks."""
import argparse
import importlib.util
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('reference',
                                             ROOT/'tools/generate_censored_reference.py')
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)
FAMILIES = ('truncexpon', 'johnsonsb', 'kappa3')


def softplus(q):
    return max(q, 0)+mp.log1p(mp.exp(-abs(q)))


def logs(x, family, shape, second, loc, scale):
    """Direct arbitrary-precision mathematical density and complementary tails."""
    x, a, b, center, s = map(mp.mpf, (x, shape, second, loc, scale))
    distance = x-center
    if family == 0:
        width = mp.mpf(float(shape*scale))
        upper_gap = center+width-x
        normalizer = -mp.expm1(-width/s)
        if distance < 0 or upper_gap < 0:
            density = -mp.inf
        else:
            density = -distance/s-mp.log(s)-mp.log(normalizer)
        if distance <= 0:
            return density, -mp.inf, mp.mpf(0)
        if upper_gap <= 0:
            return density, mp.mpf(0), -mp.inf
        lp = mp.log(-mp.expm1(-distance/s))-mp.log(normalizer)
        lq = -distance/s+mp.log(-mp.expm1(-upper_gap/s))-mp.log(normalizer)
    elif family == 1:
        upper_gap = center+s-x
        if distance <= 0:
            return -mp.inf, -mp.inf, mp.mpf(0)
        if upper_gap <= 0:
            return -mp.inf, mp.mpf(0), -mp.inf
        t = mp.log(distance/upper_gap)
        q = a+b*t
        lp, lq = reference.transformed_normal_logs(q)
        density = -q*q/2-mp.log(2*mp.pi)/2+mp.log(b*s/(distance*upper_gap))
    else:
        if distance < 0:
            return -mp.inf, -mp.inf, mp.mpf(0)
        if distance == 0:
            return -mp.log(a)/a-mp.log(s), -mp.inf, mp.mpf(0)
        lz = mp.log(distance/s)
        q = mp.log(a)-a*lz
        logh = mp.log(softplus(q))-mp.log(a) if q >= -3000 else q-mp.log(a)
        # Retain the log-SF even when the lower log tail rounds to zero.
        if logh < -1000:
            lp, lq = mp.mpf(0), logh
            h = mp.mpf(0)
        else:
            h = mp.exp(logh)
            lp = -h
            lq = mp.log(-mp.expm1(-h))
        if q >= 40:
            density = -mp.log(a)/a-mp.log(s)-(1+1/a)*softplus(-q)
        else:
            density = -h-mp.log(distance)-softplus(-q)
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
        inputs.append(f'{57+family} {shape:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        if family == 1:
            inputs.append(f'{second:.17g}')
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    for family in range(3):
        for _ in range(args.cases):
            shape = rng.uniform(-3, 3) if family == 1 else 10**rng.uniform(-2, 2)
            second = 10**rng.uniform(-1, 1)
            scale = 10**rng.uniform(-100, 100)
            width = (shape if family == 0 else 1.)*scale
            loc = width*rng.uniform(-3, 3)
            exact = [loc+width*rng.uniform(.01, .95) for _ in range(3)]
            left = [loc+width*rng.uniform(.01, .4)]
            right = [loc+width*rng.uniform(.6, .99)]
            a = loc+width*rng.uniform(.01, .9)
            b = a+width*10**rng.uniform(-14, -2)
            if b == a:
                b = math.nextafter(a, math.inf)
            intervals = [(a, b), (loc+.1*width, loc+.8*width)]
            add(family, shape, second, loc, scale, exact, left, right, intervals)
    with mp.workdps(800):
        for family in range(3):
            for shape in ((1e-300, 1e-100, .5, 1., 2., 1e100, 1e300)
                          if family != 1 else (-5., 0., 5.)):
                if family == 0:
                    x = shape/2
                elif family == 1:
                    x = .3
                else:
                    x = 1.
                add(family, shape, exact=[x], left=[x], right=[x], label='shape stress')
                add(family, shape, intervals=[(x, math.nextafter(x, math.inf))],
                    label='adjacent shape stress')
            for scale in (1e-300, 1e300):
                x = scale/2
                add(family, 1., loc=0., scale=scale, exact=[x], left=[x], right=[x],
                    label='scale stress')
                add(family, 1., scale=scale,
                    intervals=[(x, math.nextafter(x, math.inf))], label='adjacent scale stress')
                add(family, 1., scale=scale,
                    intervals=[(0., math.nextafter(0., math.inf))], label='subnormal width')
        for second in (1e-300, 1e-100, 1e100, 1e300):
            add(1, 0., second, exact=[.5], label='Johnson second shape')
            add(1, 0., second, intervals=[(.5, math.nextafter(.5, math.inf))],
                label='Johnson second shape interval')
        add(1, 0., 1e300, loc=-.5, left=[1e-300], label='Johnson amplified center')
        add(1, 0., 1e300, loc=-.5, intervals=[(0., 1e-300)], label='Johnson amplified interval')
        add(2, 1e300, loc=-1., right=[1e-300], label='kappa amplified scale point')
        add(2, 1e300, loc=-1., intervals=[(0., 1e-300)], label='kappa amplified interval')
        add(2, 2., loc=1e308, scale=1e308, exact=[1.5e308], left=[1.5e308],
            right=[1.5e308], label='kappa overflowing scale point')
        add(2, 1e300, loc=1e308, scale=1e308, exact=[1.5e308],
            label='kappa large shape and overflowing scale point')
        for family in (0, 1):
            for x in (math.nextafter(0., math.inf), 1e-300, 1e-15,
                      math.nextafter(1., -math.inf)):
                add(family, 1., exact=[x], left=[x], right=[x], label='endpoint gaps')
                if x < .5:
                    add(family, 1., intervals=[(x, math.nextafter(x, math.inf))],
                        label='endpoint adjacent gaps')
            add(family, 1., intervals=[(1e-100, math.nextafter(1., -math.inf))],
                label='nearly full support')
        add(0, 1e-300, scale=1e300, exact=[.25], label='nearly uniform log density')
        add(0, 1e300, scale=1e-300, exact=[1e-300], left=[1e-300],
            right=[1e-300], label='opposing cutoff and scale')
        add(0, math.nextafter(0., math.inf), scale=1.49, exact=[0.],
            intervals=[(0., math.nextafter(0., math.inf))], label='rounded subnormal cutoff')
        smallest = math.nextafter(0., math.inf)
        add(0, 1., scale=21*smallest, exact=[10*smallest], label='odd subnormal width')
        add(1, 0., 100., scale=21*smallest, exact=[10*smallest], left=[10*smallest],
            label='odd subnormal Johnson midpoint')
    for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow']),
                        ('optimized', ['-O3', '-fcheck=bounds'])]:
        out = ROOT/'build'/'bounded_shape_reference'/mode
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
    sample = stats.CensoredData(uncensored=[.2, .4, .6], left=[.1], right=[.7],
                               interval=[[.25, .35], [.5, .55]])
    for name in FAMILIES:
        distribution = getattr(stats, name)
        fit = optimize.minimize_scalar(lambda shape: distribution._penalized_nnlf(
            (shape, 1.2, 0., 1.3) if name == 'johnsonsb' else (shape, 0., 1.3), sample),
            bounds=(-3., 3.) if name == 'johnsonsb' else (.6 if name == 'truncexpon' else .1, 3.),
            method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: fixed-support SciPy shape={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
