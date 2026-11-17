#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent inverse-Gaussian-family censored references, using mpmath."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('invgauss', 'wald', 'recipinvgauss')


def density(x, family, shape, loc, scale):
    x, mu, loc, scale = map(mp.mpf, (x, 1. if family == 1 else shape, loc, scale))
    original = (x-loc)/scale
    z = 1/original if family == 2 else original
    penalty = (z-mu)**2/(2*mu*mu*z)
    return (-penalty-mp.log(2*mp.pi)/2
            -(mp.mpf('.5') if family == 2 else mp.mpf('1.5'))*mp.log(original)-mp.log(scale))


def logs(x, family, shape, loc, scale):
    x, mu, loc, scale = map(mp.mpf, (x, 1. if family == 1 else shape, loc, scale))
    original = (x-loc)/scale
    z = 1/original if family == 2 else original
    a = (z/mu-1)/mp.sqrt(z)
    b = (z/mu+1)/mp.sqrt(z)
    phi = -a*a/2-mp.log(2*mp.pi)/2
    # Scaled erfc evaluates the second normal term without exp(2/mu).
    if b > mp.mpf('1e154'):
        # mpmath 1.3.0's erfc series selector converts b*b to a float.
        # Use the defining Mills integral where that selector overflows.
        with mp.workdps(100):
            mills = mp.quad(lambda u: mp.exp(-u-u*u/(2*b*b)), [0, mp.inf])/b
        term = mp.exp(phi)*mills
    else:
        term = mp.exp(phi)*mp.sqrt(mp.pi/2)*mp.exp(b*b/2)*mp.erfc(b/mp.sqrt(2))
    p = mp.erfc(-a/mp.sqrt(2))/2+term
    q = mp.erfc(a/mp.sqrt(2))/2-term
    lp = mp.log(p) if p < mp.mpf('.5') else mp.log1p(-q)
    lq = mp.log(q) if q < mp.mpf('.5') else mp.log1p(-p)
    density = phi-(mp.mpf('.5') if family == 2 else mp.mpf('1.5'))*mp.log(original)-mp.log(scale)
    return (density, lq, lp) if family == 2 else (density, lp, lq)


def interval(a, b, family, shape, loc, scale):
    # Integrate a rescaled density on [0,1] for narrow intervals.
    # The scaling also avoids probability underflow in far tails.
    z = (mp.mpf(a)-loc)/scale
    mu = mp.mpf(1. if family == 1 else shape)
    bound = 2+z/(mu*mu)+1/z if family != 2 else 2+z+1/(mu*mu*z)
    if (mp.mpf(b)-a)/(mp.mpf(a)-loc)*bound < mp.mpf('.001'):
        a, b = mp.mpf(a), mp.mpf(b)
        ref = density(a, family, shape, loc, scale)
        return ref+mp.log(b-a)+mp.log(mp.quad(
            lambda t: mp.exp(density(a+t*(b-a), family, shape, loc, scale)-ref),
            [0, 1], method='gauss-legendre', maxdegree=3))
    _, lp, lq = logs(a, family, shape, loc, scale)
    _, up, uq = logs(b, family, shape, loc, scale)
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

    def add(family, shape, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed'):
        value = sum(logs(x, family, shape, loc, scale)[0] for x in exact)
        value += sum(logs(x, family, shape, loc, scale)[1] for x in left)
        value += sum(logs(x, family, shape, loc, scale)[2] for x in right)
        value += sum(interval(a, b, family, shape, loc, scale) for a, b in intervals)
        expected = float(value)
        assert math.isfinite(expected), (family, label, value)
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{75+family} {shape:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    with mp.workdps(100):
        for family in range(3):
            for _ in range(args.cases):
                shape = 10**rng.uniform(-2, 3)
                scale = 10**rng.uniform(-200, 200)
                loc = scale*rng.uniform(-3, 3)
                point = lambda p: loc+scale*p
                lo = point(rng.uniform(.05, 3))
                hi = lo+scale*10**rng.uniform(-15, -2)
                if hi == lo:
                    hi = math.nextafter(lo, math.inf)
                add(family, shape, loc, scale, exact=[point(rng.uniform(.1, 3))],
                    left=[point(.1)], right=[point(4)],
                    intervals=[(lo, hi), (point(.1), point(3))])
            print(f'{FAMILIES[family]}: mixed references ready', flush=True)
    with mp.workdps(750):
        for family in range(3):
            for shape in (1e-100, 1e-6, .01, 1., 100., 1e6, 1e100, 1e300):
                add(family, shape, exact=[.7], left=[.7], right=[.7], label='shape')
                add(family, shape, intervals=[(.7, math.nextafter(.7, math.inf))],
                    label='adjacent shape')
            for x in (1e-300, 1e-10, 144., 145., 1e10, 1e100, 1e300):
                add(family, 1.3, exact=[x], left=[x], right=[x], label='tail')
                add(family, 1.3, intervals=[(x, math.nextafter(x, math.inf))], label='adjacent tail')
            for scale in (1e-300, 1e300):
                x = .7*scale
                add(family, 1.3, scale=scale, exact=[x], left=[x], right=[x], label='scale')
                add(family, 1.3, scale=scale, intervals=[(x, math.nextafter(x, math.inf))],
                    label='adjacent scale')
            for x in (1e5, 1e10, 1e100, 1e300):
                add(family, 1e300, left=[x], right=[x], label='nearly coincident normal coordinates')
            print(f'{FAMILIES[family]}: stresses ready', flush=True)
        for family in (0, 2):
            for shape in (1e-100, 1e-300):
                center = shape if family == 0 else 1/shape
                for x in (center, math.nextafter(center, math.inf)):
                    add(family, shape, exact=[x], left=[x], right=[x],
                        label='concentrated center')
            shape = 1e-100
            for scale in (1e-100, 1e100):
                center = scale*shape if family == 0 else scale/shape
                loc = -.3*center
                x = loc+center
                add(family, shape, loc, scale, exact=[x], left=[x], right=[x],
                    label='concentrated affine center')
        shape = math.nextafter(0., 1.)
        add(0, shape, exact=[shape], left=[shape], right=[shape], label='minimum subnormal center')
    cache = ROOT/'build'/'inverse_gaussian_reference'
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
    sample = stats.CensoredData(uncensored=[.2, .4, .6], left=[.1], right=[.7],
                                interval=[[.25, .35], [.5, .55]])
    for family, name in enumerate(FAMILIES):
        distribution = getattr(stats, name)
        parameters = lambda a: (0., a) if family == 1 else (a, 0., 1.3)
        fit = optimize.minimize_scalar(lambda a: distribution._penalized_nnlf(parameters(a), sample),
                                       bounds=(.3, 5.), method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: SciPy free parameter={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
