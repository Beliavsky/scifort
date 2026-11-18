#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-precision bounded censored likelihood references.

GH: beta/hypergeometric normalization and density quadrature.
Irwin-Hall: inclusion-exclusion polynomial (not the Fortran B-spline recurrence).
Line von Mises: Bessel normalization and density quadrature.
"""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('gausshyper', 'irwinhall', 'vonmises_line')


def oracle(family, shapes):
    s = tuple(map(mp.mpf, shapes))
    lo, hi = (mp.mpf(0), mp.mpf(1))
    if family == 0:
        a, b, c, z = s
        norm = mp.beta(a, b)*mp.hyp2f1(c, a, a+b, -z)

        def pdf(x):
            return x**(a-1)*(1-x)**(b-1)*(1+z*x)**(-c)/norm

        def mass(x, y):
            if c == 0:
                return mp.betainc(a, b, x, y, regularized=True)
            if x == 0 and y == 1:
                return mp.mpf(1)
            if x == 0:
                # Remove the left endpoint power singularity analytically.
                return y**a/(a*norm)*mp.quad(
                    lambda t: (1-y*t**(1/a))**(b-1)*(1+z*y*t**(1/a))**(-c), [0, 1])
            if y == 1:
                width = 1-x
                return width**b/(b*norm)*mp.quad(
                    lambda t: (1-width*t**(1/b))**(a-1)
                    *(1+z*(1-width*t**(1/b)))**(-c), [0, 1])
            # Factor the density and width so absolute quadrature stopping
            # criteria remain meaningful for probabilities below 1e-300.
            ref = pdf((x+y)/2)
            return ref*(y-x)*mp.quad(
                lambda t: pdf(x+(y-x)*t)/ref, [0, mp.mpf('.5'), 1])
    elif family == 1:
        n = int(s[0])
        hi = mp.mpf(n)

        def pdf(x):
            return mp.fsum((-1)**k*mp.binomial(n, k)*(x-k)**(n-1)
                           for k in range(min(n, int(mp.floor(x)))+1))/mp.factorial(n-1)

        def cdf(x):
            if x <= 0:
                return mp.mpf(0)
            if x >= n:
                return mp.mpf(1)
            if x > mp.mpf(n)/2:
                return 1-cdf(n-x)
            return mp.fsum((-1)**k*mp.binomial(n, k)*(x-k)**n
                           for k in range(int(mp.floor(x))+1))/mp.factorial(n)

        def mass(x, y):
            if x > mp.mpf(n)/2:
                return cdf(n-x)-cdf(n-y)
            return cdf(y)-cdf(x)
    else:
        kappa = s[0]
        # The finite support contract uses SciFort's rounded binary64 pi.
        lo, hi = -mp.mpf(math.pi), mp.mpf(math.pi)
        norm = 2*mp.pi*mp.besseli(0, kappa)

        def pdf(x):
            return mp.exp(kappa*mp.cos(x))/norm

        def mass(x, y):
            return mp.quad(pdf, [x, (x+y)/2, y])

    return lo, hi, pdf, mass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--module-dir', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--cases', type=int, default=20, help='random mixed cases per family')
    parser.add_argument('--reuse', action='store_true', help='reuse build/ reference cache')
    args = parser.parse_args()
    cache = ROOT/'build'/'bounded_final_reference'
    cache.mkdir(parents=True, exist_ok=True)
    cases, inputs = [], []
    rng = random.Random(20261118)

    def add(family, shapes, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed'):
        lo, hi, pdf, mass = oracle(family, shapes)
        mu, sigma = mp.mpf(loc), mp.mpf(scale)
        coord = lambda x: (mp.mpf(x)-mu)/sigma
        value = mp.fsum(mp.log(pdf(coord(x)))-mp.log(sigma) for x in exact)
        value += mp.fsum(mp.log(mass(lo, coord(x))) for x in left)
        value += mp.fsum(mp.log(mass(coord(x), hi)) for x in right)
        value += mp.fsum(mp.log(mass(coord(a), coord(b))) for a, b in intervals)
        expected = float(value)
        assert math.isfinite(expected), (family, shapes, label, value)
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{87+family} {shapes[0]:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        inputs.extend(f'{s:.17g}' for s in shapes[1:])
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    if args.reuse:
        cases = json.loads((cache/'cases.json').read_text())
        inputs = (cache/'input.txt').read_text().splitlines()
    else:
        with mp.workdps(180):
            for family in range(3):
                for _ in range(args.cases):
                    shapes = ((rng.uniform(.2, 10), rng.uniform(.2, 10),
                               rng.uniform(-10, 10), rng.uniform(-.5, 10)) if family == 0 else
                              (float(rng.randint(1, 32)),) if family == 1 else
                              (rng.uniform(0, 100),))
                    scale = 10**rng.uniform(-200, 200)
                    loc = scale*rng.uniform(-.3, .3)
                    point = lambda x: loc+scale*x
                    a = point(rng.uniform(.2, .6))
                    b = a+scale*10**rng.uniform(-15, -2)
                    if b == a:
                        b = math.nextafter(a, math.inf)
                    add(family, shapes, loc, scale, exact=[point(.4)], left=[point(.1)],
                        right=[point(.8)], intervals=[(a, b), (point(.2), point(.7))])
                print(f'{FAMILIES[family]}: mixed references ready', flush=True)
            for a in (.1, 1., 20.):
                for b in (.1, 1., 20.):
                    for c, z in ((-10., -.5), (10., 10.), (0., 1.)):
                        add(0, (a, b, c, z), exact=[.4], left=[1e-100], right=[.999],
                            intervals=[(.4, math.nextafter(.4, math.inf))], label='shape limits')
            x = math.nextafter(1., 0.)
            for shapes in ((1., .1, 0., 1.), (.1, 20., 10., 10.), (20., .1, -10., -.5)):
                add(0, shapes, right=[x], intervals=[(math.nextafter(x, 0.), x)],
                    label='last upper endpoint interval')
            for n in (1., 2., 3., 16., 32.):
                for x in (1e-200, .4, n-.1):
                    add(1, (n,), exact=[x], left=[x], right=[x], label='sum tails')
                    add(1, (n,), intervals=[(x, math.nextafter(x, math.inf))], label='adjacent sum')
                x = math.nextafter(0., 1.)
                add(1, (n,), left=[x], intervals=[(x, 2*x)], label='subnormal sum')
                for fraction in (.3, .5, .7):
                    x = n*fraction
                    add(1, (n,), exact=[x], left=[x], right=[x],
                        intervals=[(.49*n, .51*n)], label='interior polynomial pieces')
            for kappa in (0., 1e-10, 1., 100.):
                for x in (-math.pi+1e-10, 0., math.pi-1e-10):
                    add(2, (kappa,), exact=[x], left=[x], right=[x], label='line tails')
                    if x != 0:
                        add(2, (kappa,), intervals=[(x, math.nextafter(x, math.inf))],
                            label='adjacent line')
            x = math.nextafter(32., 0.)
            add(1, (32.,), right=[x], intervals=[(math.nextafter(x, 0.), x)],
                label='last upper spline interval')
            for family, shapes in enumerate(((2., 3., 1., 2.), (32.,), (100.,))):
                for scale in (1e-300, 1e300):
                    x = .4*scale
                    add(family, shapes, scale=scale, exact=[x], left=[x], right=[x],
                        intervals=[(x, math.nextafter(x, math.inf))], label='extreme scale')
        (cache/'cases.json').write_text(json.dumps(cases), encoding='ascii')
        (cache/'input.txt').write_text('\n'.join(inputs)+'\n', encoding='ascii')

    for mode, flags in [('checked', ['-O0', '-g', '-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow']),
                        ('optimized', ['-O3', '-fcheck=bounds'])]:
        out = cache/mode
        out.mkdir(parents=True, exist_ok=True)
        cmd = [args.compiler, *flags, '-J'+str(out), '-I'+str(out),
               '-I'+str(args.module_dir), 'src/stats/gausshyper.f90',
               'src/special/log_gamma.f90', 'src/empirical_distribution.f90',
               'src/censored.f90', 'validation/censored_reference_driver.f90',
               str(args.library), '-o', str(out/'driver.exe')]
        print(subprocess.list2cmdline(cmd), flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)
        result = subprocess.run([str(out/'driver.exe')], input='\n'.join(inputs)+'\n',
                                text=True, capture_output=True)
        if result.returncode:
            n = len(result.stdout.split())
            raise RuntimeError((n, cases[n:n+1], result.stderr))
        actual = [float(s) for s in result.stdout.split()]
        assert len(actual) == len(cases)
        failures, worst = [], 0.
        for value, (family, label, expected) in zip(actual, cases):
            error = abs(value-expected)/max(1., abs(expected))
            if not math.isfinite(value) or error > 3e-11:
                failures.append((family, label, value, expected, error))
            worst = max(worst, error)
        assert not failures, (mode, failures)
        print(f'{mode}: {len(cases)} comparisons passed; max scaled error={worst}', flush=True)
    sample = stats.CensoredData(uncensored=[.2, .4, .6], left=[.1], right=[.8],
                                interval=[[.3, .35], [.5, .55]])
    for name, shapes in zip(FAMILIES, ((1.7, 2.7, .5, 1.), (3.,), (2.,))):
        fit = optimize.minimize_scalar(
            lambda v: getattr(stats, name)._penalized_nnlf((*shapes, 0., v), sample),
            bounds=(.3, 3.), method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: SciPy scale={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
