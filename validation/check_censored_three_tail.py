#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-precision Kolmogorov, wrapped-Cauchy and Breit-Wigner checks."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess
import warnings

import mpmath as mp
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('kstwobign', 'wrapcauchy', 'rel_breitwigner')


def kt(z):
    if z < mp.mpf('.82'):
        a = mp.pi**2/(8*z*z)
        terms = [mp.exp(-((2*j-1)**2-1)*a) for j in range(1, 25)]
        lp = mp.log(2*mp.pi)/2-mp.log(z)-a+mp.log(sum(terms))
        density = sum(((2*a*(2*j-1)**2-1)*v for j, v in enumerate(terms, 1)))
        ld = mp.log(2*mp.pi)/2-2*mp.log(z)-a+mp.log(density)
        lq = mp.log1p(-mp.exp(lp)) if lp > -10000 else mp.mpf(0)
    else:
        a = 2*z*z
        terms = [(-1)**(j-1)*mp.exp(-(j*j-1)*a) for j in range(1, 25)]
        lq = mp.log(2)-a+mp.log(sum(terms))
        ld = mp.log(8*z)-a+mp.log(sum(j*j*v for j, v in enumerate(terms, 1)))
        lp = mp.log1p(-mp.exp(lq)) if lq > -10000 else mp.mpf(0)
    return ld, lp, lq


def logdensity(x, family, shapes, loc, scale):
    x, loc, scale = map(mp.mpf, (x, loc, scale))
    z = (x-loc)/scale
    if z < 0:
        return -mp.inf
    if family == 0:
        return kt(z)[0]-mp.log(scale) if z else -mp.inf
    if family == 1:
        c = mp.mpf(shapes[0])
        # The public finite branch uses its computed binary64 upper endpoint.
        upper = mp.mpf(float(loc)+2*math.pi*float(scale))
        if x > upper:
            return -mp.inf
        d = min(x-loc, upper-x)/scale
        return mp.log(1-c*c)-mp.log(2*mp.pi)-mp.log((1-c)**2+4*c*mp.sin(d/2)**2)-mp.log(scale)
    rho = mp.mpf(shapes[0])
    norm = 2*mp.sqrt(2)*rho*rho*mp.sqrt(rho*rho+1)/(
        mp.pi*mp.sqrt(rho*rho+rho*mp.sqrt(rho*rho+1)))
    return mp.log(norm)-mp.log((z*z-rho*rho)**2+rho*rho)-mp.log(scale)


def oracle(x, family, shapes, loc, scale):
    x, loc, scale = map(mp.mpf, (x, loc, scale))
    z = (x-loc)/scale
    if z <= 0:
        return logdensity(x, family, shapes, loc, scale), -mp.inf, mp.mpf(0)
    if family == 0:
        ld, lp, lq = kt(z)
        return ld-mp.log(scale), lp, lq
    if family == 1:
        upper = mp.mpf(float(loc)+2*math.pi*float(scale))
        if x >= upper:
            return logdensity(x, family, shapes, loc, scale), mp.mpf(0), -mp.inf
        c = mp.mpf(shapes[0])
        d = min(x-loc, upper-x)/scale
        small = mp.atan((1+c)/(1-c)*mp.tan(d/2))/mp.pi
        if x-loc <= upper-x:
            lp, lq = mp.log(small), mp.log1p(-small)
        else:
            lp, lq = mp.log1p(-small), mp.log(small)
        return logdensity(x, family, shapes, loc, scale), lp, lq
    rho = mp.mpf(shapes[0])
    ld = logdensity(x, family, shapes, loc, scale)
    if z <= rho:
        f0 = logdensity(loc, family, shapes, loc, scale)+mp.log(scale)
        def relative(t):
            return mp.exp(logdensity(loc+scale*z*t, family, shapes, loc, scale)+mp.log(scale)-f0)
        lp = f0+mp.log(z)+mp.log(mp.quad(relative, [0, .5, .9, .99, 1]))
        lq = mp.log1p(-mp.exp(lp))
    else:
        f0 = logdensity(loc, family, shapes, loc, scale)+mp.log(scale)
        logk = f0+mp.log(rho**4+rho*rho)
        def reciprocal(t):
            return t*t/((1-(rho/z)**2*t*t)**2+(rho/z**2)**2*t**4)
        lq = logk-3*mp.log(z)+mp.log(mp.quad(reciprocal, [0, .5, .9, .99, 1]))
        lp = mp.log1p(-mp.exp(lq))
    return ld, lp, lq


def interval(a, b, family, shapes, loc, scale):
    a, b, loc, scale = map(mp.mpf, (a, b, loc, scale))
    if family == 0 and (b-loc)/scale <= mp.mpf('.1'):
        left, right = kt((a-loc)/scale)[1], kt((b-loc)/scale)[1]
        return right+mp.log(-mp.expm1(left-right))
    if family == 0 and (a-loc)/scale >= 4:
        left, right = kt((a-loc)/scale)[2], kt((b-loc)/scale)[2]
        return left+mp.log(-mp.expm1(right-left))
    ref = logdensity((a+b)/2, family, shapes, loc, scale)
    def relative(t):
        return mp.exp(logdensity(a+(b-a)*t, family, shapes, loc, scale)-ref)
    return ref+mp.log(b-a)+mp.log(mp.quad(relative, [0, 1], method='gauss-legendre', maxdegree=6))


def scipy_stress(driver):
    """Use multiprecision to adjudicate discrepancies with SciPy."""
    rng = random.Random(2026111805)
    inputs, refs = [], []
    for family, name in enumerate(FAMILIES):
        for _ in range(200):
            shapes = (() if family == 0 else (rng.uniform(.001, .99),) if family == 1
                      else (10**rng.uniform(-1, 2),))
            if family == 0:
                x = 10**rng.uniform(-3, 1)
            elif family == 1:
                x = 10**rng.uniform(-100, math.log10(2*math.pi-1e-10))
            else:
                x = 10**rng.uniform(-100, 100)
            dist = getattr(stats, name)(*shapes)
            for mode, method in enumerate(('logpdf', 'logcdf', 'logsf')):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    expected = float(getattr(dist, method)(x))
                if not math.isfinite(expected) or (mode and abs(expected) < 1e-7):
                    continue
                counts = [0, 0, 0]
                counts[mode] = 1
                inputs.append(f'{96+family} {shapes[0] if shapes else 0.:.17g} 0 1 '
                              + ' '.join(map(str, counts))+' 0')
                inputs.extend(f'{v:.17g}' for v in shapes[1:])
                inputs.append(f'{x:.17g}')
                refs.append((family, shapes, x, mode, expected))
    result = subprocess.run([str(driver)], input='\n'.join(inputs)+'\n',
                            text=True, capture_output=True, check=True)
    actual = list(map(float, result.stdout.split()))
    assert len(actual) == len(refs)
    adjudicated = 0
    for value, (family, shapes, x, mode, expected) in zip(actual, refs):
        error = abs(value-expected)/max(1., abs(expected))
        if not math.isfinite(value) or error > 3e-10:
            with mp.workdps(110):
                expected = float(logdensity(x, family, shapes, 0., 1.) if mode == 0 else
                                 oracle(x, family, shapes, 0., 1.)[mode])
            error = abs(value-expected)/max(1., abs(expected))
            assert math.isfinite(value) and error <= 3e-11, (
                FAMILIES[family], shapes, x, mode, value, expected, error)
            adjudicated += 1
    print(f'SciPy stress: {len(refs)} comparisons passed; '
          f'{adjudicated} discrepancies resolved against multiprecision', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--module-dir', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--cases', type=int, default=10)
    parser.add_argument('--reuse', action='store_true')
    args = parser.parse_args()
    cache = ROOT/'build'/'three_tail_reference'
    cache.mkdir(parents=True, exist_ok=True)
    cases, inputs = [], []
    rng = random.Random(2026111804)

    def add(family, shapes, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed'):
        value = sum(logdensity(x, family, shapes, loc, scale) for x in exact)
        value += sum(oracle(x, family, shapes, loc, scale)[1] for x in left)
        value += sum(oracle(x, family, shapes, loc, scale)[2] for x in right)
        value += sum(interval(a, b, family, shapes, loc, scale) for a, b in intervals)
        expected = float(value)
        assert not math.isnan(expected), (family, shapes, label, value)
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{96+family} {shapes[0] if shapes else 0.:.17g} {loc:.17g} {scale:.17g} '
                      f'{len(exact)} {len(left)} {len(right)} {len(intervals)}')
        inputs.extend(f'{s:.17g}' for s in shapes[1:])
        for group in (exact, left, right, *intervals):
            if group:
                inputs.append(' '.join(f'{x:.17g}' for x in group))

    if args.reuse:
        cases = json.loads((cache/'cases.json').read_text())
        inputs = (cache/'input.txt').read_text().splitlines()
    else:
        with mp.workdps(110):
            for family in range(3):
                for _ in range(args.cases):
                    shapes = (() if family == 0 else (rng.uniform(.001, .99),) if family == 1
                              else (10**rng.uniform(-1, 2),))
                    scale = 10**rng.uniform(-200, 200)
                    loc = scale*rng.uniform(-.3, .3)
                    point = lambda z: loc+scale*z
                    a = point(rng.uniform(.2, .6))
                    b = math.nextafter(a, math.inf)
                    add(family, shapes, loc, scale, exact=[point(.4)], left=[point(.1)],
                        right=[point(.8)], intervals=[(a, b), (point(.2), point(.7))])
                print(FAMILIES[family], 'mixed references ready', flush=True)
            for x in (1e-150, 1e-10, .01, .1, .81, .82, .83, 1., 4., 20., 1e100, 1e150):
                add(0, (), exact=[x], left=[x], right=[x], label='Kolmogorov tails')
                add(0, (), intervals=[(x, math.nextafter(x, math.inf))], label='Kolmogorov adjacent')
            for c in (1e-300, .001, .4, .99):
                for x in (1e-300, .001, 1., 3., math.nextafter(2*math.pi, -math.inf)):
                    add(1, (c,), exact=[x], left=[x], right=[x], label='angular tails')
                    add(1, (c,), intervals=[(x, math.nextafter(x, math.inf))], label='angular adjacent')
            for rho in (.1, 1., 100.):
                cutoff = 2*max(rho, math.sqrt(rho))
                for x in (1e-300, .01*rho, rho, math.nextafter(cutoff, 0.), cutoff, 1e100, 1e300):
                    add(2, (rho,), exact=[x], left=[x], right=[x], label='resonance tails')
                    add(2, (rho,), intervals=[(x, math.nextafter(x, math.inf))], label='resonance adjacent')
            for family, shapes in enumerate(((), (.4,), (1.,))):
                for scale in (1e-300, 1e300):
                    x = .4*scale
                    add(family, shapes, scale=scale, exact=[x], left=[x], right=[x],
                        intervals=[(x, math.nextafter(x, math.inf))], label='extreme scale')
                x = math.nextafter(0., 1.)
                add(family, shapes, exact=[x], left=[x], right=[x],
                    intervals=[(x, 2*x)], label='subnormal')
            for family, shapes, z in ((0, (), .8), (2, (1.,), 1.), (2, (100.,), 100.4)):
                for scale in (1e-200, 1., 1e200):
                    for ratio in (1e8, 1e10, 1e12):
                        for sign in (-1, 1):
                            loc = sign*ratio*scale
                            a = loc+z*scale
                            b = math.nextafter(a, math.inf)
                            add(family, shapes, loc, scale, intervals=[(a, b)],
                                label='large-location interval')
            add(0, (), left=[18.], label='tiny complement')
            add(0, (), right=[.05], label='tiny complement')
            add(1, (.4,), right=[1e-300], label='tiny complement')
            add(2, (1.,), right=[1e-300], label='tiny complement')
            add(2, (1.,), left=[1e100], label='tiny complement')
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
            if not math.isfinite(expected):
                error = 0. if value == expected else math.inf
            elif label == 'tiny complement':
                error = abs(value/expected-1)
            else:
                error = abs(value-expected)/max(1., abs(expected))
            if math.isnan(value) or error > 3e-11:
                failures.append((family, label, value, expected, error))
            if label == 'tiny complement':
                worst_tiny = max(worst_tiny, error)
            else:
                worst = max(worst, error)
        assert not failures, (mode, failures)
        print(f'{mode}: {len(cases)} comparisons passed; max scaled log error={worst}; '
              f'max relative tiny-log error={worst_tiny}', flush=True)
        scipy_stress(out/'driver.exe')
    sample = stats.CensoredData(uncensored=[.2, .4, .6], left=[.1], right=[.8],
                                interval=[[.3, .35], [.5, .55]])
    for name, shapes in zip(FAMILIES, ((), (.4,), (1.,))):
        fit = optimize.minimize_scalar(
            lambda v: getattr(stats, name)._penalized_nnlf((*shapes, 0., v), sample),
            bounds=(.3, 3.), method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: SciPy scale={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
