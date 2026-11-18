#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""High-precision noncentral censoring checks against Bessel/hypergeometric densities."""
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
FAMILIES = ('ncx2', 'rice', 'ncf')


def logdensity(x, family, shapes, loc, scale):
    x, loc, scale = map(mp.mpf, (x, loc, scale))
    z = (x-loc)/scale
    s = tuple(map(mp.mpf, shapes))
    if z <= 0:
        return -mp.inf
    if family == 1:
        b, = s
        return mp.log(z)-(z*z+b*b)/2+mp.log(mp.besseli(0, b*z))-mp.log(scale)
    if family == 0:
        df, nc = s
        a = df/2
        if nc == 0:
            return (a-1)*mp.log(z)-z/2-a*mp.log(2)-mp.loggamma(a)-mp.log(scale)
        return (-mp.log(2)-(nc+z)/2+(a-1)/2*mp.log(z/nc)
                +mp.log(mp.besseli(a-1, mp.sqrt(nc*z)))-mp.log(scale))
    df, dfd, nc = s
    a, b, lam = df/2, dfd/2, nc/2
    t = z*df/dfd
    u = t/(1+t)
    return (-lam+(a-1)*mp.log(t)-(a+b)*mp.log1p(t)-mp.log(mp.beta(a, b))
            +mp.log(mp.hyp1f1(a+b, a, lam*u))+mp.log(df/dfd)-mp.log(scale))


def oracle(x, family, shapes, loc, scale):
    z = (mp.mpf(x)-mp.mpf(loc))/mp.mpf(scale)
    s = tuple(map(mp.mpf, shapes))
    if family == 1:
        a, lam, t = mp.mpf(1), s[0]**2/2, z*z/2
    elif family == 0:
        a, lam, t = s[0]/2, s[1]/2, z/2
    else:
        a, b, lam = s[0]/2, s[1]/2, s[2]/2
        t = z*s[0]/s[1]
    if z <= 0:
        return -mp.inf, -mp.inf, mp.mpf(0)
    sums = [mp.mpf(0), mp.mpf(0)]
    weight = mp.exp(-lam)
    for j in range(3000):
        aj = a+j
        if family == 2:
            if t <= 1:
                lower = mp.betainc(aj, b, 0, t/(1+t), regularized=True)
                upper = 1-lower if lower < mp.mpf('.5') else mp.betainc(
                    b, aj, 0, 1/(1+t), regularized=True)
            else:
                upper = mp.betainc(b, aj, 0, 1/(1+t), regularized=True)
                lower = 1-upper if upper < mp.mpf('.5') else mp.betainc(
                    aj, b, 0, t/(1+t), regularized=True)
        else:
            lower = mp.gammainc(aj, 0, t, regularized=True)
            upper = (1-lower if lower < mp.mpf('.5') else
                     mp.gammainc(aj, t, mp.inf, regularized=True))
        terms = (weight*lower, weight*upper)
        sums = [v+w for v, w in zip(sums, terms)]
        if lam == 0 or (j > lam+50 and all(w <= v*mp.mpf('1e-85')
                                                         for v, w in zip(sums, terms))):
            break
        weight *= lam/(j+1)
    else:
        raise RuntimeError('Reference mixture did not converge')
    p, q = sums
    lp = mp.log(p) if p < mp.mpf('.5') else mp.log1p(-q)
    lq = mp.log(q) if q < mp.mpf('.5') else mp.log1p(-p)
    return logdensity(x, family, shapes, loc, scale), lp, lq


def interval(a, b, family, shapes, loc, scale):
    a, b = mp.mpf(a), mp.mpf(b)
    mid = (a+b)/2
    ref = logdensity(mid, family, shapes, loc, scale)
    def relative(t):
        return mp.exp(logdensity(a+(b-a)*t, family, shapes, loc, scale)-ref)
    return ref+mp.log(b-a)+mp.log(mp.quad(relative, [0, 1], method='gauss-legendre'))


def scipy_stress(driver):
    """Use multiprecision to adjudicate discrepancies with SciPy 1.15.3."""
    rng = random.Random(2026111803)
    inputs, refs = [], []
    for family, name in enumerate(FAMILIES):
        for _ in range(200):
            shapes = ((10**rng.uniform(-1, 3.3), rng.uniform(0, 100)) if family == 0 else
                      (rng.uniform(0, 10),) if family == 1 else
                      (10**rng.uniform(-1, 3.3), 10**rng.uniform(-1, 3.3), rng.uniform(0, 100)))
            x = 10**rng.uniform(-100, 3 if family < 2 else 100)
            dist = getattr(stats, name)(*shapes)
            for mode, method in enumerate(('logpdf', 'logcdf', 'logsf')):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    expected = float(getattr(dist, method)(x))
                if not math.isfinite(expected) or (mode and abs(expected) < 1e-7):
                    continue
                counts = [0, 0, 0]
                counts[mode] = 1
                inputs.append(f'{93+family} {shapes[0]:.17g} 0 1 '
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
    cache = ROOT/'build'/'noncentral_reference'
    cache.mkdir(parents=True, exist_ok=True)
    cases, inputs = [], []
    rng = random.Random(2026111802)

    def add(family, shapes, loc=0., scale=1., exact=(), left=(), right=(),
            intervals=(), label='mixed'):
        value = sum(logdensity(x, family, shapes, loc, scale) for x in exact)
        value += sum(oracle(x, family, shapes, loc, scale)[1] for x in left)
        value += sum(oracle(x, family, shapes, loc, scale)[2] for x in right)
        value += sum(interval(a, b, family, shapes, loc, scale) for a, b in intervals)
        expected = float(value)
        assert math.isfinite(expected), (family, shapes, label, value)
        cases.append((FAMILIES[family], label, expected))
        inputs.append(f'{93+family} {shapes[0]:.17g} {loc:.17g} {scale:.17g} '
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
                    shapes = ((10**rng.uniform(-1, 2), rng.uniform(0, 100)) if family == 0 else
                              (rng.uniform(0, 10),) if family == 1 else
                              (10**rng.uniform(-1, 2), 10**rng.uniform(-1, 2), rng.uniform(0, 100)))
                    scale = 10**rng.uniform(-200, 200)
                    loc = scale*rng.uniform(-.3, .3)
                    point = lambda z: loc+scale*z
                    a = point(rng.uniform(.2, .6))
                    b = math.nextafter(a, math.inf)
                    add(family, shapes, loc, scale, exact=[point(.4)], left=[point(.1)],
                        right=[point(.8)], intervals=[(a, b), (point(.2), point(.7))])
                print(FAMILIES[family], 'mixed references ready', flush=True)
            for family, shapes_set in enumerate((
                    ((.1, 0.), (.1, 100.), (2., 1.), (2000., 100.)),
                    ((0.,), (.01,), (10.,)),
                    ((.1, .1, 0.), (.1, 2000., 100.), (2000., .1, 100.), (2., 3., 1.)))):
                for shapes in shapes_set:
                    for x in (1e-300, .01, 1., 40., 100.):
                        add(family, shapes, exact=[x], left=[x], right=[x], label='tails')
                        if x != 100.:
                            add(family, shapes, intervals=[(x, math.nextafter(x, math.inf))],
                                label='adjacent')
                    print(FAMILIES[family], shapes, 'stress ready', flush=True)
            for family, shapes in enumerate(((2., 1.), (1.,), (2., 3., 1.))):
                for scale in (1e-300, 1e300):
                    x = .4*scale
                    add(family, shapes, scale=scale, exact=[x], left=[x], right=[x],
                        intervals=[(x, math.nextafter(x, math.inf))], label='extreme scale')
                x = math.nextafter(0., 1.)
                add(family, shapes, exact=[x], left=[x], right=[x],
                    intervals=[(x, 2*x)], label='subnormal')
            add(0, (2., 0.), left=[1400.], label='tiny complement')
            add(0, (2., 1.), right=[1e-300], label='tiny complement')
            add(1, (1.,), left=[35.], label='tiny complement')
            add(1, (1.,), right=[1e-150], label='tiny complement')
            add(2, (2., 3., 1.), left=[1e200], label='tiny complement')
            add(2, (2., 3., 1.), right=[1e-300], label='tiny complement')
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
        scipy_stress(out/'driver.exe')
    sample = stats.CensoredData(uncensored=[.2, .4, .6], left=[.1], right=[.8],
                                interval=[[.3, .35], [.5, .55]])
    for name, shapes in zip(FAMILIES, ((2., 1.), (1.,), (2., 3., 1.))):
        fit = optimize.minimize_scalar(
            lambda v: getattr(stats, name)._penalized_nnlf((*shapes, 0., v), sample),
            bounds=(.3, 3.), method='bounded', options={'xatol': 1e-13})
        assert fit.success
        print(f'{name}: SciPy scale={fit.x:.17g}; NNLF={fit.fun:.17g}')


if __name__ == '__main__':
    main()
