# SPDX-License-Identifier: MIT
"""Independent mpmath density/quadrature checks for censored Bessel families."""
import argparse
import json
import math
from pathlib import Path
import random
import subprocess
import mpmath as mp
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
NAMES = ('geninvgauss', 'norminvgauss', 'genhyperbolic')


def model(family, shapes):
    shapes = list(map(mp.mpf, shapes))
    if family == 0:
        order, a = shapes
        norm = -mp.log(2*mp.besselk(order, a))
        def logpdf(z):
            return norm+(order-1)*mp.log(z)-a*(z+1/z)/2
        return logpdf, order, a, mp.mpf(0)
    order, a, b = ([-mp.mpf('.5'), *shapes] if family == 1 else shapes)
    delta = mp.sqrt((a-b)*(a+b))
    q = order-mp.mpf('.5')
    norm = order*mp.log(delta)-mp.log(2*mp.pi)/2-q*mp.log(a)-mp.log(mp.besselk(order, delta))
    def logpdf(z):
        r = mp.sqrt(1+z*z)
        return norm+mp.log(mp.besselk(q, a*r))+b*z+q*mp.log(r)
    return logpdf, order, a, b


def tail(z, family, shapes, upper):
    extra = int(abs(mp.log10(abs(z))))+80 if z else 80
    with mp.workdps(max(mp.mp.dps, extra)):
        return integrate_tail(z, family, shapes, upper)


def integrate_tail(z, family, shapes, upper):
    logpdf, order, a, b = model(family, shapes)
    if family == 0 and z <= 0:
        return mp.mpf(0) if upper else -mp.inf
    direction = 1 if upper else -1
    if family != 0 and abs(z) < 1000:
        # Independent normal variance-mean mixture, not the Fortran y integral.
        delta = mp.sqrt((a-b)*(a+b))
        norm = order*mp.log(delta)-mp.log(2*mp.besselk(order, delta))
        def logintegrand(u):
            w = mp.exp(u)
            logmix = norm+order*u-(delta*delta*w+1/w)/2
            arg = (z-b*w)/mp.sqrt(2*w)
            return logmix+mp.log(mp.erfc(arg if upper else -arg)/2)
        center = max((mp.mpf(j)/4 for j in range(-56, 57)), key=logintegrand)
        peak = logintegrand(center)
        breaks = sorted(set([mp.mpf(t) for t in (-14, -7, -3, 0, 3, 7, 14)] +
                            [center+t for t in (-1, -.5, 0, .5, 1)]))
        return peak+mp.log(mp.quad(lambda u: mp.exp(logintegrand(u)-peak),
                                   breaks, maxdegree=5))
    y = mp.log(z) if family == 0 else mp.asinh(z)
    def h(u):
        if family == 0:
            return logpdf(mp.exp(u))+u
        return logpdf(mp.sinh(u))+mp.log(mp.cosh(u))
    # Adaptive tanh-sinh quadrature, independently evaluated at high precision.
    # Scale remote tails to avoid integrating a very thin boundary layer.
    slope = max(1, (a*mp.cosh(y)+abs(order)) if family == 0
                else (a-direction*b)*mp.cosh(y)+abs(order)+2)
    if abs(y) < 3:
        lo, hi = (-14, y) if not upper else (y, 14)
        peak = max(h(mp.mpf(t)) for t in range(math.floor(float(lo)), math.ceil(float(hi))+1))
        return peak+mp.log(mp.quad(lambda u: mp.exp(h(u)-peak),
                                   [lo, (lo+hi)/2, hi], maxdegree=5))
    base = h(y)
    integral = mp.quad(lambda t: mp.exp(h(y+direction*t/slope)-base),
                       [0, 1, 4, 16, 64, 256], maxdegree=5)
    return base-mp.log(slope)+mp.log(integral)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse', action='store_true')
    args = parser.parse_args()
    cache = ROOT/'build'/'bessel-reference'
    cache.mkdir(exist_ok=True)
    cases, inputs = [], []
    def add(family, shapes, x, part, loc=0., scale=1., end=None, label='reference'):
        z = (mp.mpf(x)-mp.mpf(loc))/mp.mpf(scale)
        logpdf, _, _, _ = model(family, shapes)
        if part == 0:
            expected = logpdf(z)-mp.log(scale)
        elif part in (1, 2):
            # Integrate the smaller tail and use a stable complement.
            small_upper = z > 0 if family != 0 else z > 1
            small = tail(z, family, shapes, small_upper)
            if small > -mp.log(2):
                small_upper = not small_upper
                small = tail(z, family, shapes, small_upper)
            expected = small if (part == 2) == small_upper else mp.log1p(-mp.exp(small))
        else:
            width = (mp.mpf(end)-mp.mpf(x))/mp.mpf(scale)
            mid = z+width/2
            base = logpdf(mid)
            if abs(width) < mp.mpf('1e-12'):
                expected = base+mp.log(width)
            elif abs(width) < mp.mpf('.001'):
                # Independent composite Simpson rule; O(width**4) relative error.
                weights = (1, 4, 2, 4, 1)
                total = sum(w*mp.exp(logpdf(z+width*j/4)-base)
                            for j, w in enumerate(weights))/12
                expected = base+mp.log(width)+mp.log(total)
            else:
                expected = base+mp.log(width)+mp.log(mp.quad(
                    lambda t: mp.exp(logpdf(z+width*t)-base), [0, .5, 1], maxdegree=5))
        counts = [0]*4
        counts[part] = 1
        inputs.append(f'{99+family} {shapes[0]:.17g} {loc:.17g} {scale:.17g} '+
                      ' '.join(map(str, counts)))
        inputs.extend(f'{v:.17g}' for v in shapes[1:])
        inputs.append(f'{x:.17g}'+(' '+f'{end:.17g}' if end is not None else ''))
        cases.append((NAMES[family], label, float(expected)))
    if args.reuse:
        cases = json.loads((cache/'cases.json').read_text())
        inputs = (cache/'input.txt').read_text().splitlines()
    else:
        with mp.workdps(65):
            for family, shapes in enumerate(((.7, 1.2), (1.5, .3), (1., 1.5, .3))):
                for x in ((.03, .7, 30.) if family == 0 else (-30., .7, 30.)):
                    for part in (range(3) if family != 1 or x in (-30., .7) else (0,)):
                        add(family, shapes, x, part, label='ordinary and remote tails')
                    add(family, shapes, x, 3, end=math.nextafter(x, math.inf), label='adjacent')
                for ratio in (1e8, 1e12):
                    for sign in (-1, 1):
                        for scale in (1e-200, 1., 1e200):
                            loc = sign*ratio*scale
                            x = loc+.7*scale
                            add(family, shapes, x, 3, loc, scale, math.nextafter(x, math.inf),
                                'large-location adjacent')
                for x in ((1e-100, 1e-9, 1e9, 1e100) if family == 0 else (-1e100, -1e9, 1e9, 1e100)):
                    for part in (0, 1 if x < 1 else 2):
                        add(family, shapes, x, part, label='extreme log tail')
                x = 400. if family == 0 else 200.
                add(family, shapes, x, 1, label='tiny complement')
                for x in ((1e-9, 1e9) if family == 0 else (-1e9, 1e9)):
                    add(family, shapes, x, 3, end=math.nextafter(x, math.inf),
                        label='remote adjacent')
                print(NAMES[family], 'multiprecision references ready', flush=True)
            add(0, (1.3957936504673851, 6.9007004756020205), .00644270263257127,
                1, label='concentrated lower-tail regression')
            # Shape boundary and half-integer Bessel orders.
            for family, shapes in ((0, (-5., .5)), (0, (5., 10.)),
                                    (1, (.5, -.4)), (1, (10., 8.)),
                                    (2, (-5., .5, .4)), (2, (5., 10., -8.)),
                                    (2, (0., 1.5, .3)), (2, (1., 1.5, -.3))):
                x = .7 if family == 0 else -.7
                add(family, shapes, x, 0, label='shape boundary')
                add(family, shapes, x, 3, end=math.nextafter(x, math.inf), label='shape boundary adjacent')
        (cache/'cases.json').write_text(json.dumps(cases), encoding='ascii')
        (cache/'input.txt').write_text('\n'.join(inputs)+'\n', encoding='ascii')
    for name, flags in (('checked', ['-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow']),
                         ('optimized', ['-O3', '-fcheck=bounds'])):
        out = cache/name
        out.mkdir(exist_ok=True)
        cmd = ['gfortran', *flags, '-J'+str(out), '-I'+str(out),
               '-Ibuild/gfortran_2806A042DAA8BBF7', 'src/special/log_gamma.f90',
               'src/empirical_distribution.f90', 'src/censored.f90',
               'validation/censored_reference_driver.f90',
               'build/gfortran_3C159EB802C22FB9/scifort/libscifort.a', '-o', str(out/'driver.exe')]
        print(subprocess.list2cmdline(cmd), flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)
        result = subprocess.run([str(out/'driver.exe')], input='\n'.join(inputs)+'\n',
                                text=True, capture_output=True, check=True)
        actual = list(map(float, result.stdout.split()))
        assert len(actual) == len(cases)
        worst = 0.
        for v, case in zip(actual, cases):
            family, label, expected = case
            error = abs(v/expected-1) if label == 'tiny complement' else abs(v-expected)/max(1., abs(expected))
            assert math.isfinite(v) and error < 3e-11, (name, case, v, error)
            worst = max(worst, error)
        print(name, len(cases), 'MP comparisons passed; max error', worst, flush=True)
        # Ordinary values only: SciPy's generic quadrature can lose deep tails.
        rng = random.Random(20261118)
        stress, expected = [], []
        for family in range(3):
            for _ in range(20):
                a = rng.uniform(.5, 10)
                b = rng.uniform(-.8, .8)*a
                p = rng.uniform(-5, 5)
                sh = (p, a) if family == 0 else (a, b) if family == 1 else (p, a, b)
                x = rng.uniform(.5, 2) if family == 0 else rng.uniform(-1, 1)
                stress.append(f'{99+family} {sh[0]} 0 1 1 0 0 0')
                stress.extend(map(str, sh[1:]))
                stress.append(str(x))
                expected.append(getattr(stats, NAMES[family]).logpdf(x, *sh))
        result = subprocess.run([str(out/'driver.exe')], input='\n'.join(stress)+'\n',
                                text=True, capture_output=True, check=True)
        for v, expected_value in zip(map(float, result.stdout.split()), expected):
            assert abs(v-expected_value) < 3e-11*max(1., abs(expected_value))
        print(name, len(expected), 'seeded SciPy density comparisons passed', flush=True)


if __name__ == '__main__':
    main()
