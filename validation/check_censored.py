#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Randomized high-precision likelihood checks in checked and optimized builds."""
import argparse
import importlib.util
import math
from pathlib import Path
import random
import subprocess

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('censored_reference',
                                             ROOT / 'tools/generate_censored_reference.py')
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module-dir', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--seed', type=int, default=20261112)
    parser.add_argument('--cases', type=int, default=100, help='cases per family')
    args = parser.parse_args()
    rng = random.Random(args.seed)
    mp.mp.dps = 100
    cases, inputs = [], []
    for number, family in enumerate(['normal', 'exponential', 'weibull', 'lognormal', 'gamma'], 1):
        for _ in range(args.cases):
            loc, scale = rng.uniform(-3., 3.), 10**rng.uniform(-4., 4.)
            c = 10**rng.uniform(-1., 1.) if number >= 3 else 1.
            exact = [loc+scale*rng.uniform(.1, 4.) for _ in range(3)]
            left = [loc+scale*rng.uniform(.1, 5.)]
            right = [loc+scale*rng.uniform(.1, 5.)]
            intervals = []
            for _ in range(2):
                a = loc+scale*rng.uniform(-30. if number == 1 else .01, 30.)
                b = a+scale*10**rng.uniform(-12., .5)
                if b == a:
                    b = math.nextafter(a, math.inf)
                intervals.append((a,b))
            value = sum(reference.logpdf(x, family, c, loc, scale) for x in exact)
            value += sum(reference.logcdf(x, family, c, loc, scale) for x in left)
            value += sum(reference.logsf(x, family, c, loc, scale) for x in right)
            value += sum(reference.loginterval(a,b,family,c,loc,scale) for a,b in intervals)
            cases.append((family, float(value)))
            inputs.append(f'{number} {c:.17g} {loc:.17g} {scale:.17g} 3 1 1 2')
            inputs.extend(' '.join(f'{x:.17g}' for x in group)
                          for group in [exact,left,right,*intervals])
    for _ in range(args.cases):
        loc, scale = rng.uniform(-3.,3.), 10**rng.uniform(-4.,4.)
        shape, second = 10**rng.uniform(-1.,1.), 10**rng.uniform(-1.,1.)
        exact=[loc+scale*rng.uniform(.05,.95) for _ in range(3)]
        left=[loc+scale*rng.uniform(.01,.2)]
        right=[loc+scale*rng.uniform(.8,.99)]
        intervals=[]
        for _ in range(2):
            a=loc+scale*rng.uniform(.01,.95)
            b=min(a+scale*10**rng.uniform(-12.,-.3),loc+scale*1.1)
            if b == a: b=math.nextafter(a,math.inf)
            intervals.append((a,b))
        value=sum(reference.logpdf(x,'beta',shape,loc,scale,second) for x in exact)
        value+=sum(reference.logcdf(x,'beta',shape,loc,scale,second) for x in left)
        value+=sum(reference.logsf(x,'beta',shape,loc,scale,second) for x in right)
        value+=sum(reference.loginterval(a,b,'beta',shape,loc,scale,second) for a,b in intervals)
        cases.append(('beta',float(value)))
        inputs.append(f'6 {shape:.17g} {loc:.17g} {scale:.17g} 3 1 1 2')
        inputs.append(f'{second:.17g}')
        inputs.extend(' '.join(f'{x:.17g}' for x in group)
                      for group in [exact,left,right,*intervals])
    for _ in range(args.cases):
        loc, scale = rng.uniform(-3.,3.), 10**rng.uniform(-4.,4.)
        exact = [loc+scale*rng.uniform(-20.,20.) for _ in range(3)]
        left = [loc+scale*rng.uniform(-40.,40.)]
        right = [loc+scale*rng.uniform(-40.,40.)]
        intervals = []
        for _ in range(2):
            a=loc+scale*rng.uniform(-1000.,1000.)
            b=a+scale*10**rng.uniform(-14.,3.)
            if b == a: b=math.nextafter(a,math.inf)
            intervals.append((a,b))
        value=sum(reference.logpdf(x,'logistic',1.,loc,scale) for x in exact)
        value+=sum(reference.logcdf(x,'logistic',1.,loc,scale) for x in left)
        value+=sum(reference.logsf(x,'logistic',1.,loc,scale) for x in right)
        value+=sum(reference.loginterval(a,b,'logistic',1.,loc,scale) for a,b in intervals)
        cases.append(('logistic',float(value)))
        inputs.append(f'7 1 {loc:.17g} {scale:.17g} 3 1 1 2')
        inputs.extend(' '.join(f'{x:.17g}' for x in group)
                      for group in [exact,left,right,*intervals])
    stdin = '\n'.join(inputs)+'\n'
    a = math.nextafter(1e100, math.inf)
    b = math.nextafter(a, math.inf)
    cases.append(('weibull_large_shape', float(reference.loginterval(
        a, b, 'weibull', 1e16, 0., 1e100))))
    stdin += f'3 1e16 0 1e100 0 0 0 1\n{a:.17g} {b:.17g}\n'
    for a, b, shape, loc, scale in [
        (1e-300, math.nextafter(1e-300, math.inf), .8, 0., 1e300),
        (1e100, math.nextafter(1e100, math.inf), .8, 0., 1e100),
        (-1e300, 1e300, 2., -1.1e300, 1e-300),
    ]:
        cases.append(('lognormal_stress', float(reference.loginterval(
            a, b, 'lognormal', shape, loc, scale))))
        stdin += f'4 {shape:.17g} {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
    for a, b, shape, scale in [
        (0., 1e-300, .4, 1e300),
        (0., 1e-300, 1000., 1e300),
        (1e-300, math.nextafter(1e-300, math.inf), .4, 1e300),
        (1., math.nextafter(1., math.inf), .4, 1.),
        (1000., 1000.+1e-10, 3.5, 1.),
        (1000., math.nextafter(1000.,math.inf), 1000., 1.),
        (900., 1000., 1000., 1.),
        (1e6, 1e6+100., 1e6, 1.),
        (1e-300, 2e-300, 1e-20, 1e300),
    ]:
        cases.append(('gamma_stress', float(reference.loginterval(
            a, b, 'gamma', shape, 0., scale))))
        stdin += f'5 {shape:.17g} 0 {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
    with mp.workdps(700):
        for a,b,shape,second,loc,scale in [
            (0.,1e-300,.4,.6,0.,1.),
            (math.nextafter(1.,0.),1.,.4,.6,0.,1.),
            (.5,math.nextafter(.5,1.),.4,.6,0.,1.),
            (1e-300,2e-300,.4,.6,0.,1e300),
            (0.,1e-300,1e-20,.6,0.,1e300),
            (.49,.51,1000.,1000.,0.,1.),
            (.5,math.nextafter(.5,1.),1000.,1000.,0.,1.),
            (.99,.999,.4,3.5,0.,1.),
            (1e300,math.nextafter(1e300,math.inf),.4,.6,1e-300,1e300),
        ]:
            cases.append(('beta_stress',float(reference.loginterval(
                a,b,'beta',shape,loc,scale,second))))
            stdin+=f'6 {shape:.17g} {loc:.17g} {scale:.17g} 0 0 0 1\n{second:.17g}\n{a:.17g} {b:.17g}\n'
    with mp.workdps(700):
        for a,b,loc,scale in [
            (0.,math.nextafter(0.,1.),0.,1.),
            (1e-300,2e-300,0.,1e300),
            (1000.,math.nextafter(1000.,math.inf),0.,1.),
            (-1000.,math.nextafter(-1000.,math.inf),0.,1.),
            (-1e308,1e308,0.,1.),
            (1e308,math.nextafter(1e308,math.inf),-1e308,1e308),
            (-1e-300,1e-300,0.,1e-300),
        ]:
            cases.append(('logistic_stress',float(reference.loginterval(
                a,b,'logistic',1.,loc,scale))))
            stdin+=f'7 1 {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
    for mode, flags in [('checked', ['-O0','-g','-fcheck=all',
                                     '-ffpe-trap=invalid,zero,overflow','-fbacktrace']),
                        ('optimized', ['-O3','-fcheck=bounds'])]:
        directory = ROOT / 'build/censored_reference' / mode
        directory.mkdir(parents=True, exist_ok=True)
        exe = directory / 'driver.exe'
        source = [ROOT / path for path in ['src/special/log_gamma.f90',
                  'src/empirical_distribution.f90',
                  'src/censored.f90', 'validation/censored_reference_driver.f90']]
        subprocess.run([args.compiler,*flags,'-J'+str(directory),'-I'+str(directory),
                        '-I'+str(args.module_dir.resolve()), *map(str,source),
                        str(args.library.resolve()),'-o',str(exe)], cwd=ROOT, check=True)
        run = subprocess.run([str(exe)], input=stdin, text=True, capture_output=True, check=True)
        values = list(map(float, run.stdout.split()))
        assert len(values) == len(cases)
        worst = (0., '')
        for value, (family, expected) in zip(values, cases):
            error = abs(value-expected)/max(1., abs(expected))
            if not math.isfinite(expected) or not math.isfinite(value) or error > 2e-11:
                raise AssertionError((mode, family, value, expected, error))
            worst = max(worst, (error, family))
        print(f'{mode}: {len(cases)} comparisons passed; '
              f'max scaled log-likelihood error={worst[0]:.3g} ({worst[1]})', flush=True)


if __name__ == '__main__':
    main()
