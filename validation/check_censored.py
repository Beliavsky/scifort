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
    for number,family in [(8,'gumbel_r'),(9,'gumbel_l')]:
        for _ in range(args.cases):
            loc,scale=rng.uniform(-3.,3.),10**rng.uniform(-4.,4.)
            exact=[loc+scale*rng.uniform(-5.,8.) for _ in range(3)]
            left=[loc+scale*rng.uniform(-7.,8.)]
            right=[loc+scale*rng.uniform(-7.,8.)]
            intervals=[]
            for _ in range(2):
                a=loc+scale*rng.uniform(-30.,1000.)
                b=a+scale*10**rng.uniform(-14.,3.)
                if b == a: b=math.nextafter(a,math.inf)
                intervals.append((a,b))
            # Reflect the sampling range as well as the family for comparable tails.
            if family == 'gumbel_l':
                exact=[2*loc-x for x in exact]
                left=[2*loc-x for x in left]
                right=[2*loc-x for x in right]
                intervals=[(2*loc-b,2*loc-a) for a,b in intervals]
                intervals=[(a,math.nextafter(a,math.inf) if b == a else b)
                           for a,b in intervals]
            value=sum(reference.logpdf(x,family,1.,loc,scale) for x in exact)
            value+=sum(reference.logcdf(x,family,1.,loc,scale) for x in left)
            value+=sum(reference.logsf(x,family,1.,loc,scale) for x in right)
            value+=sum(reference.loginterval(a,b,family,1.,loc,scale) for a,b in intervals)
            cases.append((family,float(value)))
            inputs.append(f'{number} 1 {loc:.17g} {scale:.17g} 3 1 1 2')
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
    with mp.workdps(700):
        for a,b,loc,scale in [
            (0.,math.nextafter(0.,1.),0.,1.),
            (1e-300,2e-300,0.,1e300),
            (1000.,math.nextafter(1000.,math.inf),0.,1.),
            (-5.,math.nextafter(-5.,math.inf),0.,1.),
            (-700.,math.nextafter(-700.,math.inf),0.,1.),
            (-1e308,1e308,0.,1e308),
            (1e308,math.nextafter(1e308,math.inf),-1e308,1e308),
            (-1e-300,1e-300,0.,1e-300),
        ]:
            for number,family,lo,hi,center in [
                (8,'gumbel_r',a,b,loc),(9,'gumbel_l',-b,-a,-loc)]:
                cases.append((family+'_stress',float(reference.loginterval(
                    lo,hi,family,1.,center,scale))))
                stdin+=f'{number} 1 {center:.17g} {scale:.17g} 0 0 0 1\n{lo:.17g} {hi:.17g}\n'
    with mp.workdps(700):
        for _ in range(args.cases):
            loc, scale = rng.uniform(-3., 3.), 10**rng.uniform(-4., 4.)
            exact = [loc+scale*rng.uniform(-20.,20.) for _ in range(3)]
            left = [loc+scale*rng.uniform(-40.,40.)]
            right = [loc+scale*rng.uniform(-40.,40.)]
            intervals = []
            for _ in range(3):
                a = loc+scale*rng.uniform(-1000.,1000.)
                b = a+scale*10**rng.uniform(-14.,3.)
                if b == a: b = math.nextafter(a,math.inf)
                intervals.append((a,b))
            # Force a center-crossing interval as well as same-tail intervals.
            intervals.append((loc-scale*10**rng.uniform(-12.,2.),
                              loc+scale*10**rng.uniform(-12.,2.)))
            expected = sum(reference.logpdf(x,'laplace',1.,loc,scale) for x in exact)
            expected += sum(reference.logcdf(x,'laplace',1.,loc,scale) for x in left)
            expected += sum(reference.logsf(x,'laplace',1.,loc,scale) for x in right)
            expected += sum(reference.loginterval(a,b,'laplace',1.,loc,scale)
                            for a,b in intervals)
            cases.append(('laplace',float(expected)))
            stdin += f'10 1 {loc:.17g} {scale:.17g} 3 1 1 4\n'
            stdin += '\n'.join(' '.join(f'{x:.17g}' for x in group)
                               for group in [exact,left,right,*intervals])+'\n'
        for a,b,loc,scale in [
            (0.,math.nextafter(0.,1.),0.,1.),
            (-math.nextafter(0.,1.),math.nextafter(0.,1.),0.,1.),
            (1e-300,2e-300,0.,1e300),
            (-1e-300,1e-300,0.,1e300),
            (1000.,math.nextafter(1000.,math.inf),0.,1.),
            (-math.nextafter(1000.,math.inf),-1000.,0.,1.),
            (-1e308,1e308,0.,1e308),
            (1e308,math.nextafter(1e308,math.inf),-1e308,1e308),
            (-1e-300,1e-300,0.,1e-300),
        ]:
            expected = reference.loginterval(a,b,'laplace',1.,loc,scale)
            cases.append(('laplace_stress',float(expected)))
            stdin += f'10 1 {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
    with mp.workdps(700):
        for _ in range(args.cases):
            loc, scale = rng.uniform(-3., 3.), 10**rng.uniform(-4., 4.)
            exact = [loc+scale*rng.uniform(-20.,20.) for _ in range(3)]
            left = [loc+scale*rng.uniform(-40.,40.)]
            right = [loc+scale*rng.uniform(-40.,40.)]
            intervals = []
            for _ in range(3):
                a = loc+scale*rng.uniform(-1000.,1000.)
                b = a+scale*10**rng.uniform(-14.,3.)
                if b == a: b = math.nextafter(a,math.inf)
                intervals.append((a,b))
            # Force a center-crossing interval as well as same-tail intervals.
            intervals.append((loc-scale*10**rng.uniform(-12.,2.),
                              loc+scale*10**rng.uniform(-12.,2.)))
            expected = sum(reference.logpdf(x,'cauchy',1.,loc,scale) for x in exact)
            expected += sum(reference.logcdf(x,'cauchy',1.,loc,scale) for x in left)
            expected += sum(reference.logsf(x,'cauchy',1.,loc,scale) for x in right)
            expected += sum(reference.loginterval(a,b,'cauchy',1.,loc,scale)
                            for a,b in intervals)
            cases.append(('cauchy',float(expected)))
            stdin += f'11 1 {loc:.17g} {scale:.17g} 3 1 1 4\n'
            stdin += '\n'.join(' '.join(f'{x:.17g}' for x in group)
                               for group in [exact,left,right,*intervals])+'\n'
        for a,b,loc,scale in [
            (0.,math.nextafter(0.,1.),0.,1.),
            (-math.nextafter(0.,1.),math.nextafter(0.,1.),0.,1.),
            (1e-300,2e-300,0.,1e300),
            (-1e-300,1e-300,0.,1e300),
            (1e300,math.nextafter(1e300,math.inf),0.,1.),
            (-math.nextafter(1e300,math.inf),-1e300,0.,1.),
            (1e308,math.nextafter(1e308,math.inf),0.,math.nextafter(0.,1.)),
            (-1e308,1e308,0.,1e308),
            (1e308,math.nextafter(1e308,math.inf),-1e308,1e308),
            (-1e-300,1e-300,0.,1e-300),
            (-1e300,1e300,0.,1.),
            (-1.,1.,0.,math.nextafter(0.,1.)),
            (0.,1e300,0.,1e-300),
            (0.,math.nextafter(0.,1.),0.,1e308),
            (1e-300,2e-300,0.,math.nextafter(0.,1.)),
            (1.,math.nextafter(1.,math.inf),1.,1e-300),
            (0.,1.,0.,1.),
            (-1.,1.,0.,1.),
            (-1e-300,1e300,0.,1.),
        ]:
            expected = reference.loginterval(a,b,'cauchy',1.,loc,scale)
            cases.append(('cauchy_stress',float(expected)))
            stdin += f'11 1 {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
    with mp.workdps(700):
        for number,family in [(12,'halfcauchy'),(13,'halflogistic')]:
            for _ in range(args.cases):
                loc, scale = rng.uniform(-3.,3.), 10**rng.uniform(-4.,4.)
                exact = [loc+scale*rng.uniform(.01,20.) for _ in range(3)]
                left = [loc+scale*rng.uniform(.01,40.)]
                right = [loc+scale*rng.uniform(-2.,40.)]
                intervals=[]
                for _ in range(2):
                    a=loc+scale*rng.uniform(.01,1000.)
                    b=a+scale*10**rng.uniform(-14.,2.)
                    if b == a: b=math.nextafter(a,math.inf)
                    intervals.append((a,b))
                b=loc+scale*10**rng.uniform(-12.,1.)
                if b <= loc: b=math.nextafter(loc,math.inf)
                intervals.append((loc-scale,b))
                value=sum(reference.logpdf(x,family,1.,loc,scale) for x in exact)
                value+=sum(reference.logcdf(x,family,1.,loc,scale) for x in left)
                value+=sum(reference.logsf(x,family,1.,loc,scale) for x in right)
                value+=sum(reference.loginterval(a,b,family,1.,loc,scale) for a,b in intervals)
                cases.append((family,float(value)))
                stdin+=f'{number} 1 {loc:.17g} {scale:.17g} 3 1 1 3\n'
                stdin+='\n'.join(' '.join(f'{x:.17g}' for x in group)
                                  for group in [exact,left,right,*intervals])+'\n'
            stress=[
                (0.,math.nextafter(0.,1.),0.,1.),
                (-math.nextafter(0.,1.),math.nextafter(0.,1.),0.,1.),
                (1e-300,2e-300,0.,1e300),
                (-1e-300,1e-300,0.,1e300),
                (1000.,math.nextafter(1000.,math.inf),0.,1.),
                (-1e308,1e308,0.,1e308),
                (1e308,math.nextafter(1e308,math.inf),-1e308,1e308),
                (-1e-300,1e-300,0.,1e-300),
                (0.,1.,0.,1.),
                (1.,math.nextafter(1.,math.inf),1.,1.),
                (-1.,1.,0.,1.),
                (0.,1000.,0.,1.),
            ]
            if family == 'halfcauchy':
                stress += [
                    (1e-300,1e300,0.,1.),
                    (1e300,math.nextafter(1e300,math.inf),0.,1.),
                    (1e308,math.nextafter(1e308,math.inf),0.,math.nextafter(0.,1.)),
                    (1e-300,2e-300,0.,math.nextafter(0.,1.)),
                    (0.,math.nextafter(0.,1.),0.,1e308),
                ]
            if family == 'halflogistic':
                stress += [(1e-300,700.,0.,1.)]
            for a,b,loc,scale in stress:
                value=reference.loginterval(a,b,family,1.,loc,scale)
                cases.append((family+'_stress',float(value)))
                stdin+=f'{number} 1 {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
    with mp.workdps(800):
        for number,family in [(14,'rayleigh'),(15,'weibull_max')]:
            for _ in range(args.cases):
                loc,scale=rng.uniform(-3.,3.),10**rng.uniform(-4.,4.)
                shape=10**rng.uniform(-1.,1.) if number == 15 else 1.
                exact=[loc+scale*rng.uniform(.01,4.) for _ in range(3)]
                left=[loc+scale*rng.uniform(.01,5.)]
                right=[loc+scale*rng.uniform(-1.,5.)]
                intervals=[]
                for _ in range(2):
                    a=loc+scale*rng.uniform(.01,20.)
                    b=a+scale*10**rng.uniform(-14.,.5)
                    if b == a: b=math.nextafter(a,math.inf)
                    intervals.append((a,b))
                b=loc+scale*10**rng.uniform(-12.,1.)
                if b <= loc: b=math.nextafter(loc,math.inf)
                intervals.append((loc-scale,b))
                if number == 15:
                    loc=-loc
                    exact=[-x for x in exact]
                    left,right=[-x for x in right],[-x for x in left]
                    intervals=[(-b,-a) for a,b in intervals]
                expected=sum(reference.logpdf(x,family,shape,loc,scale) for x in exact)
                expected+=sum(reference.logcdf(x,family,shape,loc,scale) for x in left)
                expected+=sum(reference.logsf(x,family,shape,loc,scale) for x in right)
                expected+=sum(reference.loginterval(a,b,family,shape,loc,scale) for a,b in intervals)
                cases.append((family,float(expected)))
                stdin+=f'{number} {shape:.17g} {loc:.17g} {scale:.17g} 3 1 1 3\n'
                stdin+='\n'.join(' '.join(f'{x:.17g}' for x in group)
                                  for group in [exact,left,right,*intervals])+'\n'
        for a,b,loc,scale in [
            (0.,math.nextafter(0.,1.),0.,1.),
            (-math.nextafter(0.,1.),math.nextafter(0.,1.),0.,1.),
            (1e-300,2e-300,0.,1e300),
            (-1e-300,1e-300,0.,1e300),
            (40.,math.nextafter(40.,math.inf),0.,1.),
            (-1e308,1e308,0.,1e308),
            (1e308,math.nextafter(1e308,math.inf),-1e308,1e308),
            (1e-300,40.,0.,1.),
            (1e-150,37.,0.,1.),
            (0.,1e308,0.,1.5e308),
        ]:
            expected=reference.loginterval(a,b,'rayleigh',1.,loc,scale)
            cases.append(('rayleigh_stress',float(expected)))
            stdin+=f'14 1 {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
        for a,b,shape,loc,scale in [
            (0.,math.nextafter(0.,1.),.3,0.,1.),
            (1e-300,2e-300,.3,0.,1e300),
            (1e-300,2e-300,2.,0.,1e300),
            (40.,math.nextafter(40.,math.inf),2.,0.,1.),
            (-1e308,1e308,2.,0.,1e308),
            (1e308,math.nextafter(1e308,math.inf),2.,-1e308,1e308),
            (1e-300,40.,2.,0.,1.),
            (0.,1e308,1.,0.,1.5e308),
        ]:
            a,b,loc=-b,-a,-loc
            expected=reference.loginterval(a,b,'weibull_max',shape,loc,scale)
            cases.append(('weibull_max_stress',float(expected)))
            stdin+=f'15 {shape:.17g} {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
    for number, family in [(16,'erlang'),(17,'gibrat'),(18,'halfnorm')]:
        for _ in range(args.cases):
            loc, scale = rng.uniform(-3., 3.), 10**rng.uniform(-4., 4.)
            c = 10**rng.uniform(-1., 1.) if number >= 3 else 1.
            if family in ('gibrat','halfnorm'): c = 1.
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
            stdin += f'{number} {c:.17g} {loc:.17g} {scale:.17g} 3 1 1 2\n'
            stdin += '\n'.join(' '.join(f'{x:.17g}' for x in group)
                                for group in [exact,left,right,*intervals])+'\n'
    with mp.workdps(800):
        for a,b,loc,scale in [
            (1.,math.nextafter(1.,math.inf),0.,1.),
            (math.exp(40.),math.nextafter(math.exp(40.),math.inf),0.,1.),
            (1e-300,2e-300,0.,1e300),
            (1e308,math.nextafter(1e308,math.inf),-1e308,1e308),
        ]:
            expected=reference.loginterval(a,b,'gibrat',1.,loc,scale)
            cases.append(('gibrat_stress',float(expected)))
            stdin+=f'17 1 {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
    with mp.workdps(800):
        for a,b,loc,scale in [
            (0.,math.nextafter(0.,1.),0.,1.),
            (-math.nextafter(0.,1.),math.nextafter(0.,1.),0.,1.),
            (1e-300,2e-300,0.,1e300),
            (-1e-300,1e-300,0.,1e300),
            (40.,math.nextafter(40.,math.inf),0.,1.),
            (-1e308,1e308,-1.1e308,1e308),
            (1e308,math.nextafter(1e308,math.inf),-1e308,1e308),
            (1e-20,10.,0.,1.),
            (0.,1e308,0.,1.5e308),
            (.01,math.nextafter(.01,math.inf),0.,1.),
            (1.,math.nextafter(1.,math.inf),0.,1.),
            (-1.,1.,0.,1.),
        ]:
            expected=reference.loginterval(a,b,'halfnorm',1.,loc,scale)
            cases.append(('halfnorm_stress',float(expected)))
            stdin+=f'18 1 {loc:.17g} {scale:.17g} 0 0 0 1\n{a:.17g} {b:.17g}\n'
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
        laplace_errors = []
        cauchy_errors = []
        half_errors = []
        hazard_errors = []
        erlang_errors = []
        gibrat_errors = []
        halfnorm_errors = []
        for value, (family, expected) in zip(values, cases):
            error = abs(value-expected)/max(1., abs(expected))
            if not math.isfinite(expected) or not math.isfinite(value) or error > 2e-11:
                raise AssertionError((mode, family, value, expected, error))
            worst = max(worst, (error, family))
            if family.startswith('halfnorm'):
                halfnorm_errors.append(error)
            if family.startswith('gibrat'):
                gibrat_errors.append(error)
            if family.startswith('erlang'):
                erlang_errors.append(error)
            if family.startswith(('rayleigh','weibull_max')):
                hazard_errors.append(error)
            if family.startswith(('halfcauchy','halflogistic')):
                half_errors.append(error)
            if family.startswith('cauchy'):
                cauchy_errors.append(error)
            if family.startswith('laplace'):
                laplace_errors.append(error)
        print(f'{mode}: {len(cases)} comparisons passed; '
              f'max scaled log-likelihood error={worst[0]:.3g} ({worst[1]})', flush=True)
        print(f'  Erlang: {len(erlang_errors)} comparisons; '
              f'max scaled error={max(erlang_errors):.3g}', flush=True)
        print(f'  Half-normal: {len(halfnorm_errors)} comparisons; '
              f'max scaled error={max(halfnorm_errors):.3g}', flush=True)
        print(f'  Gibrat: {len(gibrat_errors)} comparisons; '
              f'max scaled error={max(gibrat_errors):.3g}', flush=True)
        print(f'  Laplace: {len(laplace_errors)} comparisons; '
              f'max scaled error={max(laplace_errors):.3g}', flush=True)
        print(f'  Cauchy: {len(cauchy_errors)} comparisons; '
              f'max scaled error={max(cauchy_errors):.3g}', flush=True)
        print(f'  Half-Cauchy/logistic: {len(half_errors)} comparisons; '
              f'max scaled error={max(half_errors):.3g}', flush=True)
        print(f'  Rayleigh/Weibull max: {len(hazard_errors)} comparisons; '
              f'max scaled error={max(hazard_errors):.3g}', flush=True)


if __name__ == '__main__':
    main()
