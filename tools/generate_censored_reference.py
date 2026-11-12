#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Print independent censored-likelihood references and SciPy fit comparisons."""
import json
import math

import mpmath as mp
import scipy
from scipy import stats

mp.mp.dps = 100


def logpdf(x, family, shape, loc, scale):
    x, shape, loc, scale = map(mp.mpf, (x, shape, loc, scale))
    z = (x - loc) / scale
    if family == 'normal':
        return -z*z/2 - mp.log(scale) - mp.log(2*mp.pi)/2
    if family == 'lognormal':
        if z <= 0:
            return mp.ninf
        return -mp.log(z)**2/(2*shape**2)-mp.log(shape)-mp.log(x-loc)-mp.log(2*mp.pi)/2
    if family == 'gamma':
        return (shape-1)*mp.log(z)-z-mp.loggamma(shape)-mp.log(scale)
    if family == 'exponential':
        return -z - mp.log(scale)
    return mp.log(shape/scale) + (shape-1)*mp.log(z) - z**shape


def logcdf(x, family, shape, loc, scale):
    x, shape, loc, scale = map(mp.mpf, (x, shape, loc, scale))
    z = (x-loc)/scale
    if family == 'normal':
        return mp.log(mp.erfc(-z/mp.sqrt(2))/2)
    if z <= 0:
        return mp.ninf
    if family == 'gamma':
        return mp.log(mp.gammainc(shape,0,z,regularized=True))
    if family == 'lognormal':
        return mp.log(mp.erfc(-mp.log(z)/shape/mp.sqrt(2))/2)
    return mp.log(-mp.expm1(-z**shape))


def logsf(x, family, shape, loc, scale):
    x, shape, loc, scale = map(mp.mpf, (x, shape, loc, scale))
    z = (x-loc)/scale
    if family == 'normal':
        return mp.log(mp.erfc(z/mp.sqrt(2))/2)
    if family == 'gamma':
        return mp.log(mp.gammainc(shape,max(z,mp.mpf(0)),mp.inf,regularized=True))
    if family == 'lognormal':
        return mp.log(mp.erfc(mp.log(z)/shape/mp.sqrt(2))/2) if z > 0 else mp.mpf(0)
    return -max(z, mp.mpf(0))**shape


def loginterval(a, b, family, shape=1., loc=0., scale=1.):
    if logsf(a, family, shape, loc, scale) < logcdf(b, family, shape, loc, scale):
        large = logsf(a, family, shape, loc, scale)
        small = logsf(b, family, shape, loc, scale)
    else:
        large = logcdf(b, family, shape, loc, scale)
        small = logcdf(a, family, shape, loc, scale)
    return large + mp.log(-mp.expm1(small-large))


def main():
    data = stats.CensoredData(uncensored=[.5, 1.2, 2.5], left=[.3], right=[3.5],
                              interval=[[1.1, 1.4], [4., 4.4]])
    result = {'scipy_version': scipy.__version__, 'mpmath_version': mp.__version__}
    for family, distribution, shape, loc, scale in [
        ('normal', stats.norm, 1., .4, 1.3),
        ('exponential', stats.expon, 1., 0., 2.3),
        ('weibull', stats.weibull_min, 1.7, 0., 2.3),
        ('lognormal', stats.lognorm, .8, 0., 2.3),
        ('gamma', stats.gamma, .6, 0., 2.3),
    ]:
        value = sum(logpdf(x, family, shape, loc, scale) for x in [.5, 1.2, 2.5])
        value += logcdf(.3, family, shape, loc, scale)
        value += logsf(3.5, family, shape, loc, scale)
        value += sum(loginterval(a, b, family, shape, loc, scale)
                     for a, b in [(1.1, 1.4), (4., 4.4)])
        result[family+'_mixed'] = float(value)
        kwargs = {} if family == 'normal' else {'floc': 0.}
        result[family+'_fit'] = list(distribution.fit(data, **kwargs))
    intervals = [
        ('lognormal_adjacent', 1e100, math.nextafter(1e100, math.inf), 'lognormal', .8, 1e100),
        ('lognormal_tail', math.exp(40), math.exp(40.01), 'lognormal', 1., 1.),
        ('gamma_wide_left', .1, .2, 'gamma', .4, 1.),
        ('gamma_wide_tail', 100., 101., 'gamma', 2.5, 1.),
        ('gamma_boundary', 0., 1e-300, 'gamma', .4, 1e300),
        ('gamma_adjacent', 1., math.nextafter(1., math.inf), 'gamma', .4, 1.),
        ('gamma_tail', 1000., 1000.+1e-10, 'gamma', 3.5, 1.),
        ('normal_narrow', 0., 1e-20, 'normal', 1., 1.),
        ('normal_tail', 40., 40.01, 'normal', 1., 1.),
        ('exponential_tail', 1000., 1000.+1e-10, 'exponential', 1., 1.),
        ('exponential_underflow', 0., 1e-300, 'exponential', 1., 1e300),
        ('weibull_tail', 100., 100.+1e-8, 'weibull', 2., 1.),
        ('weibull_narrow', 1., math.nextafter(1., math.inf), 'weibull', .01, 1.),
    ]
    # Resolving a 1e-600 interval probability needs more than 100 digits when
    # subtracting high-precision logarithms near zero.
    with mp.workdps(700):
        for name, a, b, family, shape, scale in intervals:
            result[name] = float(loginterval(a, b, family, shape, scale=scale))
    with mp.workdps(400):
        result['gamma_small_sf'] = float(logsf(1e-300, 'gamma', 1e-20, 0., 1e300))
        result['gamma_small_pdf'] = float(logpdf(1e-300, 'gamma', .4, 0., 1e300))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
