#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Print independent censored-likelihood references and SciPy fit comparisons."""
import json
import math

import mpmath as mp
import scipy
from scipy import stats

mp.mp.dps = 100


def logpdf(x, family, shape, loc, scale, second_shape=1.):
    x, shape, loc, scale = map(mp.mpf, (x, shape, loc, scale))
    z = (x - loc) / scale
    if family == 'beta':
        second_shape=mp.mpf(second_shape)
        return ((shape-1)*mp.log(z)+(second_shape-1)*mp.log1p(-z)
                -mp.log(mp.beta(shape,second_shape))-mp.log(scale))
    if family in ('gumbel_r','gumbel_l'):
        if family == 'gumbel_l': z=-z
        return -z-mp.exp(-z)-mp.log(scale)
    if family == 'logistic':
        return -abs(z)-2*mp.log1p(mp.exp(-abs(z)))-mp.log(scale)
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


def logcdf(x, family, shape, loc, scale, second_shape=1.):
    x, shape, loc, scale = map(mp.mpf, (x, shape, loc, scale))
    z = (x-loc)/scale
    if family == 'gumbel_r':
        return -mp.exp(-z)
    if family == 'gumbel_l':
        return mp.log(-mp.expm1(-mp.exp(z)))
    if family == 'logistic':
        return min(z,0)-mp.log1p(mp.exp(-abs(z)))
    if family == 'normal':
        return mp.log(mp.erfc(-z/mp.sqrt(2))/2)
    if z <= 0:
        return mp.ninf
    if family == 'beta':
        return mp.log(mp.betainc(shape,mp.mpf(second_shape),0,min(z,mp.mpf(1)),regularized=True))
    if family == 'gamma':
        return mp.log(mp.gammainc(shape,0,z,regularized=True))
    if family == 'lognormal':
        return mp.log(mp.erfc(-mp.log(z)/shape/mp.sqrt(2))/2)
    return mp.log(-mp.expm1(-z**shape))


def logsf(x, family, shape, loc, scale, second_shape=1.):
    x, shape, loc, scale = map(mp.mpf, (x, shape, loc, scale))
    z = (x-loc)/scale
    if family == 'gumbel_r':
        return mp.log(-mp.expm1(-mp.exp(-z)))
    if family == 'gumbel_l':
        return -mp.exp(z)
    if family == 'logistic':
        return -max(z,0)-mp.log1p(mp.exp(-abs(z)))
    if family == 'normal':
        return mp.log(mp.erfc(z/mp.sqrt(2))/2)
    if family == 'beta':
        z=max(mp.mpf(0),min(z,mp.mpf(1)))
        return mp.log(mp.betainc(shape,mp.mpf(second_shape),z,1,regularized=True))
    if family == 'gamma':
        return mp.log(mp.gammainc(shape,max(z,mp.mpf(0)),mp.inf,regularized=True))
    if family == 'lognormal':
        return mp.log(mp.erfc(mp.log(z)/shape/mp.sqrt(2))/2) if z > 0 else mp.mpf(0)
    return -max(z, mp.mpf(0))**shape


def loginterval(a, b, family, shape=1., loc=0., scale=1., second_shape=1.):
    if logsf(a, family, shape, loc, scale, second_shape) < logcdf(b, family, shape, loc, scale, second_shape):
        large = logsf(a, family, shape, loc, scale, second_shape)
        small = logsf(b, family, shape, loc, scale, second_shape)
    else:
        large = logcdf(b, family, shape, loc, scale, second_shape)
        small = logcdf(a, family, shape, loc, scale, second_shape)
    return large + mp.log(-mp.expm1(small-large))


def main():
    data = stats.CensoredData(uncensored=[.5, 1.2, 2.5], left=[.3], right=[3.5],
                              interval=[[1.1, 1.4], [4., 4.4]])
    result = {'scipy_version': scipy.__version__, 'mpmath_version': mp.__version__}
    for family, distribution, shape, loc, scale in [
        ('normal', stats.norm, 1., .4, 1.3),
        ('logistic', stats.logistic, 1., .4, 1.3),
        ('gumbel_r', stats.gumbel_r, 1., .4, 1.3),
        ('gumbel_l', stats.gumbel_l, 1., .4, 1.3),
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
        kwargs = {} if family in ('normal','logistic','gumbel_r','gumbel_l') else {'floc': 0.}
        result[family+'_fit'] = list(distribution.fit(data, **kwargs))
    intervals = [
        ('gumbel_r_tail',1000.,math.nextafter(1000.,math.inf),'gumbel_r',1.,1.),
        ('gumbel_r_left',-5.,math.nextafter(-5.,math.inf),'gumbel_r',1.,1.),
        ('gumbel_r_subnormal',0.,math.nextafter(0.,1.),'gumbel_r',1.,1.),
        ('gumbel_r_underflow',1e-300,2e-300,'gumbel_r',1.,1e300),
        ('logistic_tail',1000.,math.nextafter(1000.,math.inf),'logistic',1.,1.),
        ('logistic_subnormal',0.,math.nextafter(0.,1.),'logistic',1.,1.),
        ('logistic_underflow',1e-300,2e-300,'logistic',1.,1e300),
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
    exact, left, right, intervals = [.1,.3,.6], [.05], [.8], [(.2,.25),(.7,.75)]
    beta_data=stats.CensoredData(uncensored=exact,left=left,right=right,interval=intervals)
    result['beta_fit']=list(stats.beta.fit(beta_data,floc=0.,fscale=1.))
    with mp.workdps(700):
        result['beta_mixed']=float(sum(logpdf(x,'beta',.4,0.,1.,.6) for x in exact)
            +sum(logcdf(x,'beta',.4,0.,1.,.6) for x in left)
            +sum(logsf(x,'beta',.4,0.,1.,.6) for x in right)
            +sum(loginterval(a,b,'beta',.4,0.,1.,.6) for a,b in intervals))
        for name,a,b in [('beta_adjacent',.5,math.nextafter(.5,1.)),
                         ('beta_left',0.,1e-300),
                         ('beta_right',math.nextafter(1.,0.),1.),
                         ('beta_tail',.99,.999)]:
            result[name]=float(loginterval(a,b,'beta',.4,0.,1.,.6))
        result['beta_scaled_left']=float(logcdf(1e-300,'beta',.4,0.,1e300,.6))
        result['beta_compensated_right']=float(logsf(1e300,'beta',.4,1e-300,1e300,.6))
        result['beta_small_shape_sf']=float(logsf(1e-300,'beta',1e-20,0.,1e300,.6))
        result['beta_small_pdf']=float(logpdf(1e-300,'beta',.4,0.,1e300,.6))
        result['beta_large_pdf']=float(logpdf(.5,'beta',1000.,0.,1.,1000.))
        result['beta_large_adjacent']=float(loginterval(.5,math.nextafter(.5,1.),'beta',1000.,0.,1.,1000.))
    with mp.workdps(700):
        a,b=map(mp.mpf,(1e-20,1e308))
        z=mp.mpf(1e-300)/mp.mpf(1e300)
        # Leading-term relative correction is bounded by b*z < 1e-291.
        lp=a*mp.log(z)-mp.log(a)-mp.log(mp.beta(a,b))
        result['beta_tiny_huge_sf']=float(mp.log(-mp.expm1(lp)))
        b=mp.mpf(1e6)
        lp=a*mp.log(z)-mp.log(a)-mp.log(mp.beta(a,b))
        result['beta_tiny_large_sf']=float(mp.log(-mp.expm1(lp)))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
