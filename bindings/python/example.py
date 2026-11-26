# SPDX-License-Identifier: MIT
"""Run from any directory: python bindings/python/example.py (gfortran required)."""
from build_library import build
from scifort import Scifort
import numpy as np


def main():
    library, metadata = build()
    api = Scifort(library, dll_directories=metadata['dll_directories'])
    x = np.array([-1., 0., 1.])
    print('normal CDF:', api.normal.cdf(x))
    print('normal quantile:', api.normal.ppf(.975))
    print('Weibull survival:', api.weibull.sf(x, c=1.7))
    print('chi-square CDF:', api.chi2.cdf([0., 1., 4.], df=5.))
    print('gamma CDF:', api.gamma.cdf([0., 1., 4.], a=2.3))
    print('beta quantiles:', api.beta.ppf([.1, .5, .9], a=2.3, b=3.1))
    print('Student t normal limit:', api.t.cdf(x, df=float('inf')))
    print('beta broadcast CDF:', api.beta.cdf([[.25], [.75]], a=[1., 2., 3.], b=2.))
    print('Student t varying df:', api.t.cdf(x, df=[5., float('inf'), 10.]))
    output = np.empty_like(x)
    api.normal.logpdf(x, loc=.5, scale=2., out=output)
    print('normal log-PDF in caller buffer:', output)
    matrix = np.array([[4., 1.], [1., 3.]])
    factorization = api.linalg.cholesky_lower(matrix)
    if factorization.status:
        raise RuntimeError(f'Cholesky failed: {factorization.status}')
    solution = api.linalg.solve_spd_matrix(factorization.outputs['l'], np.array([[1.], [2.]]))
    if solution.status:
        raise RuntimeError(f'SPD solve failed: {solution.status}')
    print('SPD solution:', solution.outputs['x'])
    psd = api.linalg.psd_decompose(np.diag([4., 0., 1.]))
    print('PSD status/rank:', psd.status, psd.outputs['rank'])
    print('multivariate normal PDF:', api.multivariate_normal.pdf([0., 0.]))
    print('multivariate normal log-PDF:',
          api.multivariate_normal.logpdf([.5, -.2], cov=[[2., .3], [.3, 1.]]))
    print('multivariate Student t PDF:', api.multivariate_t.pdf([0., 0.], df=5.))
    print('multivariate Student t log-PDF:',
          api.multivariate_t.logpdf([.5, -.2], shape=[[2., .3], [.3, 1.]], df=5.))
    print('Dirichlet log-PDF, inferred coordinate:', api.dirichlet.logpdf([.2, .3], [2., 3., 4.]))
    print('Dirichlet PDF:', api.dirichlet.pdf([.2, .3, .5], [2., 3., 4.]))
    print('Dirichlet entropy:', api.dirichlet.entropy([2., 3., 4.]))
    print('Dirichlet mean/variance:', api.dirichlet.mean([2., 3., 4.]),
          api.dirichlet.var([2., 3., 4.]))
    covariance = np.empty((3, 3), order='F')
    api.dirichlet.cov([2., 3., 4.], out=covariance)
    print('Dirichlet covariance:', covariance)
    print('multinomial PMF/log-PMF:', api.multinomial.pmf([1., 3.], 4, [.2, .8]),
          api.multinomial.logpmf([1., 3.], 4, [.2, .8]))
    print('multinomial entropy/mean:', api.multinomial.entropy(4, [.2, .8]),
          api.multinomial.mean(4, [.2, .8]))
    print('multinomial covariance:', api.multinomial.cov(4, [.2, .8]))
    print('multivariate hypergeometric PMF/log-PMF:',
          api.multivariate_hypergeom.pmf([1, 3], [3, 5], 4),
          api.multivariate_hypergeom.logpmf([1, 3], [3, 5], 4))
    print('multivariate hypergeometric mean/variance:',
          api.multivariate_hypergeom.mean([3, 5], 4),
          api.multivariate_hypergeom.var([3, 5], 4))
    print('multivariate hypergeometric covariance:',
          api.multivariate_hypergeom.cov([3, 5], 4))

    print('Dirichlet-multinomial PMF/log-PMF:',
          api.dirichlet_multinomial.pmf([1, 3], [2., 3.], 4),
          api.dirichlet_multinomial.logpmf([1, 3], [2., 3.], 4))
    print('Dirichlet-multinomial mean/variance/covariance:',
          api.dirichlet_multinomial.mean([2., 3.], 4),
          api.dirichlet_multinomial.var([2., 3.], 4),
          api.dirichlet_multinomial.cov([2., 3.], 4))

    scale = [[2., .4], [.4, 1.]]
    point = [[1.5, .2], [.2, 1.1]]
    for name in ('wishart', 'invwishart'):
        distribution = getattr(api, name)
        print(name, 'PDF/log-PDF/entropy:', distribution.pdf(point, 7., scale),
              distribution.logpdf(point, 7., scale), distribution.entropy(7., scale))
        print(name, 'mean/mode/variance:', distribution.mean(7., scale),
              distribution.mode(7., scale), distribution.var(7., scale))

    point = [[.1, .2, .3], [.4, .5, .6]]
    print('matrix normal PDF/log-PDF/entropy:', api.matrix_normal.pdf(point),
          api.matrix_normal.logpdf(point), api.matrix_normal.entropy(np.eye(2), np.eye(3)))
    print('matrix normal entropy with defaults:', api.matrix_normal.entropy())
    print('matrix t PDF/log-PDF:', api.matrix_t.pdf(point, df=3.),
          api.matrix_t.logpdf(point, df=3.))


if __name__ == '__main__':
    main()
