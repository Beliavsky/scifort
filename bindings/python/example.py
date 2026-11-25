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


if __name__ == '__main__':
    main()
