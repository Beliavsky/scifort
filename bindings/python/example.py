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
    output = np.empty_like(x)
    api.normal.logpdf(x, loc=.5, scale=2., out=output)
    print('normal log-PDF in caller buffer:', output)


if __name__ == '__main__':
    main()
