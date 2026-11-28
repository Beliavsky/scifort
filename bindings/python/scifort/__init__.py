# SPDX-License-Identifier: MIT
"""SciFort distributions and matrix procedures through its binary64 C interface."""
from ._generated import API
from .censored import CensoredLikelihoods
from .random import RandomState
from .kde import GaussianKDE
from .survival import logrank, LogrankResult
from .ecdf import ecdf, EmpiricalDistributionFunction, ECDFResult, ECDFConfidenceInterval

from ._runtime import ArrayResult

class Scifort(API):
    """Configured library interfaces with explicit, independent RNG streams."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.censored = CensoredLikelihoods(self)

    def rng(self, seed=12345):
        """Create a caller-owned native stream bound to this library instance."""
        return RandomState(self, seed)

    def logrank(self, x, y, *, x_right=None, y_right=None, alternative='two-sided'):
        """Compare ordinary/right-censored survival samples with explicit alternatives."""
        return logrank(self, x, y, x_right=x_right, y_right=y_right, alternative=alternative)

    def ecdf(self, sample, *, right=None):
        """Create independent empirical CDF/SF owners; optionally include right censoring."""
        return ecdf(self, sample, right=right)

    def gaussian_kde(self, dataset, bw_method=None, weights=None):
        """Create an owning KDE; use close() or a context manager to release it."""
        return GaussianKDE(self, dataset, bw_method=bw_method, weights=weights)


__all__ = ['Scifort', 'CensoredLikelihoods', 'ArrayResult', 'RandomState', 'GaussianKDE',
           'EmpiricalDistributionFunction', 'ECDFResult', 'ECDFConfidenceInterval', 'LogrankResult']
