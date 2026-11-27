# SPDX-License-Identifier: MIT
"""SciFort distributions and matrix procedures through its binary64 C interface."""
from ._generated import API
from .random import RandomState
from .kde import GaussianKDE

from ._runtime import ArrayResult

class Scifort(API):
    """Configured library interfaces with explicit, independent RNG streams."""
    def rng(self, seed=12345):
        """Create a caller-owned native stream bound to this library instance."""
        return RandomState(self, seed)

    def gaussian_kde(self, dataset, bw_method=None, weights=None):
        """Create an owning KDE; use close() or a context manager to release it."""
        return GaussianKDE(self, dataset, bw_method=bw_method, weights=weights)


__all__ = ['Scifort', 'ArrayResult', 'RandomState', 'GaussianKDE']
