# SPDX-License-Identifier: MIT
"""SciFort distributions and matrix procedures through its binary64 C interface."""
from ._generated import API
from .random import RandomState

from ._runtime import ArrayResult

class Scifort(API):
    """Configured library interfaces with explicit, independent RNG streams."""
    def rng(self, seed=12345):
        """Create a caller-owned native stream bound to this library instance."""
        return RandomState(self, seed)


__all__ = ['Scifort', 'ArrayResult', 'RandomState']
