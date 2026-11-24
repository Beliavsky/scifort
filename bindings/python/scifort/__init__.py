# SPDX-License-Identifier: MIT
"""SciFort distributions and matrix procedures through its binary64 C interface."""
from ._generated import API as Scifort

from ._runtime import ArrayResult

__all__ = ['Scifort', 'ArrayResult']
