# SPDX-License-Identifier: MIT
"""Stateless survival comparisons through reviewed flat array adapters."""
from typing import NamedTuple
import numpy as np
from ._runtime import Distribution


class LogrankResult(NamedTuple):
    """Signed standard-normal logrank statistic and its tail probability."""
    statistic: float
    pvalue: float


def logrank(api, x, y, *, x_right=None, y_right=None, alternative='two-sided'):
    """Compare finite event/censor vectors; undefined comparisons return NaNs."""
    alternatives = ('two-sided', 'less', 'greater')
    if not isinstance(alternative, str) or alternative not in alternatives:
        raise ValueError('alternative must be two-sided, less, or greater')
    arrays = []
    for name, values in (('x', x), ('x_right', [] if x_right is None else x_right),
                         ('y', y), ('y_right', [] if y_right is None else y_right)):
        array = Distribution._real_array(values, name)
        if array.ndim != 1:
            raise ValueError(f'{name} must be a vector')
        if any(np.may_share_memory(array, previous) for previous in arrays):
            array = array.copy()
        arrays.append(array)
    result = api.survival.logrank(*arrays, alternative=alternatives.index(alternative))
    if result.status:
        raise ValueError(f'Logrank comparison failed: status {result.status}')
    return LogrankResult(result.outputs['statistic'], result.outputs['pvalue'])
