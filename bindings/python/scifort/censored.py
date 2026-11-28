# SPDX-License-Identifier: MIT
"""Copied-data likelihoods with the native SciFort censoring conventions."""
import numpy as np
from ._runtime import Distribution


class CensoredLikelihoods:
    """Censored likelihoods for six reviewed distribution families."""

    def __init__(self, api):
        self._bridge = api.censored_bridge

    def _evaluate(self, method, uncensored, left, right, interval, parameters):
        arrays = []
        for name, values in (('uncensored', uncensored), ('left', left), ('right', right)):
            array = Distribution._real_array([] if values is None else values, name)
            if array.ndim != 1:
                raise ValueError(f'{name} must be a vector')
            arrays.append(array)
        bounds = Distribution._real_array(np.empty((0, 2)) if interval is None else interval,
                                          'interval')
        if bounds.shape == (0,):
            bounds = bounds.reshape(0, 2)
        if bounds.ndim != 2 or bounds.shape[1] != 2:
            raise ValueError('interval must have shape (n_interval, 2)')
        arrays.extend((bounds[:, 0], bounds[:, 1]))
        # Read-only aliases are harmless mathematically; supply disjoint C buffers.
        for index, array in enumerate(arrays):
            if any(np.may_share_memory(array, previous) for previous in arrays[:index]):
                arrays[index] = array.copy()
        result = getattr(self._bridge, method)(*arrays, **parameters)
        if result.status:
            raise ValueError(f'Censored likelihood failed: status {result.status}')
        return float(result.outputs['value'])

    def normal_loglikelihood(self, uncensored=None, *, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate normal log likelihood; empty data contributes zero."""
        return self._evaluate('normal_loglikelihood', uncensored, left, right, interval,
                              dict(loc=loc, scale=scale))

    def normal_nnlf(self, uncensored=None, *, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate normal negative log likelihood; empty data contributes zero."""
        return self._evaluate('normal_nnlf', uncensored, left, right, interval,
                              dict(loc=loc, scale=scale))

    def exponential_loglikelihood(self, uncensored=None, *, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate exponential log likelihood; empty data contributes zero."""
        return self._evaluate('exponential_loglikelihood', uncensored, left, right, interval,
                              dict(loc=loc, scale=scale))

    def exponential_nnlf(self, uncensored=None, *, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate exponential negative log likelihood; empty data contributes zero."""
        return self._evaluate('exponential_nnlf', uncensored, left, right, interval,
                              dict(loc=loc, scale=scale))

    def weibull_loglikelihood(self, uncensored=None, *, c, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate weibull log likelihood; empty data contributes zero."""
        return self._evaluate('weibull_loglikelihood', uncensored, left, right, interval,
                              dict(c=c, loc=loc, scale=scale))

    def weibull_nnlf(self, uncensored=None, *, c, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate weibull negative log likelihood; empty data contributes zero."""
        return self._evaluate('weibull_nnlf', uncensored, left, right, interval,
                              dict(c=c, loc=loc, scale=scale))

    def lognormal_loglikelihood(self, uncensored=None, *, s, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate lognormal log likelihood; empty data contributes zero."""
        return self._evaluate('lognormal_loglikelihood', uncensored, left, right, interval,
                              dict(s=s, loc=loc, scale=scale))

    def lognormal_nnlf(self, uncensored=None, *, s, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate lognormal negative log likelihood; empty data contributes zero."""
        return self._evaluate('lognormal_nnlf', uncensored, left, right, interval,
                              dict(s=s, loc=loc, scale=scale))

    def gamma_loglikelihood(self, uncensored=None, *, a, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate gamma log likelihood; empty data contributes zero."""
        return self._evaluate('gamma_loglikelihood', uncensored, left, right, interval,
                              dict(a=a, loc=loc, scale=scale))

    def gamma_nnlf(self, uncensored=None, *, a, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate gamma negative log likelihood; empty data contributes zero."""
        return self._evaluate('gamma_nnlf', uncensored, left, right, interval,
                              dict(a=a, loc=loc, scale=scale))

    def beta_loglikelihood(self, uncensored=None, *, a, b, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate beta log likelihood; empty data contributes zero."""
        return self._evaluate('beta_loglikelihood', uncensored, left, right, interval,
                              dict(a=a, b=b, loc=loc, scale=scale))

    def beta_nnlf(self, uncensored=None, *, a, b, left=None, right=None,
            interval=None, loc=0.0, scale=1.0):
        """Evaluate beta negative log likelihood; empty data contributes zero."""
        return self._evaluate('beta_nnlf', uncensored, left, right, interval,
                              dict(a=a, b=b, loc=loc, scale=scale))

