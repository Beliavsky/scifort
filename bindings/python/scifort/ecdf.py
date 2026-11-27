# SPDX-License-Identifier: MIT
"""Owned empirical CDF/SF functions, including Kaplan-Meier estimates."""
import numpy as np
from ._runtime import Distribution


class EmpiricalDistributionFunction:
    """An independently owned empirical step function; use close or a with block."""
    def __init__(self, group, handle):
        self._group = group
        self._handle = handle
        result = group.info(handle)
        if result.status:
            handle.close()
            raise ValueError(f'EDF information query failed: status {result.status}')
        self._info = result.outputs

    @property
    def closed(self):
        return self._handle.closed

    def close(self):
        self._handle.close()

    def __enter__(self):
        if self.closed:
            raise ValueError('EDF is closed')
        return self

    def __exit__(self, *exc):
        self.close()

    def __copy__(self):
        raise TypeError('EDF ownership cannot be copied')

    def __deepcopy__(self, memo):
        return self.__copy__()

    def __reduce_ex__(self, protocol):
        raise TypeError('EDF objects cannot be pickled')

    def _live_info(self, name):
        if self.closed:
            raise ValueError('EDF is closed')
        return self._info[name]

    def get_values(self, *, out=None):
        """Copy knots and probabilities; out may supply disjoint float64 vectors."""
        with self._handle._lock:
            result = self._group.values(self._handle, m=self._live_info('count'), out=out)
            if result.status:
                raise ValueError(f'EDF array retrieval failed: status {result.status}')
            return result.outputs['quantiles'], result.outputs['probabilities']

    @property
    def quantiles(self):
        return self.get_values()[0]

    @property
    def probabilities(self):
        return self.get_values()[1]

    def evaluate(self, x, *, out=None):
        """Evaluate a scalar or vector; out is a float64 vector of matching length."""
        x = Distribution._real_array(x, 'x')
        scalar = x.ndim == 0
        if x.ndim > 1:
            raise ValueError('x must be a scalar or vector')
        result = self._group.evaluate(self._handle, x.reshape(-1),
                                      out=None if out is None else {'values': out})
        if result.status:
            raise ValueError(f'EDF evaluation failed: status {result.status}')
        values = result.outputs['values']
        return float(values[0]) if scalar else values

    __call__ = evaluate

    def confidence_interval(self, confidence_level=.95, method='linear'):
        """Return independent Greenwood bounds; nested intervals are unsupported."""
        if method not in ('linear', 'log-log'):
            raise ValueError('method must be linear or log-log')
        result = self._group.interval(self._handle, confidence_level, int(method == 'log-log'))
        return ECDFConfidenceInterval._from_result(self._group, result,
                                                  ('low_handle', 'high_handle'), ('low', 'high'))


class _OwnedPair:
    @classmethod
    def _from_result(cls, group, result, names, attributes):
        handles = [result.outputs[name] for name in names]
        try:
            if result.status:
                raise ValueError(f'ECDF operation failed: status {result.status}')
            if any(handle is None for handle in handles):
                raise RuntimeError('ECDF operation returned a null owner')
            pair = cls()
            pair._functions = []
            for attribute, handle in zip(attributes, handles):
                function = EmpiricalDistributionFunction(group, handle)
                pair._functions.append(function)
                setattr(pair, attribute, function)
            return pair
        except Exception:
            for handle in handles:
                if handle is not None:
                    handle.close()
            raise

    def close(self):
        """Close both functions; retained references to either then become closed."""
        for function in self._functions:
            function.close()

    def __enter__(self):
        if any(function.closed for function in self._functions):
            raise ValueError('ECDF pair contains a closed function')
        return self

    def __exit__(self, *exc):
        self.close()

    def __copy__(self):
        raise TypeError('ECDF pair ownership cannot be copied')

    def __deepcopy__(self, memo):
        return self.__copy__()

    def __reduce_ex__(self, protocol):
        raise TypeError('ECDF pairs cannot be pickled')


class ECDFResult(_OwnedPair):
    """Independent cdf and sf owners. Closing one does not close the other."""


class ECDFConfidenceInterval(_OwnedPair):
    """Independent low and high owners, which survive closure of their source."""


def ecdf(api, sample, *, right=None):
    """Construct from finite events and optional finite right-censoring times."""
    events = Distribution._real_array(sample, 'sample')
    censoring = Distribution._real_array([] if right is None else right, 'right')
    if events.ndim != 1 or censoring.ndim != 1:
        raise ValueError('sample and right must be vectors')
    if np.may_share_memory(events, censoring):
        censoring = censoring.copy()
    result = api.ecdf_bridge.create(events, censoring)
    return ECDFResult._from_result(api.ecdf_bridge, result,
                                  ('cdf_handle', 'sf_handle'), ('cdf', 'sf'))
