# SPDX-License-Identifier: MIT
"""Owned Gaussian KDE objects backed by reviewed opaque-handle adapters."""
import numpy as np
from ._runtime import Distribution


class GaussianKDE:
    """Copy observations into an estimator; close explicitly or use a with block."""
    def __init__(self, api, dataset, bw_method=None, weights=None):
        data = Distribution._real_array(dataset, 'dataset')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.ndim != 2:
            raise ValueError('dataset must have shape (n,) or (d,n)')
        w = (np.ones(data.shape[1], dtype=np.float64) if weights is None else weights)
        if bw_method is None:
            mode, factor = 0, 1.
        elif isinstance(bw_method, str):
            if bw_method not in ('scott', 'silverman'):
                raise ValueError('bw_method must be scott, silverman, or a positive scalar factor')
            mode, factor = int(bw_method == 'silverman'), 1.
        else:
            factor_array = Distribution._real_array(bw_method, 'bw_method')
            if factor_array.ndim or isinstance(bw_method, (bool, np.bool_)):
                raise TypeError('bw_method must be a scalar bandwidth factor')
            mode, factor = 2, float(factor_array)
        result = api.kde_bridge.create(data, w, mode, factor)
        self._adopt_result(api.kde_bridge, result, 'handle', 'initialization')

    def _adopt_result(self, group, result, name, operation):
        self._group = group
        handle = result.outputs[name]
        if result.status:
            if handle is not None:
                handle.close()
            raise ValueError(f'Gaussian KDE {operation} failed: status {result.status}')
        if handle is None:
            raise RuntimeError(f'Gaussian KDE {operation} returned a null handle')
        self._handle = handle
        try:
            self._info = self._query_info()
        except Exception:
            self.close()
            raise

    def _query_info(self):
        result = self._group.info(self._handle)
        if result.status:
            raise RuntimeError(f'Gaussian KDE information query failed: status {result.status}')
        return result.outputs

    def _get_info(self, name):
        with self._handle._lock:
            if self.closed:
                raise ValueError('Gaussian KDE is closed')
            return self._info[name]

    @property
    def closed(self):
        return self._handle.closed

    def close(self):
        """Release native storage; repeated close calls are harmless."""
        self._handle.close()

    def __enter__(self):
        if self.closed:
            raise ValueError('Gaussian KDE is closed')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __copy__(self):
        raise TypeError('Gaussian KDE ownership cannot be copied')

    def __deepcopy__(self, memo):
        return self.__copy__()

    def __reduce_ex__(self, protocol):
        raise TypeError('Gaussian KDE objects cannot be pickled')

    @property
    def d(self):
        return self._get_info('d')

    @property
    def n(self):
        return self._get_info('n')

    @property
    def neff(self):
        return self._get_info('neff')

    @property
    def factor(self):
        return self._get_info('factor')

    def scotts_factor(self):
        return self._get_info('scott')

    def silverman_factor(self):
        return self._get_info('silverman')

    def _evaluate(self, method, points, out):
        with self._handle._lock:
            d = self.d
            x = Distribution._real_array(points, 'points')
            if x.ndim == 0 and d == 1:
                x = x.reshape(1, 1)
            elif x.ndim == 1:
                x = x.reshape(1, -1) if d == 1 else x.reshape(-1, 1)
            if x.ndim != 2 or x.shape[0] != d:
                raise ValueError(f'points must have shape ({d},m)')
            buffers = None if out is None else {'values': out}
            result = getattr(self._group, method)(self._handle, x, out=buffers)
            if result.status:
                raise ValueError(f'Gaussian KDE {method} failed: status {result.status}')
            return result.outputs['values']

    def pdf(self, points, *, out=None):
        """Evaluate columns; always return an array of length m."""
        return self._evaluate('pdf', points, out)

    evaluate = pdf
    __call__ = pdf

    def logpdf(self, points, *, out=None):
        """Evaluate native log densities directly, without logging the PDF."""
        return self._evaluate('logpdf', points, out)

    @staticmethod
    def _scalar_result(result, operation):
        if result.status:
            raise ValueError(f'Gaussian KDE {operation} failed: status {result.status}')
        return result.outputs['value']

    def integrate_gaussian(self, mean, cov):
        """Integrate the product with a Gaussian density of matching dimension."""
        mean = np.atleast_1d(Distribution._real_array(mean, 'mean'))
        cov = np.atleast_2d(Distribution._real_array(cov, 'cov'))
        result = self._group.integrate_gaussian(self._handle, mean, cov)
        return self._scalar_result(result, 'Gaussian integration')

    def integrate_box_1d(self, low, high):
        """Integrate an interval for a univariate KDE; reversed limits change the sign."""
        result = self._group.integrate_box_1d(self._handle, low, high)
        return self._scalar_result(result, 'interval integration')

    def integrate_box(self, low, high, *, state, maxpts=0):
        """Integrate a box with explicit QMC state; maxpts=0 selects native defaults."""
        low = np.atleast_1d(Distribution._real_array(low, 'low'))
        high = np.atleast_1d(Distribution._real_array(high, 'high'))
        if np.may_share_memory(low, high):
            high = np.array(high, copy=True, order='F')
        result = self._group.integrate_box(handle=self._handle, low=low, high=high,
                                           state=state, maxpts=maxpts)
        return self._scalar_result(result, 'box integration')

    def integrate_kde(self, other):
        """Integrate two KDEs; acquire both handle locks in the runtime's common order."""
        if not isinstance(other, GaussianKDE):
            raise TypeError('other must be a GaussianKDE from the same library instance')
        # Do not hold a facade lock first: the runtime orders/deduplicates both locks.
        result = self._group.integrate_kde(self._handle, other._handle)
        return self._scalar_result(result, 'KDE product integration')

    def resample(self, state, size=None, *, out=None):
        """Draw columns using explicit RNG state; default size is int(neff)."""
        with self._handle._lock:
            if size is None:
                size = int(self.neff)
            result = self._group.resample(handle=self._handle, state=state, d=self.d,
                                          n_samples=size, out=None if out is None else {'samples': out})
            if result.status:
                raise ValueError(f'Gaussian KDE resampling failed: status {result.status}')
            return result.outputs['samples']

    def marginal(self, dimensions):
        """Return an independent KDE for distinct zero-based (possibly negative) axes."""
        with self._handle._lock:
            d = self.d
            axes = np.atleast_1d(np.asarray(dimensions))
            if axes.ndim != 1 or not axes.size:
                raise ValueError('dimensions must be a nonempty coordinate vector')
            if axes.dtype.kind not in 'iu':
                raise TypeError('dimensions must contain integers')
            if np.any(axes < -d) or np.any(axes >= d):
                raise ValueError('dimensions must lie in [-d,d)')
            axes = axes.astype(np.int64)
            axes = np.where(axes < 0, axes+d, axes)
            if np.unique(axes).size != axes.size:
                raise ValueError('dimensions must be distinct after normalizing negative indices')
            result = self._group.marginal(self._handle, axes+1)
        child = type(self).__new__(type(self))
        child._adopt_result(self._group, result, 'marginal_handle', 'marginalization')
        return child
