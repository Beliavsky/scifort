# SPDX-License-Identifier: MIT
"""Test the generated package against SciPy and checked/optimized native libraries."""
import ctypes as C
import inspect
import math
import os
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from build_library import build, FAMILIES
import build_library as builder
from scifort import Scifort

SCIPY_NAMES = dict(normal='norm', exponential='expon', weibull='weibull_min', lognormal='lognorm')
SHAPES = dict(weibull=(1.7,), lognormal=(.8,), gamma=(2.3,), beta=(2.3, 3.1), t=(5.,), chi2=(5.,))
METHODS = ('pdf', 'logpdf', 'cdf', 'sf', 'logcdf', 'logsf', 'ppf', 'isf')


class BindingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran is required')
        cls.apis = []
        for checked in (True, False):
            path, metadata = build(checked=checked)
            cls.apis.append((checked, Scifort(path, dll_directories=metadata['dll_directories'])))

    def test_all_methods_scalar_array_and_scipy_agreement(self):
        for checked, api in self.apis:
            for family in FAMILIES:
                group = getattr(api, family)
                scipy_group = getattr(stats, SCIPY_NAMES.get(family, family))
                shape = SHAPES.get(family, ())
                for method in METHODS:
                    with self.subTest(checked=checked, family=family, method=method):
                        func, reference = getattr(group, method), getattr(scipy_group, method)
                        values = (np.array([0., 1e-12, .1, .5, .9, 1-1e-12, 1., -.1, 1.1, np.nan])
                                  if method in ('ppf', 'isf') else
                                  np.array([-np.inf, -8., -.1, 0., .25, 1., 8., 40., np.inf, np.nan]))
                        actual = func(values, *shape, loc=.25, scale=1.5)
                        with np.errstate(invalid='ignore'):
                            expected = reference(values, *shape, loc=.25, scale=1.5)
                        # SciPy 1.15.3 returns NaN for gamma/chi-square/Weibull density at +inf.
                        # Preserve SciFort's established mathematical endpoint limits.
                        if family in ('weibull', 'gamma', 'chi2') and method in ('pdf', 'logpdf'):
                            endpoint = np.isposinf(values)
                            limit = 0. if method == 'pdf' else -np.inf
                            np.testing.assert_array_equal(actual[endpoint], limit)
                            expected[endpoint] = limit
                        # SciPy 1.15.3 t(5) quantiles at .1/.9 miss ~3e-11.
                        # Independent 80-digit mpmath inversion of the incomplete-beta
                        # CDF gives |t| = 1.4758840488244810785464589478088632618.
                        # Retain the strict tolerances using this reference at those points.
                        if family == 't' and method in ('ppf', 'isf'):
                            quantile = 1.4758840488244810785464589478088632618
                            for probability, sign in ((.1, -1.), (.9, 1.)):
                                if method == 'isf':
                                    sign = -sign
                                index = values == probability
                                expected[index] = .25 + sign*1.5*quantile
                                np.testing.assert_allclose(actual[index], expected[index],
                                                           rtol=2e-15, atol=0)
                        # Density/probability comparisons use relative errors, including small tails.
                        if method in ('pdf', 'cdf', 'sf'):
                            np.testing.assert_allclose(actual, expected, rtol=4e-12, atol=0, equal_nan=True)
                        else:
                            np.testing.assert_allclose(actual, expected, rtol=4e-12, atol=3e-13, equal_nan=True)
                        scalars = [func(float(x), *shape, loc=.25, scale=1.5) for x in values]
                        np.testing.assert_array_equal(actual, scalars)
                        output = np.empty_like(actual)
                        self.assertIs(func(values, *shape, loc=.25, scale=1.5, out=output), output)
                        np.testing.assert_array_equal(actual, output)
                        self.assertIsInstance(func(.5, *shape), float)
                        np.testing.assert_allclose(func(.5, *shape), reference(.5, *shape),
                                                   rtol=4e-12, atol=3e-13)

    def test_shapes_conversion_alignment_and_empty_inputs(self):
        for _, api in self.apis:
            func = api.normal.cdf
            for x in ([0, 1, 2], np.arange(12).reshape(3, 4)[:, ::2],
                      np.arange(6, dtype=np.float32).reshape(2, 3),
                      np.array([1., 2.], dtype='>f8'), np.zeros((2, 0, 3))):
                y = func(x)
                self.assertEqual(y.shape, np.shape(x))
                self.assertEqual(y.dtype, np.float64)
                np.testing.assert_allclose(y, stats.norm.cdf(x), rtol=3e-13, atol=0)
            unaligned = np.ndarray((2,), dtype=np.float64, buffer=bytearray(17), offset=1)
            unaligned[:] = [0., 1.]
            self.assertFalse(unaligned.flags.aligned)
            np.testing.assert_allclose(func(unaligned), stats.norm.cdf(unaligned), rtol=3e-13)
            self.assertIn('s', inspect.signature(api.lognormal.cdf).parameters)
            self.assertEqual(inspect.signature(func).parameters['scale'].default, 1.)

    def test_invalid_parameters_status_and_output_safety(self):
        for _, api in self.apis:
            x = np.array([0., 1., 2.])
            for family in FAMILIES:
                func = getattr(api, family).cdf
                shape = SHAPES.get(family, ())
                for params in (dict(scale=0.), dict(scale=-1.), dict(scale=np.inf),
                               dict(loc=np.nan), dict(loc=np.inf)):
                    out = np.full_like(x, 99.)
                    with self.assertRaises(ValueError):
                        func(x, *shape, out=out, **params)
                    np.testing.assert_array_equal(out, 99.)
                    self.assertTrue(math.isnan(func(.5, *shape, **params)))
                    with self.assertRaises(ValueError):
                        func([], *shape, **params)
                for index in range(len(shape)):
                    for bad in (0., -1., np.nan, -np.inf, *(() if family == 't' else (np.inf,))):
                        invalid = list(shape)
                        invalid[index] = bad
                        with self.assertRaises(ValueError):
                            func(x, *invalid)
            np.testing.assert_array_equal(api.beta.cdf(x, 2., 3., scale=[1.]),
                                          api.beta.cdf(x, 2., 3.))
            for bad in ([1+2j], ['1'], [object()]):
                with self.assertRaises(TypeError):
                    api.normal.cdf(bad)
            invalid_outputs = [np.empty(3, dtype=np.float32), np.empty(4),
                               np.empty(6)[::2], x, x[::-1]]
            readonly = np.empty(3)
            readonly.flags.writeable = False
            invalid_outputs += [readonly, np.ndarray((3,), dtype=np.float64,
                                                    buffer=bytearray(25), offset=1),
                                np.empty(3, dtype='>f8')]
            for out in invalid_outputs:
                with self.assertRaises((TypeError, ValueError)):
                    api.normal.cdf(x, out=out)
            with self.assertRaises(ValueError):
                api.normal.cdf(0., out=np.empty(()))

    def test_broadcast_parameters(self):
        for checked, api in self.apis:
            for family, reference in (('normal', stats.norm), ('gamma', stats.gamma)):
                group = getattr(api, family)
                for method in METHODS:
                    with self.subTest(checked=checked, family=family, method=method):
                        x = (np.array([[.01], [.5], [.99]]) if method in ('ppf', 'isf')
                             else np.array([[-3.], [.2], [8.]]))
                        loc = np.array([-.5, 0., .5], dtype='>f8')
                        scale = np.array([1., 2., 3.], dtype=np.float32)
                        args = (np.array([1., 2., 3.]),) if family == 'gamma' else ()
                        func = getattr(group, method)
                        actual = func(x, *args, loc=loc, scale=scale)
                        # Match this binary64 interface's documented conversion. SciPy
                        # otherwise keeps float32 precision in log(scale) for logpdf.
                        expected = getattr(reference, method)(x, *args, loc=loc.astype(np.float64),
                                                               scale=scale.astype(np.float64))
                        np.testing.assert_allclose(actual, expected, rtol=4e-12, atol=3e-13)
                        # Every element uses the unchanged native scalar function.
                        for i, j in np.ndindex(actual.shape):
                            shapes = (args[0][j],) if args else ()
                            self.assertEqual(actual[i,j], func(x[i,0], *shapes,
                                                              loc=float(loc[j]), scale=float(scale[j])))
                        out = np.empty_like(actual)
                        self.assertIs(func(x, *args, loc=loc, scale=scale, out=out), out)
                        np.testing.assert_array_equal(out, actual)
                        np.testing.assert_allclose(func(.5, *args, scale=scale),
                                                   getattr(reference, method)(.5, *args, scale=scale.astype(np.float64)),
                                                   rtol=4e-12, atol=3e-13)
                # Array parameters follow native per-element behavior, including invalids.
                for method in METHODS:
                    func = getattr(group, method)
                    scales = np.array([1., 0., -1., np.nan, np.inf])
                    shapes = (2.,) if family == 'gamma' else ()
                    actual = func(.5, *shapes, scale=scales)
                    expected = [func(.5, *shapes, scale=float(v)) for v in scales]
                    np.testing.assert_array_equal(actual, expected)
                for method in METHODS:
                    func = getattr(group, method)
                    points = ([0., 1., -.1, 1.1, np.nan] if method in ('ppf', 'isf')
                              else [-np.inf, 0., np.inf, np.nan, 1e100])
                    scales = np.arange(1., 6.)
                    shapes = (2.,) if family == 'gamma' else ()
                    np.testing.assert_array_equal(func(points, *shapes, scale=scales),
                        [func(v, *shapes, scale=float(t)) for v, t in zip(points, scales)])
                probabilities = np.array([[1e-12], [.5], [1-1e-12]])
                shapes = ([1., 2., 3.],) if family == 'gamma' else ()
                # Zero gamma location avoids rounding a tiny quantile into an offset.
                locations = [0., 0., 0.] if family == 'gamma' else [-.5, 0., .5]
                quantiles = group.ppf(probabilities, *shapes, loc=locations, scale=[1., 2., 3.])
                recovered = group.cdf(quantiles, *shapes, loc=locations, scale=[1., 2., 3.])
                np.testing.assert_allclose(recovered, np.broadcast_to(probabilities, (3,3)),
                                           rtol=2e-10, atol=0)
                empty = group.cdf(np.empty((0,1)), *(([2., -1., np.nan],) if family == 'gamma' else ()),
                                  scale=[1., 0., np.nan])
                self.assertEqual(empty.shape, (0,3))
            for label in ('loc', 'scale'):
                a = np.array([1., 2., 3.])
                with self.assertRaises(ValueError):
                    api.normal.cdf(0., out=a, **{label:a})
                np.testing.assert_array_equal(a, [1., 2., 3.])
            a = np.array([1., 2., 3.])
            with self.assertRaises(ValueError):
                api.gamma.cdf(.5, a, out=a)
            # Zero-dimensional parameter views also count as inputs for alias checks.
            with self.assertRaises(ValueError):
                api.normal.cdf(a, loc=a[:1].reshape(()), out=a)
            for params in (dict(loc=np.ones(4)), dict(scale=[1+2j])):
                with self.assertRaises((ValueError, TypeError)):
                    api.normal.cdf(np.ones(3), **params)
            for out in (np.empty((3,1)), np.empty(3, dtype='f4'), np.empty(6)[::2]):
                with self.assertRaises((ValueError, TypeError)):
                    api.normal.cdf(.5, loc=[0., 1., 2.], out=out)
            unaligned = np.ndarray((3,), dtype=np.float64, buffer=bytearray(25), offset=1)
            unaligned[:] = [1., 2., 3.]
            np.testing.assert_allclose(api.normal.cdf(.5, scale=unaligned),
                                       stats.norm.cdf(.5, scale=unaligned), rtol=3e-13)

    def test_remaining_family_broadcasts(self):
        shape_arrays = dict(weibull=([1.3, 2., 3.],), lognormal=([.5, .8, 1.2],),
                            beta=([2., 3., 4.], [3., 4., 5.]), t=([8., 12., 20.],),
                            chi2=([2., 4., 6.],))
        for checked, api in self.apis:
            for family in FAMILIES:
                if family in ('normal', 'gamma'):
                    continue
                group = getattr(api, family)
                reference = getattr(stats, SCIPY_NAMES.get(family, family))
                args = tuple(np.asarray(a) for a in shape_arrays.get(family, ()))
                for method in METHODS:
                    with self.subTest(checked=checked, family=family, method=method):
                        func = getattr(group, method)
                        points = (np.array([[.2], [.5], [.8]]) if method in ('ppf', 'isf')
                                  else np.array([[-.5], [.5], [3.]]))
                        loc = np.array([-.25, 0., .25], dtype='>f8')
                        scale = np.array([1., 2., 3.], dtype=np.float32)
                        actual = func(points, *args, loc=loc, scale=scale)
                        expected = np.array([[func(float(v), *(float(a[j]) for a in args),
                                                   loc=float(loc[j]), scale=float(scale[j]))
                                              for j in range(3)] for v in points[:,0]])
                        np.testing.assert_array_equal(actual, expected)
                        if family == 't' and method in ('ppf', 'isf'):
                            # SciPy 1.15.3 inverse-t has known reference inaccuracies.
                            # Check its forward probability and exact native agreement.
                            forward = reference.cdf if method == 'ppf' else reference.sf
                            np.testing.assert_allclose(forward(actual, *args, loc=loc.astype('f8'),
                                                              scale=scale.astype('f8')),
                                                       np.broadcast_to(points, (3,3)),
                                                       rtol=4e-12, atol=0)
                        else:
                            expected = getattr(reference, method)(points, *args,
                                loc=loc.astype('f8'), scale=scale.astype('f8'))
                            np.testing.assert_allclose(actual, expected, rtol=4e-12, atol=3e-13)
                        out = np.empty_like(actual)
                        self.assertIs(func(points, *args, loc=loc, scale=scale, out=out), out)
                        np.testing.assert_array_equal(out, actual)
                        # All shape arrays and loc/scale broadcast from scalar points too.
                        actual = func(.5, *args, loc=loc, scale=scale)
                        expected = [func(.5, *(float(a[j]) for a in args),
                                         loc=float(loc[j]), scale=float(scale[j])) for j in range(3)]
                        np.testing.assert_array_equal(actual, expected)
                        edge = np.array([0., 1., -.1, 1.1, np.nan]) if method in ('ppf', 'isf') else np.array(
                            [-np.inf, 0., 1., np.inf, np.nan])
                        shapes = SHAPES.get(family, ())
                        scales = np.arange(1., 6.)
                        np.testing.assert_array_equal(func(edge, *shapes, scale=scales),
                            [func(float(v), *shapes, scale=float(t)) for v, t in zip(edge, scales)])
                self.assertEqual(group.cdf(np.empty((0,1)), *args, scale=[1., 0., np.nan]).shape, (0,3))
                # Ordinary round trips avoid representational loss at shifted support endpoints.
                probabilities = np.array([[.2], [.5], [.8]])
                for inverse, forward in (('ppf', 'cdf'), ('isf', 'sf')):
                    points = getattr(group, inverse)(probabilities, *args, scale=[1., 2., 3.])
                    actual = getattr(group, forward)(points, *args, scale=[1., 2., 3.])
                    np.testing.assert_allclose(actual, np.broadcast_to(probabilities, (3,3)),
                                               rtol=4e-12, atol=0)

    def test_remaining_family_invalid_elements_and_buffers(self):
        for checked, api in self.apis:
            for family in FAMILIES:
                if family in ('normal', 'gamma'):
                    continue
                group = getattr(api, family)
                shape = SHAPES.get(family, ())
                for method in METHODS:
                    func = getattr(group, method)
                    for name, bad in (('loc', [-np.inf, np.inf, np.nan, 0., 1.]),
                                      ('scale', [0., -1., np.nan, np.inf, 1.])):
                        with self.subTest(checked=checked, family=family, method=method, parameter=name):
                            np.testing.assert_array_equal(func(.5, *shape, **{name:bad}),
                                [func(.5, *shape, **{name:v}) for v in bad])
                    for index in range(len(shape)):
                        bad = [0., -1., np.nan, np.inf, shape[index]]
                        args = list(shape)
                        args[index] = bad
                        expected = []
                        for v in bad:
                            scalar_args = list(shape)
                            scalar_args[index] = v
                            expected.append(func(.5, *scalar_args))
                        np.testing.assert_array_equal(func(.5, *args), expected)
                for name in ('loc', 'scale'):
                    a = np.array([1., 2., 3.])
                    with self.assertRaises(ValueError):
                        group.cdf(.5, *shape, out=a, **{name:a})
                    np.testing.assert_array_equal(a, [1., 2., 3.])
                for index in range(len(shape)):
                    a = np.array([1., 2., 3.])
                    args = list(shape)
                    args[index] = a
                    with self.assertRaises(ValueError):
                        group.cdf(.5, *args, out=a)
                    np.testing.assert_array_equal(a, [1., 2., 3.])
                with self.assertRaises(ValueError):
                    group.cdf(np.ones(2), *shape, scale=np.ones(3))
                with self.assertRaises(ValueError):
                    group.cdf(.5, *shape, scale=np.ones(3), out=np.empty((3,1)))
                with self.assertRaises(TypeError):
                    group.cdf(.5, *shape, scale=[1+2j])

    def test_infinite_student_df_is_exact_normal_limit(self):
        for _, api in self.apis:
            for method in METHODS:
                func, normal = getattr(api.t, method), getattr(api.normal, method)
                x = np.array([0., 1e-12, .1, .5, .9, 1.]) if method in ('ppf', 'isf') else np.array([-8., 0., 8.])
                np.testing.assert_array_equal(func(x, np.inf), normal(x))
                self.assertEqual(func(.5, np.inf), normal(.5))
                self.assertEqual(func([], np.inf).size, 0)
                df = np.full((2,1), np.inf)
                loc, scale = np.array([-.5, 0., .5]), np.array([1., 2., 3.])
                actual = func(x[:,None,None], df, loc=loc, scale=scale)
                expected = normal(x[:,None,None], loc=np.broadcast_to(loc, (2,3)), scale=scale)
                np.testing.assert_array_equal(actual, expected)
                self.assertEqual(func(np.empty((0,1)), [np.inf]*3).shape, (0,3))

    def test_generated_build_plan_and_stale_source_guard(self):
        self.assertEqual(set(FAMILIES), set(builder.PLAN['families']))
        self.assertIn('chi2', FAMILIES)
        for _, api in self.apis:
            for family in FAMILIES:
                specs = getattr(api, family)._specs
                self.assertEqual(set(specs), set(METHODS))
                self.assertTrue(all(s.get('elementwise') for s in specs.values()))
        self.assertEqual(builder.PLAN['modules']['scifort_chi2'], 'src/stats/chi2.f90')
        sources = list(builder.SOURCES)
        self.assertLess(sources.index('src/stats/gamma.f90'), sources.index('src/stats/chi2.f90'))
        hashes = dict(builder.PLAN['source_sha256'])
        hashes[sources[0]] = 'deliberately stale'
        with patch.dict(builder.PLAN, source_sha256=hashes):
            with self.assertRaisesRegex(RuntimeError, 'regenerate'):
                build()

    def test_unexpected_status_and_missing_library(self):
        with self.assertRaises(FileNotFoundError):
            Scifort(HERE / 'missing-library.so')
        for _, api in self.apis:
            group = api.normal
            previous = group._functions['cdf']
            def rejected(*args):
                C.cast(args[-1], C.POINTER(C.c_int))[0] = 17
            try:
                group._functions['cdf'] = (previous[0], rejected)
                with self.assertRaisesRegex(RuntimeError, '17'):
                    group.cdf([0., 1.])
            finally:
                group._functions['cdf'] = previous

    def test_round_trips_and_tail_probabilities(self):
        p = np.array([1e-12, 1e-6, .1, .5, .9, 1-1e-6])
        for _, api in self.apis:
            for family in FAMILIES:
                group = getattr(api, family)
                shape = SHAPES.get(family, ())
                np.testing.assert_allclose(group.cdf(group.ppf(p, *shape), *shape), p,
                                           rtol=2e-10, atol=0)
                np.testing.assert_allclose(group.sf(group.isf(p, *shape), *shape), p,
                                           rtol=2e-10, atol=0)
            z = np.array([-40., -12., -8., 8., 12., 40.])
            np.testing.assert_allclose(api.normal.logcdf(z), stats.norm.logcdf(z), rtol=4e-13, atol=0)
            np.testing.assert_allclose(api.normal.logsf(z), stats.norm.logsf(z), rtol=4e-13, atol=0)


if __name__ == '__main__':
    unittest.main()
