# SPDX-License-Identifier: MIT
"""Directional density, fitting and explicit-state sampling interfaces."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import vonmises_fisher
from scipy.special import ive
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT = Path(__file__).resolve().parents[3]


class DirectionalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True, False):
            path, metadata = build(checked=checked)
            cls.builds.append((path, metadata, Scifort(path, dll_directories=metadata['dll_directories'])))

    def test_density_entropy_references_and_identities(self):
        for _, _, api in self.builds:
            for n in (2, 3, 5):
                mu = np.arange(1., n+1)
                mu /= np.linalg.norm(mu)
                for kappa in (.1, 3., 25.):
                    for x in (mu, -mu, np.eye(n)[0]):
                        x = x.copy()
                        for method in ('pdf', 'logpdf'):
                            value = getattr(api.vonmises_fisher, method)(x, mu, kappa)
                            expected = getattr(vonmises_fisher, method)(x, mu, kappa)
                            np.testing.assert_allclose(value, expected, rtol=5e-12, atol=3e-13)
                        self.assertAlmostEqual(np.exp(api.vonmises_fisher.logpdf(x, mu, kappa)),
                                               api.vonmises_fisher.pdf(x, mu, kappa), places=13)
                    np.testing.assert_allclose(api.vonmises_fisher.entropy(mu, kappa),
                                               vonmises_fisher.entropy(mu, kappa), rtol=5e-12, atol=3e-13)
                self.assertGreater(api.vonmises_fisher.pdf(mu.copy(), mu, 3.), api.vonmises_fisher.pdf(-mu, mu, 3.))
            kappa = 3.
            expected = np.log(kappa)-np.log(2*np.pi)-np.log1p(-np.exp(-2*kappa))
            self.assertAlmostEqual(api.vonmises_fisher.logpdf([1, 0, 0], [1, 0, 0], kappa), expected, places=13)

    def test_single_batch_split_streams_and_norms(self):
        for _, _, api in self.builds:
            for n in (2, 3, 5):
                mu = np.zeros(n); mu[-1] = 1.
                for vmf in (False, True):
                    group = api.vonmises_fisher if vmf else api.uniform_direction
                    args = (mu, 3.) if vmf else ()
                    sizes = {} if vmf else {'n': n}
                    state, repeated = api.rng(2**40+123), api.rng(2**40+123)
                    result = group.rvs_array(state, *args, n_samples=11, **sizes)
                    self.assertEqual(result.status, 0)
                    samples = result.outputs['samples']
                    singles = []
                    for _ in range(11):
                        draw = group.rvs(repeated, *args, **sizes)
                        self.assertEqual(draw.status, 0)
                        singles.append(draw.outputs['sample'])
                    np.testing.assert_array_equal(samples, np.stack(singles, axis=1))
                    np.testing.assert_array_equal(state.state, repeated.state)
                    np.testing.assert_allclose(np.linalg.norm(samples, axis=0), 1., rtol=0, atol=3e-14)
                    state, repeated = api.rng(12), api.rng(12)
                    whole = group.rvs_array(state, *args, n_samples=7, **sizes).outputs['samples']
                    a = group.rvs_array(repeated, *args, n_samples=2, **sizes).outputs['samples']
                    b = group.rvs_array(repeated, *args, n_samples=5, **sizes).outputs['samples']
                    np.testing.assert_array_equal(whole, np.concatenate([a, b], axis=1))
                    np.testing.assert_array_equal(state.state, repeated.state)
            single = api.uniform_direction.rvs(api.rng(), n=1)
            self.assertEqual(abs(single.outputs['sample'][0]), 1.)
            self.assertEqual(api.uniform_direction.rvs(api.rng()).outputs['sample'].shape, (3,))

    def test_uniform_moments_and_vmf_mean_resultant(self):
        for _, _, api in self.builds:
            uniform = api.uniform_direction.rvs_array(api.rng(11), n=3, n_samples=4000)
            self.assertEqual(uniform.status, 0)
            x = uniform.outputs['samples']
            self.assertLess(np.max(abs(x.mean(axis=1))), .05)
            np.testing.assert_allclose(np.mean(x*x, axis=1), 1/3, rtol=0, atol=.035)
            mu = np.array([.6, .8, 0.])
            result = api.vonmises_fisher.rvs_array(api.rng(13), mu, 4., n_samples=4000)
            self.assertEqual(result.status, 0)
            expected = (1/np.tanh(4.)-1/4)*mu
            np.testing.assert_allclose(result.outputs['samples'].mean(axis=1), expected, rtol=0, atol=.04)

    def test_fit_direction_concentration_and_degenerate_outputs(self):
        for _, _, api in self.builds:
            data = np.array([[1., .9, .8, .7, .6], [0., .1, -.2, .3, -.1], [0., .2, .3, -.1, .4]])
            data /= np.linalg.norm(data, axis=0)
            result = api.vonmises_fisher.fit(data)
            self.assertEqual(result.status, 0)
            mu, kappa = result.outputs['mu'], result.outputs['kappa']
            expected_mu, expected_kappa = vonmises_fisher.fit(data.T)
            np.testing.assert_allclose(mu, expected_mu, rtol=5e-14, atol=2e-15)
            np.testing.assert_allclose(kappa, expected_kappa, rtol=4e-10, atol=2e-11)
            r = np.linalg.norm(data.mean(axis=1))
            np.testing.assert_allclose(ive(1.5, kappa)/ive(.5, kappa), r, rtol=3e-12, atol=1e-13)
            for bad, code in ((np.ones((3, 4)), 2), (np.full((3, 4), np.nan), 1),
                              (np.array([[1., -1.], [0., 0.], [0., 0.]]), 3)):
                result = api.vonmises_fisher.fit(bad)
                self.assertEqual(result.status, code)
                self.assertTrue(np.isnan(result.outputs['kappa']))
                self.assertTrue(np.all(np.isnan(result.outputs['mu'])))
            result = api.vonmises_fisher.fit([[1.], [0.], [0.]])
            self.assertEqual(result.status, 3)
            np.testing.assert_array_equal(result.outputs['mu'], [1., 0., 0.])
            self.assertTrue(np.isnan(result.outputs['kappa']))

    def test_invalid_parameters_shapes_and_empty_batches(self):
        for _, _, api in self.builds:
            state = api.rng(12); saved = state.state
            mu = np.array([1., 0., 0.])
            for kappa in (0., -1., np.nan, np.inf, -np.inf):
                for method in ('pdf', 'logpdf'):
                    self.assertTrue(np.isnan(getattr(api.vonmises_fisher, method)(mu.copy(), mu, kappa)))
                self.assertTrue(np.isnan(api.vonmises_fisher.entropy(mu, kappa)))
                draw = api.vonmises_fisher.rvs_array(state, mu, kappa, n_samples=3)
                self.assertEqual(draw.status, 2)
                np.testing.assert_array_equal(draw.outputs['samples'], 0.)
            for bad in ([0., 0., 0.], [2., 0., 0.], [np.nan, 0., 0.], [np.inf, 0., 0.]):
                self.assertTrue(np.isnan(api.vonmises_fisher.pdf(bad, mu, 3.)))
                self.assertTrue(np.isnan(api.vonmises_fisher.pdf(mu, bad, 3.)))
                self.assertEqual(api.vonmises_fisher.rvs(state, bad, 3.).status, 2)
            empty = api.vonmises_fisher.rvs_array(state, [np.nan, 0., 0.], np.nan, n_samples=0)
            self.assertEqual(empty.status, 0)
            self.assertEqual(empty.outputs['samples'].shape, (3, 0))
            self.assertEqual(api.uniform_direction.rvs_array(state, n=3, n_samples=0).status, 0)
            np.testing.assert_array_equal(state.state, saved)
            for call in (lambda: api.vonmises_fisher.pdf([1., 0.], mu, 3.),
                         lambda: api.vonmises_fisher.rvs(state, [1.], 3.),
                         lambda: api.vonmises_fisher.fit(np.empty((3, 0))),
                         lambda: api.uniform_direction.rvs(state, n=0)):
                with self.assertRaises(ValueError): call()

    def test_buffers_dimensions_and_raw_c_guards(self):
        ptr = lambda a: a.ctypes.data_as(C.c_void_p)
        for _, _, api in self.builds:
            state = api.rng(12); saved = state.state
            mu = np.array([1., 0., 0.])
            for vmf in (False, True):
                group = api.vonmises_fisher if vmf else api.uniform_direction
                args = (mu, 3.) if vmf else ()
                sizes = {} if vmf else {'n': 3}
                for bad in (-1, True, 1.5, 2**31):
                    with self.assertRaises((TypeError, ValueError)):
                        group.rvs_array(state, *args, n_samples=bad, **sizes)
                readonly = np.empty((3, 4), order='F'); readonly.flags.writeable = False
                for out in (np.empty((3, 4)), np.empty((3, 4), dtype=np.intc, order='F'),
                            np.empty((4, 3), order='F'), readonly):
                    with self.assertRaises((TypeError, ValueError)):
                        group.rvs_array(state, *args, n_samples=4, out={'samples': out}, **sizes)
                out = np.full((3, 4), 77., order='F'); serialized = state.state.copy(); status = C.c_int(-1)
                raw = group._functions['rvs_array']
                def call(n=3, count=4):
                    values = [n, count, ptr(serialized)]
                    if vmf: values += [ptr(mu), 3.]
                    raw(*values, ptr(out), C.byref(status))
                call(n=-1)
                self.assertEqual(status.value, 1)
                np.testing.assert_array_equal(out, 77.); np.testing.assert_array_equal(serialized, saved)
                call(count=-1)
                self.assertEqual(status.value, 1)
                np.testing.assert_array_equal(out, 77.)
                serialized[:] = 0
                for count in (0, 4):
                    call(count=count)
                    self.assertEqual(status.value, 4)
                    np.testing.assert_array_equal(out, 77.); np.testing.assert_array_equal(serialized, 0)
            with self.assertRaises(ValueError):
                api.vonmises_fisher.rvs(state, mu, 3., out={'sample': mu})
            np.testing.assert_array_equal(state.state, saved)
            out = np.empty((3, 4), order='F')
            result = api.vonmises_fisher.rvs_array(state, mu, 3., n_samples=4, out={'samples': out})
            self.assertEqual(result.status, 0); self.assertIs(result.outputs['samples'], out)

    def test_native_c_and_cpp_clients(self):
        for path, metadata, _ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder = Path(temp); env = os.environ.copy()
                env['PATH'] = str(path.parent)+os.pathsep+env.get('PATH', '')
                flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe = folder/'native.exe'
                run([metadata['compiler_path'], *flags, '-I', path.parent, ROOT/'test/test_directional_c_api.f90', path, '-o', exe])
                run([exe], env=env)
                if shutil.which('gcc'):
                    client = folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 int64_t state[6]; double mu[3]={1,0,0},x[3]={1,0,0},out[12],value,fitted[3],kappa; int status;
 scifort_rng_seed_i64(123,state,&status); if(status) return 1;
 scifort_uniform_direction_rvs_f64(3,state,out,&status); if(status) return 2;
 scifort_uniform_direction_rvs_array_f64(3,4,state,out,&status); if(status) return 3;
 scifort_vonmises_fisher_pdf_f64(3,x,mu,3,&value,&status); if(status || value<=0) return 4;
 scifort_vonmises_fisher_logpdf_f64(3,x,mu,3,&value,&status); if(status || !isfinite(value)) return 5;
 scifort_vonmises_fisher_entropy_f64(3,mu,3,&value,&status); if(status || !isfinite(value)) return 6;
 scifort_vonmises_fisher_rvs_f64(3,state,mu,3,out,&status); if(status) return 7;
 scifort_vonmises_fisher_rvs_array_f64(3,4,state,mu,3,out,&status); if(status) return 8;
 scifort_vonmises_fisher_fit_f64(3,4,out,fitted,&kappa,&status);
 return status || !isfinite(kappa) || kappa<=0;
}
''', encoding='ascii')
                    exe = folder/'client.exe'
                    run(['gcc', '-std=c11', '-Wall', '-Wextra', '-Werror', '-I', ROOT/'include', client, path, '-o', exe])
                    run([exe], env=env)
                    if shutil.which('g++'):
                        run(['g++', '-std=c++17', '-Wall', '-Wextra', '-Werror', '-I', ROOT/'include', '-fsyntax-only', client])
