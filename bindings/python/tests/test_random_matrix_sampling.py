# SPDX-License-Identifier: MIT
"""Real Haar and prescribed-spectrum matrix bindings, streams and buffer contracts."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT = Path(__file__).resolve().parents[3]


class RandomMatrixSamplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True, False):
            path, metadata = build(checked=checked)
            cls.builds.append((path, metadata, Scifort(path, dll_directories=metadata['dll_directories'])))

    def test_haar_properties_single_batch_and_streams(self):
        for _, _, api in self.builds:
            for name in ('ortho_group', 'special_ortho_group'):
                group = getattr(api, name)
                state, repeated = api.rng(2**40+123), api.rng(2**40+123)
                result = group.rvs_array(state, n=4, n_samples=32)
                self.assertEqual(result.status, 0)
                samples = result.outputs['samples']
                singles = []
                for _ in range(32):
                    single = group.rvs(repeated, n=4)
                    self.assertEqual(single.status, 0)
                    singles.append(single.outputs['sample'])
                np.testing.assert_array_equal(samples, np.stack(singles, axis=2))
                np.testing.assert_array_equal(state.state, repeated.state)
                determinants = []
                for q in samples.transpose(2, 0, 1):
                    np.testing.assert_allclose(q.T @ q, np.eye(4), rtol=0, atol=2e-14)
                    determinants.append(np.linalg.det(q))
                if name == 'special_ortho_group':
                    np.testing.assert_allclose(determinants, 1., rtol=0, atol=2e-14)
                else:
                    self.assertTrue(any(d < 0 for d in determinants))
                    self.assertTrue(any(d > 0 for d in determinants))
                first_state, split_state = api.rng(12), api.rng(12)
                whole = group.rvs_array(first_state, n=3, n_samples=7).outputs['samples']
                a = group.rvs_array(split_state, n=3, n_samples=2).outputs['samples']
                b = group.rvs_array(split_state, n=3, n_samples=5).outputs['samples']
                np.testing.assert_array_equal(whole, np.concatenate([a, b], axis=2))
                np.testing.assert_array_equal(first_state.state, split_state.state)
                self.assertEqual(group.rvs(api.rng()).outputs['sample'].shape, (1, 1))
                if name == 'special_ortho_group':
                    np.testing.assert_array_equal(group.rvs(api.rng(), n=1).outputs['sample'], [[1.]])

    def test_haar_moments_and_empty_extents(self):
        for _, _, api in self.builds:
            for group in (api.ortho_group, api.special_ortho_group):
                samples = group.rvs_array(api.rng(13), n=3, n_samples=3000).outputs['samples']
                self.assertLess(np.max(abs(samples.mean(axis=2))), .06)
                np.testing.assert_allclose(np.mean(samples**2, axis=2), 1/3, rtol=0, atol=.035)
                state = api.rng(12)
                saved = state.state
                for n, count in ((0, 0), (0, 3), (3, 0)):
                    result = group.rvs_array(state, n=n, n_samples=count)
                    self.assertEqual(result.status, 0)
                    self.assertEqual(result.outputs['samples'].shape, (n, n, count))
                self.assertEqual(group.rvs(state, n=0).status, 0)
                np.testing.assert_array_equal(state.state, saved)

    def test_correlation_spectrum_defaults_and_stream(self):
        for _, _, api in self.builds:
            for eigs in ([.2, .8, 2.], [0., 0., 3.], [1., 1., 1.], [0., 2.]):
                state, repeated = api.rng(12), api.rng(12)
                result = api.random_correlation.rvs(state, eigs)
                explicit = api.random_correlation.rvs(repeated, np.asarray(eigs)[::-1][::-1],
                                                       tol=1e-13, diag_tol=1e-7)
                self.assertEqual(result.status, 0)
                matrix = result.outputs['sample']
                np.testing.assert_array_equal(matrix, explicit.outputs['sample'])
                np.testing.assert_array_equal(state.state, repeated.state)
                np.testing.assert_array_equal(np.diag(matrix), np.ones(len(eigs)))
                np.testing.assert_allclose(matrix, matrix.T, rtol=0, atol=2e-15)
                np.testing.assert_allclose(np.linalg.eigvalsh(matrix), np.sort(eigs), rtol=0, atol=2e-13)
                advanced = api.rng(12)
                api.ortho_group.rvs(advanced, n=len(eigs))
                np.testing.assert_array_equal(state.state, advanced.state)

    def test_correlation_invalid_parameters_and_shape_guards(self):
        for _, _, api in self.builds:
            state = api.rng(12)
            saved = state.state
            for eigs in ([np.nan, 2.], [np.inf, 0.], [-1., 3.], [1., 2.]):
                out = np.full((2, 2), 77., order='F')
                result = api.random_correlation.rvs(state, eigs, out={'sample': out})
                self.assertEqual(result.status, 2)
                np.testing.assert_array_equal(out, 0.)
            for key in ('tol', 'diag_tol'):
                for bad in (-1., np.nan, np.inf, -np.inf):
                    result = api.random_correlation.rvs(state, [1., 1.], **{key: bad})
                    self.assertEqual(result.status, 2)
                    np.testing.assert_array_equal(result.outputs['sample'], 0.)
            for bad in ([], [1.], [[1., 1.]], [1+1j, 1]):
                with self.assertRaises((TypeError, ValueError)):
                    api.random_correlation.rvs(state, bad)
            np.testing.assert_array_equal(state.state, saved)

    def test_python_buffers_dimensions_and_ownership(self):
        for _, _, api in self.builds:
            state = api.rng(12)
            saved = state.state
            for group in (api.ortho_group, api.special_ortho_group):
                for key in ('n', 'n_samples'):
                    for bad in (-1, True, 1.5, 2**31):
                        with self.assertRaises((TypeError, ValueError)):
                            group.rvs_array(state, **{key: bad})
                readonly = np.empty((2, 2, 3), order='F')
                readonly.flags.writeable = False
                for bad in (np.empty((2, 2, 3)), np.empty((2, 2, 3), dtype=np.intc, order='F'),
                            np.empty((3, 2, 3), order='F'), readonly):
                    with self.assertRaises((TypeError, ValueError)):
                        group.rvs_array(state, n=2, n_samples=3, out={'samples': bad})
                with self.assertRaises(TypeError):
                    group.rvs([1, 2, 3, 4, 5, 6])
            # A valid Fortran buffer sharing storage with an input must be rejected.
            storage = np.ones(4)
            with self.assertRaises(ValueError):
                api.random_correlation.rvs(state, storage[:2], out={'sample': storage.reshape((2, 2), order='F')})
            np.testing.assert_array_equal(state.state, saved)
            out = np.empty((2, 2, 3), order='F')
            result = api.ortho_group.rvs_array(state, n=2, n_samples=3, out={'samples': out})
            self.assertEqual(result.status, 0)
            self.assertIs(result.outputs['samples'], out)

    def test_raw_c_guards_state_and_python_agreement(self):
        ptr = lambda a: a.ctypes.data_as(C.c_void_p)
        for _, _, api in self.builds:
            for name in ('ortho_group', 'special_ortho_group', 'random_correlation'):
                group = getattr(api, name)
                rng = api.rng(12)
                state = rng.state.copy()
                out = np.full((3, 3), 77., order='F')
                eigs = np.array([.2, .8, 2.])
                status = C.c_int(-1)
                raw = group._functions['rvs']
                def call(n=3):
                    args = [n, ptr(state)]
                    if name == 'random_correlation':
                        args += [ptr(eigs), ptr(out), C.byref(status), 1e-13, 1e-7]
                    else:
                        args += [ptr(out), C.byref(status)]
                    raw(*args)
                call()
                expected = group.rvs(rng, eigs) if name == 'random_correlation' else group.rvs(rng, n=3)
                self.assertEqual(status.value, 0)
                np.testing.assert_array_equal(out, expected.outputs['sample'])
                np.testing.assert_array_equal(state, rng.state)
                saved = state.copy()
                out[:] = 77.
                for n in ((-1, 1) if name == 'random_correlation' else (-1,)):
                    call(n)
                    self.assertEqual(status.value, 1)
                    np.testing.assert_array_equal(state, saved)
                    np.testing.assert_array_equal(out, 77.)
                state[:] = 0
                call()
                self.assertEqual(status.value, 4)
                np.testing.assert_array_equal(state, 0)
                np.testing.assert_array_equal(out, 77.)

    def test_raw_c_batch_guards_and_empty_state_validation(self):
        ptr = lambda a: a.ctypes.data_as(C.c_void_p)
        for _, _, api in self.builds:
            for group in (api.ortho_group, api.special_ortho_group):
                state = api.rng(12).state.copy()
                saved = state.copy()
                out = np.full((3, 3, 2), 77., order='F')
                status = C.c_int(-1)
                raw = group._functions['rvs_array']
                for n, count in ((-1, 2), (3, -1)):
                    raw(n, count, ptr(state), ptr(out), C.byref(status))
                    self.assertEqual(status.value, 1)
                    np.testing.assert_array_equal(state, saved)
                    np.testing.assert_array_equal(out, 77.)
                state[:] = 0
                for count in (0, 2):
                    raw(3, count, ptr(state), ptr(out), C.byref(status))
                    self.assertEqual(status.value, 4)
                    np.testing.assert_array_equal(state, 0)
                    np.testing.assert_array_equal(out, 77.)

    def test_native_c_and_cpp_clients(self):
        for path, metadata, _ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder = Path(temp)
                env = os.environ.copy()
                env['PATH'] = str(path.parent)+os.pathsep+env.get('PATH', '')
                flags = ['-O0', '-g', '-fcheck=all', '-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe = folder/'native.exe'
                run([metadata['compiler_path'], *flags, '-I', path.parent,
                     ROOT/'test/test_random_matrix_c_api.f90', path, '-o', exe])
                run([exe], env=env)
                if shutil.which('gcc'):
                    client = folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 int64_t state[6]; double out[18], eigs[3]={.2,.8,2}; int status;
 scifort_rng_seed_i64(123,state,&status); if(status) return 1;
 scifort_ortho_group_rvs_f64(3,state,out,&status); if(status) return 2;
 scifort_ortho_group_rvs_array_f64(3,2,state,out,&status); if(status) return 3;
 scifort_special_ortho_group_rvs_f64(3,state,out,&status); if(status) return 4;
 scifort_special_ortho_group_rvs_array_f64(3,2,state,out,&status); if(status) return 5;
 scifort_random_correlation_rvs_f64(3,state,eigs,out,&status,1e-13,1e-7);
 return status || !isfinite(out[0]) || out[0]!=1.0;
}
''', encoding='ascii')
                    exe = folder/'client.exe'
                    run(['gcc', '-std=c11', '-Wall', '-Wextra', '-Werror', '-I', ROOT/'include', client, path, '-o', exe])
                    run([exe], env=env)
                    if shutil.which('g++'):
                        run(['g++', '-std=c++17', '-Wall', '-Wextra', '-Werror', '-I', ROOT/'include', '-fsyntax-only', client])
