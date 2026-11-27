# SPDX-License-Identifier: MIT
"""Fixed-margin table probabilities, integer cubes, streams and overflow boundaries."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import random_table
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT = Path(__file__).resolve().parents[3]


class RandomTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True, False):
            path, metadata = build(checked=checked)
            cls.builds.append((path, metadata, Scifort(path, dll_directories=metadata['dll_directories'])))

    def test_probabilities_means_and_normalization(self):
        tables = [np.array([[k, 2-k], [2-k, k]]) for k in range(3)]
        for _, _, api in self.builds:
            probabilities = [api.random_table.pmf(x, [2, 2], [2, 2]) for x in tables]
            np.testing.assert_allclose(probabilities, [1/6, 4/6, 1/6], rtol=2e-14)
            self.assertAlmostEqual(sum(probabilities), 1., places=14)
            for x, row, col in ((tables[1], [2, 2], [2, 2]),
                                (np.array([[1, 1, 1], [1, 3, 1]]), [3, 5], [2, 4, 2]),
                                (np.array([[1, 3, 2]]), [6], [1, 3, 2])):
                for method in ('pmf', 'logpmf'):
                    np.testing.assert_allclose(getattr(api.random_table, method)(x, row, col),
                                               getattr(random_table, method)(x, row, col), rtol=3e-13, atol=2e-14)
                np.testing.assert_allclose(api.random_table.mean(row, col), random_table.mean(row, col),
                                           rtol=2e-15, atol=0)
            x = np.zeros((2, 3), dtype=np.intc)
            self.assertEqual(api.random_table.pmf(x, [0, 0], [0, 0, 0]), 1.)
            self.assertEqual(api.random_table.logpmf(x, [0, 0], [0, 0, 0]), 0.)
            self.assertTrue(np.all(np.isnan(api.random_table.mean([0, 0], [0, 0, 0]))))
            self.assertEqual(api.random_table.pmf([[0, 0], [0, 0]], [2, 2], [2, 2]), 0.)
            self.assertEqual(api.random_table.logpmf([[0, 0], [0, 0]], [2, 2], [2, 2]), -np.inf)

    def test_sampler_streams_and_fixed_margins(self):
        for _, _, api in self.builds:
            row, col = [3, 5], [2, 4, 2]
            state, repeated = api.rng(2**40+123), api.rng(2**40+123)
            result = api.random_table.rvs_array(state, row, col, n_samples=17)
            self.assertEqual(result.status, 0)
            samples = result.outputs['samples']
            self.assertEqual(samples.dtype, np.dtype(np.intc))
            self.assertTrue(samples.flags.f_contiguous)
            singles = []
            for _ in range(17):
                draw = api.random_table.rvs(repeated, row, col)
                self.assertEqual(draw.status, 0)
                singles.append(draw.outputs['sample'])
            np.testing.assert_array_equal(samples, np.stack(singles, axis=2))
            np.testing.assert_array_equal(state.state, repeated.state)
            np.testing.assert_array_equal(samples.sum(axis=0), np.repeat(np.array(col)[:, None], 17, axis=1))
            np.testing.assert_array_equal(samples.sum(axis=1), np.repeat(np.array(row)[:, None], 17, axis=1))
            self.assertTrue(np.all(samples >= 0))
            state, split = api.rng(12), api.rng(12)
            whole = api.random_table.rvs_array(state, row, col, n_samples=7).outputs['samples']
            first = api.random_table.rvs_array(split, row, col, n_samples=2).outputs['samples']
            rest = api.random_table.rvs_array(split, row, col, n_samples=5).outputs['samples']
            np.testing.assert_array_equal(whole, np.concatenate([first, rest], axis=2))
            np.testing.assert_array_equal(state.state, split.state)

    def test_sampler_law_and_degenerate_tables(self):
        for _, _, api in self.builds:
            result = api.random_table.rvs_array(api.rng(13), [2, 2], [2, 2], n_samples=6000)
            self.assertEqual(result.status, 0)
            freq = np.bincount(result.outputs['samples'][0, 0, :], minlength=3)/6000
            p = np.array([1/6, 4/6, 1/6])
            self.assertTrue(np.all(abs(freq-p) < 6*np.sqrt(p*(1-p)/6000)))
            for row, col in (([0, 0], [0, 0, 0]), ([6], [1, 3, 2]), ([2, 4], [6])):
                draw = api.random_table.rvs(api.rng(12), row, col)
                self.assertEqual(draw.status, 0)
                np.testing.assert_array_equal(draw.outputs['sample'].sum(axis=1), row)
                np.testing.assert_array_equal(draw.outputs['sample'].sum(axis=0), col)

    def test_invalid_counts_overflow_boundaries_and_empty_batches(self):
        maximum = np.iinfo(np.intc).max
        for _, _, api in self.builds:
            state = api.rng(12)
            saved = state.state
            for row, col in (([-1, 3], [1, 1]), ([1, 2], [1, 1]),
                             ([maximum, 1], [0]), ([maximum, maximum], [maximum, maximum]),
                             ([maximum], [maximum])):
                x = np.zeros((len(row), len(col)), dtype=np.intc)
                for method in ('pmf', 'logpmf'):
                    self.assertTrue(np.isnan(getattr(api.random_table, method)(x, row, col)))
                self.assertTrue(np.all(np.isnan(api.random_table.mean(row, col))))
                result = api.random_table.rvs_array(state, row, col, n_samples=2)
                self.assertEqual(result.status, 2)
                np.testing.assert_array_equal(result.outputs['samples'], 0)
            self.assertEqual(api.random_table.pmf([[maximum-1]], [maximum-1], [maximum-1]), 1.)
            np.testing.assert_allclose(api.random_table.mean([maximum-1], [maximum-1]), [[maximum-1.]], rtol=2e-15)
            # Off-support table rows/columns may exceed the integer range even for small margins.
            self.assertEqual(api.random_table.pmf(np.full((2, 2), maximum, dtype=np.intc), [1, 1], [1, 1]), 0.)
            self.assertTrue(np.isnan(api.random_table.pmf([[-1, 1], [1, 1]], [1, 1], [1, 1])))
            for row, col in (([], []), ([], [0]), ([0], []), ([-1], [3])):
                result = api.random_table.rvs_array(state, np.array(row, dtype=np.intc),
                                                     np.array(col, dtype=np.intc), n_samples=0)
                self.assertEqual(result.status, 0)
                self.assertEqual(result.outputs['samples'].shape, (len(row), len(col), 0))
            np.testing.assert_array_equal(state.state, saved)
            for method, args in (('mean', (np.array([], dtype=np.intc), np.array([], dtype=np.intc))),
                                 ('rvs', (state, np.array([], dtype=np.intc), np.array([], dtype=np.intc))),
                                 ('pmf', (np.empty((0, 0), dtype=np.intc), np.array([], dtype=np.intc),
                                          np.array([], dtype=np.intc)))):
                with self.assertRaises(ValueError):
                    getattr(api.random_table, method)(*args)

    def test_python_input_output_contracts(self):
        for _, _, api in self.builds:
            state = api.rng(12)
            saved = state.state
            for bad in (-1, True, 1.5, 2**31):
                with self.assertRaises((ValueError, TypeError)):
                    api.random_table.rvs_array(state, [2, 2], [2, 2], n_samples=bad)
            for bad in ([1.5, 2.5], [1+1j, 3], [2**31, 0], [True, False]):
                with self.assertRaises((ValueError, TypeError)):
                    api.random_table.mean(bad, [2, 2])
            with self.assertRaises(ValueError):
                api.random_table.pmf(np.zeros((3, 2), dtype=np.intc), [2, 2], [2, 2])
            readonly = np.empty((2, 2, 3), dtype=np.intc, order='F')
            readonly.flags.writeable = False
            for out in (np.empty((2, 2, 3), order='F'), np.empty((2, 2, 3), dtype=np.intc),
                        np.empty((3, 2, 3), dtype=np.intc, order='F'), readonly):
                with self.assertRaises((ValueError, TypeError)):
                    api.random_table.rvs_array(state, [2, 2], [2, 2], n_samples=3, out={'samples': out})
            storage = np.array([2, 2, 0, 0], dtype=np.intc)
            with self.assertRaises(ValueError):
                api.random_table.rvs(state, storage[:2], [2, 2], out={'sample': storage.reshape((2, 2), order='F')})
            np.testing.assert_array_equal(state.state, saved)
            out = np.empty((2, 2, 3), dtype=np.intc, order='F')
            result = api.random_table.rvs_array(state, [2, 2], [2, 2], n_samples=3, out={'samples': out})
            self.assertEqual(result.status, 0)
            self.assertIs(result.outputs['samples'], out)

    def test_raw_c_guards_and_state(self):
        ptr = lambda a: a.ctypes.data_as(C.c_void_p)
        for _, _, api in self.builds:
            group = api.random_table
            state = api.rng(12).state.copy()
            saved = state.copy()
            row = np.array([2, 2], dtype=np.intc)
            col = row.copy()
            out = np.full((2, 2, 3), 77, dtype=np.intc, order='F')
            status = C.c_int(-1)
            raw = group._functions['rvs_array']
            for nr, nc, count in ((-1, 2, 3), (2, -1, 3), (2, 2, -1)):
                raw(nr, nc, count, ptr(state), ptr(row), ptr(col), ptr(out), C.byref(status))
                self.assertEqual(status.value, 1)
                np.testing.assert_array_equal(state, saved)
                np.testing.assert_array_equal(out, 77)
            state[:] = 0
            for count in (0, 3):
                raw(2, 2, count, ptr(state), ptr(row), ptr(col), ptr(out), C.byref(status))
                self.assertEqual(status.value, 4)
                np.testing.assert_array_equal(state, 0)
                np.testing.assert_array_equal(out, 77)
            value = C.c_double(77.)
            x = np.ones((2, 2), dtype=np.intc, order='F')
            group._functions['pmf'](-1, 2, ptr(x), ptr(row), ptr(col), C.byref(value), C.byref(status))
            self.assertEqual(status.value, 1)
            self.assertEqual(value.value, 77.)

    def test_native_c_and_cpp_clients(self):
        for path, metadata, _ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder = Path(temp)
                env = os.environ.copy()
                env['PATH'] = str(path.parent)+os.pathsep+env.get('PATH', '')
                flags = ['-O0', '-g', '-fcheck=all', '-ftrapv'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe = folder/'native.exe'
                run([metadata['compiler_path'], *flags, '-J', folder, '-I', path.parent,
                     ROOT/'src/stats/random_table.f90', ROOT/'test/test_random_table_c_api.f90', path, '-o', exe], cwd=folder)
                run([exe], env=env)
                if shutil.which('gcc'):
                    client = folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <limits.h>
#include <math.h>
int main(void) {
 int64_t state[6]; int row[2]={2,2},col[2]={2,2},x[4]={1,1,1,1},samples[12],status;
 double p,lp,mu[4];
 scifort_random_table_pmf_f64(2,2,x,row,col,&p,&status); if(status || fabs(p-2.0/3)>1e-14) return 1;
 scifort_random_table_logpmf_f64(2,2,x,row,col,&lp,&status); if(status || fabs(exp(lp)-p)>1e-14) return 2;
 scifort_random_table_mean_f64(2,2,row,col,mu,&status); if(status || mu[0]!=1) return 3;
 scifort_rng_seed_i64(123,state,&status); if(status) return 4;
 scifort_random_table_rvs_f64(2,2,state,row,col,x,&status); if(status) return 5;
 scifort_random_table_rvs_array_f64(2,2,3,state,row,col,samples,&status); if(status) return 6;
 row[0]=INT_MAX; row[1]=1;
 scifort_random_table_pmf_f64(2,2,x,row,col,&p,&status);
 return status || !isnan(p);
}
''', encoding='ascii')
                    exe = folder/'client.exe'
                    run(['gcc', '-std=c11', '-Wall', '-Wextra', '-Werror', '-I', ROOT/'include', client, path, '-o', exe])
                    run([exe], env=env)
                    if shutil.which('g++'):
                        run(['g++', '-std=c++17', '-Wall', '-Wextra', '-Werror', '-I', ROOT/'include', '-fsyntax-only', client])
