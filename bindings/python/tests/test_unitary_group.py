# SPDX-License-Identifier: MIT
"""Unitary sampler complex transport, Haar properties and explicit stream contracts."""
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


class UnitaryGroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'): raise unittest.SkipTest('gfortran required')
        cls.builds = []
        for checked in (True, False):
            path, metadata = build(checked=checked)
            cls.builds.append((path, metadata, Scifort(path, dll_directories=metadata['dll_directories'])))

    def test_unitarity_and_identical_single_batch_streams(self):
        for _, _, api in self.builds:
            for n in (1, 2, 4):
                state, repeated = api.rng(2**40+123), api.rng(2**40+123)
                result = api.unitary_group.rvs_array(state, n=n, n_samples=16)
                self.assertEqual(result.status, 0)
                samples = result.outputs['samples']
                self.assertEqual(samples.dtype, np.dtype(np.complex128))
                self.assertTrue(samples.flags.f_contiguous)
                singles = []
                for _ in range(16):
                    draw = api.unitary_group.rvs(repeated, n=n)
                    self.assertEqual(draw.status, 0)
                    singles.append(draw.outputs['sample'])
                np.testing.assert_array_equal(samples, np.stack(singles, axis=2))
                np.testing.assert_array_equal(state.state, repeated.state)
                self.assertTrue(np.any(samples.imag != 0))
                for q in samples.transpose(2, 0, 1):
                    np.testing.assert_allclose(q.conj().T @ q, np.eye(n), rtol=0, atol=2e-14)
                    self.assertAlmostEqual(abs(np.linalg.det(q)), 1., places=13)
                state, split = api.rng(12), api.rng(12)
                whole = api.unitary_group.rvs_array(state, n=n, n_samples=7).outputs['samples']
                first = api.unitary_group.rvs_array(split, n=n, n_samples=2).outputs['samples']
                rest = api.unitary_group.rvs_array(split, n=n, n_samples=5).outputs['samples']
                np.testing.assert_array_equal(whole, np.concatenate([first, rest], axis=2))
                np.testing.assert_array_equal(state.state, split.state)
            self.assertEqual(api.unitary_group.rvs(api.rng()).outputs['sample'].shape, (1, 1))

    def test_haar_moments_and_empty_extents(self):
        for _, _, api in self.builds:
            result = api.unitary_group.rvs_array(api.rng(13), n=3, n_samples=3000)
            self.assertEqual(result.status, 0)
            x = result.outputs['samples']
            self.assertLess(np.max(abs(x.mean(axis=2))), .07)
            np.testing.assert_allclose(np.mean(abs(x)**2, axis=2), 1/3, rtol=0, atol=.035)
            np.testing.assert_allclose(np.mean(x.real**2, axis=2), 1/6, rtol=0, atol=.025)
            np.testing.assert_allclose(np.mean(x.imag**2, axis=2), 1/6, rtol=0, atol=.025)
            state = api.rng(12); saved = state.state
            for n, count in ((0, 0), (0, 3), (3, 0)):
                result = api.unitary_group.rvs_array(state, n=n, n_samples=count)
                self.assertEqual(result.status, 0)
                self.assertEqual(result.outputs['samples'].shape, (n, n, count))
                self.assertEqual(result.outputs['samples'].dtype, np.dtype(np.complex128))
            self.assertEqual(api.unitary_group.rvs(state, n=0).status, 0)
            np.testing.assert_array_equal(state.state, saved)

    def test_buffers_dimensions_and_ownership(self):
        for _, _, api in self.builds:
            state = api.rng(12); saved = state.state
            for key in ('n', 'n_samples'):
                for bad in (-1, True, 1.5, 2**31):
                    with self.assertRaises((TypeError, ValueError)):
                        api.unitary_group.rvs_array(state, **{key: bad})
            readonly = np.empty((2, 2, 3), dtype=np.complex128, order='F')
            readonly.flags.writeable = False
            unaligned = np.ndarray((2, 2, 3), dtype=np.complex128,
                                   buffer=bytearray(2*2*3*16+1), offset=1, order='F')
            for bad in (np.empty((2, 2, 3), order='F'), np.empty((2, 2, 3), dtype=np.complex64, order='F'),
                        np.empty((2, 2, 3), dtype=np.complex128),
                        np.empty((2, 2, 3), dtype='>c16', order='F'),
                        np.empty((3, 2, 3), dtype=np.complex128, order='F'), readonly, unaligned):
                with self.assertRaises((TypeError, ValueError)):
                    api.unitary_group.rvs_array(state, n=2, n_samples=3, out={'samples': bad})
            with self.assertRaises(TypeError): api.unitary_group.rvs([1, 2, 3, 4, 5, 6])
            np.testing.assert_array_equal(state.state, saved)
            out = np.empty((2, 2, 3), dtype=np.complex128, order='F')
            result = api.unitary_group.rvs_array(state, n=2, n_samples=3, out={'samples': out})
            self.assertEqual(result.status, 0); self.assertIs(result.outputs['samples'], out)

    def test_raw_interleaved_c_transport_and_guards(self):
        ptr = lambda a: a.ctypes.data_as(C.c_void_p)
        for _, _, api in self.builds:
            rng = api.rng(12); serialized = rng.state.copy()
            # Physical Fortran shape (2,n,n,n_samples), filled through a plain double pointer.
            raw_out = np.full((2, 3, 3, 2), 77., order='F')
            raw = api.unitary_group._functions['rvs_array']; status = C.c_int(-1)
            raw(3, 2, ptr(serialized), ptr(raw_out), C.byref(status))
            self.assertEqual(status.value, 0)
            expected = api.unitary_group.rvs_array(rng, n=3, n_samples=2)
            np.testing.assert_array_equal(raw_out[0]+1j*raw_out[1], expected.outputs['samples'])
            np.testing.assert_array_equal(serialized, rng.state)
            saved = serialized.copy(); raw_out[:] = 77.
            for n, count in ((-1, 2), (3, -1)):
                raw(n, count, ptr(serialized), ptr(raw_out), C.byref(status))
                self.assertEqual(status.value, 1)
                np.testing.assert_array_equal(raw_out, 77.); np.testing.assert_array_equal(serialized, saved)
            serialized[:] = 0
            for n, count in ((0, 0), (3, 0), (3, 2)):
                raw(n, count, ptr(serialized), ptr(raw_out), C.byref(status))
                self.assertEqual(status.value, 4)
                np.testing.assert_array_equal(raw_out, 77.); np.testing.assert_array_equal(serialized, 0)

    def test_native_c_and_cpp_clients(self):
        for path, metadata, _ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder = Path(temp); env = os.environ.copy()
                env['PATH'] = str(path.parent)+os.pathsep+env.get('PATH', '')
                flags = ['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe = folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_unitary_group_c_api.f90',path,'-o',exe])
                run([exe],env=env)
                client = folder/'client.c'
                client.write_text('''#include "scifort.h"
#include <math.h>
static int unitary(const double *q) {
 for(int a=0;a<3;++a) for(int b=0;b<3;++b) {
  double re=0,im=0;
  for(int r=0;r<3;++r) {
   double ar=q[2*(r+3*a)],ai=q[2*(r+3*a)+1],br=q[2*(r+3*b)],bi=q[2*(r+3*b)+1];
   re+=ar*br+ai*bi; im+=ar*bi-ai*br;
  }
  if(fabs(re-(a==b))>1e-12 || fabs(im)>1e-12) return 0;
 }
 return 1;
}
int main(void) {
 int64_t state[6]; double out[36]; int status;
 scifort_rng_seed_i64(123,state,&status); if(status) return 1;
 scifort_unitary_group_rvs_f64(3,state,out,&status); if(status || !unitary(out)) return 2;
 scifort_unitary_group_rvs_array_f64(3,2,state,out,&status);
 return status || !unitary(out) || !unitary(out+18);
}
''',encoding='ascii')
                for compiler, standard in (('gcc','c11'), ('g++','c++17')):
                    if not shutil.which(compiler): continue
                    exe = folder/(compiler+'.exe')
                    run([compiler,'-std='+standard,'-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe])
                    run([exe],env=env)
