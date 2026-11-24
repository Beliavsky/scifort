# SPDX-License-Identifier: MIT
"""Validate all eight generated array methods with checked and optimized libraries."""
from pathlib import Path
import shutil
import sys
import unittest

import numpy as np
from scipy import linalg

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort, ArrayResult


class LinalgPythonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):
            raise unittest.SkipTest('gfortran required')
        cls.apis = []
        for checked in (True,False):
            path,metadata = build(checked=checked)
            cls.apis.append((checked,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_all_eight_methods_and_reference_values(self):
        a = np.array([[4.,1.,2.],[1.,3.,.5],[2.,.5,5.]])
        for checked,api in self.apis:
            group = api.linalg
            self.assertEqual(len(group._specs),8)
            result = group.cholesky_lower(a)
            self.assertIsInstance(result,ArrayResult)
            self.assertEqual(result.status,0)
            factor = result.outputs['l']
            self.assertTrue(factor.flags.f_contiguous)
            np.testing.assert_allclose(factor,linalg.cholesky(a,lower=True),rtol=3e-15,atol=0)
            result = group.spd_factor(a)
            self.assertEqual(result.status,0)
            np.testing.assert_array_equal(result.outputs['l'],factor)
            self.assertAlmostEqual(result.outputs['logdet'],np.linalg.slogdet(a)[1],places=14)
            eigen = group.symmetric_eigen_jacobi(a)
            self.assertEqual(eigen.status,0)
            w,v = eigen.outputs['eigenvalues'],eigen.outputs['eigenvectors']
            np.testing.assert_allclose(w,linalg.eigh(a,eigvals_only=True),rtol=3e-14,atol=4e-14)
            np.testing.assert_allclose(v.T@v,np.eye(3),rtol=0,atol=4e-14)
            np.testing.assert_allclose(v@np.diag(w)@v.T,a,rtol=3e-14,atol=4e-14)
            psd = group.psd_decompose(a)
            self.assertEqual(psd.status,0)
            self.assertEqual(psd.outputs['rank'],3)
            self.assertIsInstance(psd.outputs['rank'],int)
            self.assertAlmostEqual(psd.outputs['log_pdet'],np.linalg.slogdet(a)[1],places=13)
            inverse = group.invert_lower(factor)
            self.assertEqual(inverse.status,0)
            np.testing.assert_allclose(inverse.outputs['linv'],
                linalg.solve_triangular(factor,np.eye(3),lower=True),rtol=4e-15,atol=0)
            b = np.array([2.,-1.,4.])
            solved = group.solve_lower(factor,b)
            self.assertEqual(solved.status,0)
            np.testing.assert_allclose(solved.outputs['x'],linalg.solve_triangular(factor,b,lower=True),
                                       rtol=4e-15,atol=0)
            rhs = np.column_stack([b,b+1])
            for name,reference in (('solve_lower_matrix',linalg.solve_triangular(factor,rhs,lower=True)),
                                   ('solve_spd_matrix',linalg.cho_solve((factor,True),rhs))):
                solved = getattr(group,name)(factor,rhs)
                self.assertEqual(solved.status,0)
                self.assertEqual(solved.outputs['x'].shape,(3,2))
                self.assertTrue(solved.outputs['x'].flags.f_contiguous)
                np.testing.assert_allclose(solved.outputs['x'],reference,rtol=4e-15,atol=0)

    def test_statuses_partial_outputs_and_empty_batches(self):
        for _,api in self.apis:
            group = api.linalg
            singular = np.diag([4.,0.,1.])
            allowed = group.psd_decompose(singular)
            rejected = group.psd_decompose(singular,allow_singular=False)
            self.assertEqual(allowed.status,0)
            self.assertEqual(rejected.status,5)
            self.assertEqual(allowed.outputs['rank'],2)
            self.assertEqual(rejected.outputs['rank'],2)
            self.assertAlmostEqual(allowed.outputs['log_pdet'],np.log(4.),places=14)
            self.assertEqual(rejected.outputs['log_pdet'],0.)
            np.testing.assert_array_equal(allowed.outputs['eigenvalues'],rejected.outputs['eigenvalues'])
            self.assertGreater(rejected.outputs['cutoff'],0.)
            indefinite = group.psd_decompose(np.diag([-1.,2.]))
            self.assertEqual(indefinite.status,3)
            self.assertEqual(indefinite.outputs['rank'],0)
            bad = np.diag([np.nan,1.])
            self.assertEqual(group.cholesky_lower(bad).status,2)
            self.assertEqual(group.symmetric_eigen_jacobi(bad).status,4)
            psd = group.psd_decompose(bad)
            self.assertEqual(psd.status,4)
            self.assertEqual([psd.outputs[k] for k in ('rank','log_pdet','cutoff')],[0,0.,0.])
            spd = group.spd_factor(bad)
            self.assertEqual(spd.status,2)
            np.testing.assert_array_equal(spd.outputs['l'],np.zeros((2,2)))
            self.assertEqual(spd.outputs['logdet'],0.)
            # Complete/partial native failure output is returned, including caller-owned buffers.
            out = np.full((3,3),99.,order='F')
            result = group.cholesky_lower(singular,out={'l':out})
            self.assertEqual(result.status,2)
            self.assertIs(result.outputs['l'],out)
            self.assertEqual(out[0,0],2.)
            self.assertEqual(out[1,1],0.)
            self.assertEqual(group.solve_lower(singular,np.ones(3)).status,5)
            self.assertEqual(group.invert_lower(singular).status,3)
            for name in ('solve_lower_matrix','solve_spd_matrix'):
                output = np.empty((2,0),order='F')
                result = getattr(group,name)(bad,np.empty((2,0)),out={'x':output})
                self.assertEqual(result.status,0)
                self.assertIs(result.outputs['x'],output)
            # Cholesky/eigenpairs ignore the unused upper triangle, but spd_factor validates it.
            upper_bad = np.array([[4.,np.nan],[1.,3.]])
            self.assertEqual(group.cholesky_lower(upper_bad).status,0)
            self.assertEqual(group.symmetric_eigen_jacobi(upper_bad).status,0)
            self.assertEqual(group.spd_factor(upper_bad).status,2)

    def test_input_conversion_and_native_failure_agreement(self):
        base = np.array([[4.,1.],[1.,3.]])
        views = [base.astype('f4'),base.astype('>f8'),base.astype('i4'),base.tolist(),
                 np.kron(base,np.ones((2,2)))[::2,::2],np.asfortranarray(base)]
        unaligned = np.ndarray((2,2),dtype='f8',buffer=bytearray(33),offset=1)
        unaligned[:] = base
        views.append(unaligned)
        readonly = base.copy()
        readonly.flags.writeable = False
        views.append(readonly)
        for _,api in self.apis:
            for value in views:
                result = api.linalg.cholesky_lower(value)
                self.assertEqual(result.status,0)
                np.testing.assert_allclose(result.outputs['l'],linalg.cholesky(base,lower=True),
                                           rtol=3e-15,atol=0)
                np.testing.assert_array_equal(value,base)
            for value in ([0.],np.array([[-1.]]),np.array([[np.inf]])):
                matrix = np.asarray(value).reshape(1,1)
                self.assertNotEqual(api.linalg.cholesky_lower(matrix).status,0)

    def test_shape_and_type_validation_before_writes(self):
        for _,api in self.apis:
            group = api.linalg
            for bad in (np.ones(3),np.ones((3,2)),np.empty((0,0)),np.ones((1,1,1))):
                with self.assertRaises(ValueError):
                    group.cholesky_lower(bad)
            for bad in ([[1+2j]], [['1']], [[object()]]):
                with self.assertRaises(TypeError):
                    group.cholesky_lower(bad)
            for flag in (0,1,2,None,[],np.array(True)):
                with self.assertRaises(TypeError):
                    group.psd_decompose(np.eye(2),allow_singular=flag)
            self.assertEqual(group.psd_decompose(np.eye(2),allow_singular=np.bool_(True)).status,0)
            with self.assertRaises(ValueError):
                group.solve_lower(np.eye(2),np.ones(3))
            with self.assertRaises(ValueError):
                group.solve_lower_matrix(np.eye(2),np.ones(2))
            huge = np.broadcast_to(np.array(1.),(np.iinfo(np.intc).max+1,1))
            with self.assertRaises(ValueError):
                group.cholesky_lower(huge)
            output = np.full((2,2),77.,order='F')
            for supplied in ({'unknown':output},{'logdet':np.empty(())},output):
                with self.assertRaises((TypeError,ValueError)):
                    group.spd_factor(np.eye(2),out=supplied)
            np.testing.assert_array_equal(output,77.)

    def test_output_layout_and_overlap(self):
        for _,api in self.apis:
            group = api.linalg
            original = np.eye(2,order='F')
            with self.assertRaises(ValueError):
                group.cholesky_lower(original,out={'l':original})
            np.testing.assert_array_equal(original,np.eye(2))
            with self.assertRaises(ValueError):
                group.solve_lower_matrix(original,original)
            bad = [np.empty((2,2),order='C'),np.empty((2,2),dtype='f4',order='F'),
                   np.empty((2,2),dtype='>f8',order='F'),np.empty((3,2),order='F'),
                   np.ndarray((2,2),dtype='f8',buffer=bytearray(33),offset=1,order='F')]
            readonly = np.empty((2,2),order='F')
            readonly.flags.writeable = False
            bad.append(readonly)
            for out in bad:
                with self.assertRaises((TypeError,ValueError)):
                    group.cholesky_lower(np.eye(2),out={'l':out})
            storage = np.full(4,88.)
            with self.assertRaises(ValueError):
                group.symmetric_eigen_jacobi(np.eye(2),out={
                    'eigenvalues':storage[:2], 'eigenvectors':storage.reshape((2,2),order='F')})
            np.testing.assert_array_equal(storage,88.)
            factor = np.empty((2,2),order='F')
            result = group.cholesky_lower(np.eye(2),out={'l':factor})
            self.assertEqual(result.status,0)
            self.assertIs(result.outputs['l'],factor)
