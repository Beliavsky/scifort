# SPDX-License-Identifier: MIT
"""Optional fixed parameters, column-oriented fits and native/C/Python agreement."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import multivariate_normal
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]


class MultivariateFitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_all_four_constraints_scipy_and_independent_references(self):
        data=np.array([[1.,3.,5.,4.],[2.,6.,4.,7.]])
        for _,_,api in self.builds:
            for fixed_mean in (None,np.array([2.,3.])):
                for fixed_cov in (None,np.array([[2.,.2],[.2,3.]])):
                    result=api.multivariate_normal.fit(data,fixed_mean,fixed_cov)
                    self.assertEqual(result.status,0)
                    expected=multivariate_normal.fit(data.T,fix_mean=fixed_mean,fix_cov=fixed_cov)
                    mean=data.mean(axis=1) if fixed_mean is None else fixed_mean
                    centered=data-mean[:,None]
                    cov=centered@centered.T/data.shape[1] if fixed_cov is None else fixed_cov
                    np.testing.assert_allclose(result.outputs['mean'],expected[0],rtol=2e-14,atol=1e-15)
                    np.testing.assert_allclose(result.outputs['cov'],expected[1],rtol=2e-14,atol=1e-15)
                    np.testing.assert_allclose(result.outputs['mean'],mean,rtol=2e-14,atol=1e-15)
                    np.testing.assert_allclose(result.outputs['cov'],cov,rtol=2e-14,atol=1e-15)

    def test_single_observation_singular_fixed_covariance_and_lower_triangle(self):
        for _,_,api in self.builds:
            fit=api.multivariate_normal.fit
            result=fit([[2.],[5.]])
            self.assertEqual(result.status,0)
            np.testing.assert_array_equal(result.outputs['mean'],[2.,5.])
            np.testing.assert_array_equal(result.outputs['cov'],np.zeros((2,2)))
            result=fit([[2.],[5.]],fixed_mean=[1.,3.])
            np.testing.assert_array_equal(result.outputs['cov'],[[1.,2.],[2.,4.]])
            for matrix in (np.zeros((2,2)),np.diag([2.,0.]),np.array([[2.,99.],[.2,3.]])):
                result=fit([[1.,3.],[2.,6.]],fixed_cov=matrix)
                self.assertEqual(result.status,0)
                expected=np.tril(matrix)+np.tril(matrix,-1).T
                np.testing.assert_array_equal(result.outputs['cov'],expected)
            result=fit([[1.,2.,3.]])
            np.testing.assert_allclose(result.outputs['cov'],[[2/3]],rtol=2e-14)

    def test_numerical_errors_preserve_native_partial_outputs(self):
        data=np.array([[1.,3.,5.],[2.,6.,4.]])
        for _,_,api in self.builds:
            fit=api.multivariate_normal.fit
            for bad in (np.nan,np.inf,-np.inf):
                invalid=data.copy();invalid[0,0]=bad
                result=fit(invalid)
                self.assertEqual(result.status,1)
                self.assertTrue(np.all(np.isnan(result.outputs['mean'])))
                self.assertTrue(np.all(np.isnan(result.outputs['cov'])))
                result=fit(data,fixed_mean=[bad,1.])
                self.assertEqual(result.status,2)
                self.assertTrue(np.all(np.isnan(result.outputs['mean'])))
                self.assertTrue(np.all(np.isnan(result.outputs['cov'])))
            for fixed_cov in (np.diag([1.,-1.]),np.full((2,2),np.nan),np.full((2,2),np.inf)):
                for fixed_mean in (None,[9.,8.]):
                    mean=np.full(2,77.);cov=np.full((2,2),88.,order='F')
                    result=fit(data,fixed_mean,fixed_cov,out={'mean':mean,'cov':cov})
                    self.assertEqual(result.status,2)
                    np.testing.assert_array_equal(mean,data.mean(axis=1) if fixed_mean is None else fixed_mean)
                    self.assertTrue(np.all(np.isnan(cov)))

    def test_shapes_layouts_out_and_alias_validation(self):
        data=np.array([[1.,3.,5.],[2.,6.,4.]])
        for _,_,api in self.builds:
            fit=api.multivariate_normal.fit
            for item in (data,data.astype('>f8'),data[:,::-1]):
                item.flags.writeable=False
                mean=np.empty(2);cov=np.empty((2,2),order='F')
                result=fit(item,out={'mean':mean,'cov':cov})
                self.assertEqual(result.status,0)
                self.assertIs(result.outputs['mean'],mean);self.assertIs(result.outputs['cov'],cov)
                np.testing.assert_allclose(mean,[3.,4.],rtol=2e-14)
            for args in (([],),([1.,2.],),(np.empty((2,0)),),(np.empty((0,2)),),
                         (data,[1.]),(data,None,np.eye(3)),(data,[],None),
                         (data,None,[1.,2.]),(np.ones((2,2,1)),)):
                with self.assertRaises(ValueError):fit(*args)
            with self.assertRaises(TypeError):fit(data.astype(complex))
            with self.assertRaises(TypeError):fit(data,fixed_mean=[1+1j,2.])
            with self.assertRaises(ValueError):fit(data,fixed_mean=data[:,0])
            with self.assertRaises(ValueError):fit(data,out={'mean':data[:,0]})
            with self.assertRaises(ValueError):fit(data,out={'cov':np.empty((2,2),order='C')})
            with self.assertRaises(TypeError):fit(data,out={'mean':np.empty(2,dtype=np.float32)})
            shared=np.empty((2,2),order='F')
            with self.assertRaises(ValueError):fit(data,out={'cov':shared,'mean':shared[:,0]})

    def test_raw_c_flags_extents_and_ignored_placeholders(self):
        data=np.array([[1.,3.,5.],[2.,6.,4.]],order='F')
        fixed_mean=np.full(2,np.nan);fixed_cov=np.full((2,2),np.nan,order='F')
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            raw=api.multivariate_normal._functions['fit'];mean=np.empty(2);cov=np.empty((2,2),order='F');status=C.c_int(-1)
            raw(2,3,ptr(data),ptr(mean),ptr(cov),C.byref(status),0,ptr(fixed_mean),0,ptr(fixed_cov))
            self.assertEqual(status.value,0);np.testing.assert_array_equal(mean,[3.,4.])
            for n,k,m,c in ((0,3,0,0),(2,0,0,0),(-1,3,0,0),(2,3,-1,0),(2,3,1,2),(2,3,2,1)):
                mean[:]=77.;cov[:]=88.
                raw(n,k,ptr(data),ptr(mean),ptr(cov),C.byref(status),m,ptr(fixed_mean),c,ptr(fixed_cov))
                self.assertEqual(status.value,1)
                np.testing.assert_array_equal(mean,[77.,77.]);np.testing.assert_array_equal(cov,np.full((2,2),88.))

    def test_native_and_c_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_multivariate_fit_c_api.f90',path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    client=folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 double x[6]={1,2,3,6,5,4},m[2],c[4],fm[2]={2,3},fc[4]={2,.2,.2,3}; int status;
 scifort_multivariate_normal_fit_f64(2,3,x,m,c,&status,0,fm,0,fc);
 if(status || m[0]!=3 || m[1]!=4 || fabs(c[0]-8./3)>1e-14) return 1;
 scifort_multivariate_normal_fit_f64(2,3,x,m,c,&status,1,fm,0,fc);
 if(status || m[0]!=2 || fabs(c[0]-11./3)>1e-14) return 2;
 scifort_multivariate_normal_fit_f64(2,3,x,m,c,&status,0,fm,1,fc);
 if(status || m[0]!=3 || c[0]!=2 || c[3]!=3) return 3;
 scifort_multivariate_normal_fit_f64(2,3,x,m,c,&status,1,fm,1,fc);
 return status || m[0]!=2 || m[1]!=3 || c[0]!=2 || c[3]!=3;
}
''',encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',client])
