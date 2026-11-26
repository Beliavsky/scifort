# SPDX-License-Identifier: MIT
"""Entropy defaults, marginal index conventions and native/C parameter agreement."""
import ctypes as C
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import multivariate_normal, multivariate_t, norm, t
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]


class MultivariateParameterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_entropy_references_defaults_and_normal_limit(self):
        matrix=np.array([[4.,1.,.5],[1.,3.,.2],[.5,.2,2.]])
        for _,_,api in self.builds:
            normal=api.multivariate_normal;student=api.multivariate_t
            self.assertIsInstance(normal.entropy(),float)
            np.testing.assert_allclose(normal.entropy(),norm.entropy(),rtol=2e-14)
            np.testing.assert_allclose(student.entropy(),t.entropy(1),rtol=3e-13)
            for mean,cov in (([2.],[[3.]]),([1.,2.,3.],matrix)):
                np.testing.assert_allclose(normal.entropy(mean,cov),multivariate_normal.entropy(mean,cov),rtol=3e-13)
                np.testing.assert_allclose(normal.entropy(cov=cov),normal.entropy(mean,cov),rtol=2e-14)
                np.testing.assert_allclose(normal.entropy(mean=mean),len(mean)*norm.entropy(),rtol=2e-14)
                for df in (1.,5.,100.,1e6):
                    expected=multivariate_t.entropy(loc=mean,shape=cov,df=df)
                    np.testing.assert_allclose(student.entropy(mean,cov,df),expected,rtol=4e-12,atol=2e-13)
                    np.testing.assert_allclose(student.entropy(shape=cov,df=df),expected,rtol=4e-12,atol=2e-13)
                np.testing.assert_allclose(student.entropy(mean,cov,np.inf),normal.entropy(mean,cov),rtol=2e-14)
            np.testing.assert_allclose(student.entropy(loc=[1.,2.]),
                                       multivariate_t.entropy(loc=[1.,2.],shape=np.eye(2),df=1),rtol=3e-13)

    def test_entropy_singular_native_convention_and_invalid_inputs(self):
        for _,_,api in self.builds:
            normal=api.multivariate_normal;student=api.multivariate_t
            singular=np.diag([2.,0.])
            np.testing.assert_allclose(normal.entropy(cov=singular),norm.entropy(scale=math.sqrt(2)),rtol=2e-14)
            self.assertEqual(normal.entropy(cov=np.zeros((2,2))),0.)
            # Native finite-df t entropy uses ambient dimension and pseudo-determinant.
            np.testing.assert_allclose(student.entropy(shape=singular,df=5),
                                       student.entropy(shape=np.diag([2.,1.]),df=5),rtol=2e-14)
            np.testing.assert_allclose(student.entropy(shape=singular,df=np.inf),normal.entropy(cov=singular),rtol=2e-14)
            for group,center,matrix in ((normal,'mean','cov'),(student,'loc','shape')):
                for kwargs in ({matrix:np.diag([1.,-1.])},{matrix:[[np.nan]]},
                               {matrix:[[np.inf]]},{center:[np.nan]}):
                    self.assertTrue(math.isnan(group.entropy(**kwargs)))
                for kwargs in ({center:[]},{matrix:np.empty((0,0))},{matrix:np.ones((2,3))},
                               {center:[0.],matrix:np.eye(2)},{matrix:[1.,2.]}):
                    with self.assertRaises(ValueError):group.entropy(**kwargs)
                with self.assertRaises(TypeError):group.entropy(**{matrix:[[1+1j]]})
            for df in (0.,-1.,np.nan,-np.inf):
                self.assertTrue(math.isnan(student.entropy(df=df)))

    def test_marginal_order_lower_triangle_and_density_consistency(self):
        center=np.array([1.,2.,3.]);matrix=np.array([[4.,1.,.5],[1.,3.,.2],[.5,.2,2.]])
        for _,_,api in self.builds:
            for group,cname,mname in ((api.multivariate_normal,'mean','cov'),(api.multivariate_t,'loc','shape')):
                for indices in ([2,0],[1],[0,1,2]):
                    result=group.marginal(indices,center,matrix)
                    self.assertEqual(result.status,0)
                    expected_center=center[indices];expected_matrix=matrix[np.ix_(indices,indices)]
                    np.testing.assert_array_equal(result.outputs['marginal_'+cname],expected_center)
                    np.testing.assert_array_equal(result.outputs['marginal_'+mname],expected_matrix)
                    kwargs={cname:result.outputs['marginal_'+cname],mname:result.outputs['marginal_'+mname]}
                    expected=(multivariate_normal.pdf(expected_center,expected_center,expected_matrix)
                              if cname=='mean' else multivariate_t.pdf(expected_center,expected_center,expected_matrix))
                    np.testing.assert_allclose(group.pdf(expected_center,**kwargs),expected,rtol=3e-13)
                upper=matrix.copy();upper[np.triu_indices(3,1)]=[91.,92.,93.]
                result=group.marginal([2,0],center,upper)
                np.testing.assert_array_equal(result.outputs['marginal_'+mname],matrix[np.ix_([2,0],[2,0])])
                # Marginal selection validates finite entries, not positive definiteness.
                result=group.marginal([0],center,np.diag([-1.,1.,1.]))
                self.assertEqual(result.status,0)
                self.assertEqual(result.outputs['marginal_'+mname][0,0],-1.)

    def test_marginal_conversions_buffers_and_failures(self):
        center=np.array([1.,2.,3.]);matrix=np.eye(3)
        for _,_,api in self.builds:
            for group,cname,mname in ((api.multivariate_normal,'mean','cov'),(api.multivariate_t,'loc','shape')):
                keys=('marginal_'+cname,'marginal_'+mname)
                for indices in (np.array([2,0],dtype='>i8'),np.array([2,9,0])[::2],np.array([2,0],dtype=np.uint64)):
                    indices.flags.writeable=False
                    out={keys[0]:np.empty(2),keys[1]:np.empty((2,2),order='F')}
                    result=group.marginal(indices,center,matrix,out=out)
                    self.assertEqual(result.status,0)
                    self.assertIs(result.outputs[keys[0]],out[keys[0]])
                    self.assertIs(result.outputs[keys[1]],out[keys[1]])
                    np.testing.assert_array_equal(indices,[2,0])
                for indices in ([-1,0],[3,0],[2**40,0],np.array([2**64-1,0],dtype=np.uint64),
                                [0.,1.],[True,False],[[0,1]],np.array([],dtype=int)):
                    out={keys[0]:np.array([77.,78.]),keys[1]:np.full((2,2),88.,order='F')}
                    with self.assertRaises((ValueError,TypeError)):
                        group.marginal(indices,center,matrix,out=out)
                    np.testing.assert_array_equal(out[keys[0]],[77.,78.])
                out={keys[0]:np.full(2,77.),keys[1]:np.full((2,2),88.,order='F')}
                for indices,point,spread in (([0,0],center,matrix),([0,1],[np.nan,2.,3.],matrix),
                                             ([0,1],center,np.full((3,3),np.inf))):
                    result=group.marginal(indices,point,spread,out=out)
                    self.assertEqual(result.status,2)
                    np.testing.assert_array_equal(out[keys[0]],[77.,77.])
                    np.testing.assert_array_equal(out[keys[1]],np.full((2,2),88.))
                with self.assertRaises(ValueError):group.marginal([0,1],center,matrix,out={keys[0]:center[:2]})
                with self.assertRaises(ValueError):group.marginal([0,1],center,matrix,out={keys[1]:np.empty((2,2),order='C')})
                with self.assertRaises(TypeError):group.marginal([0,1],center,matrix,out={keys[0]:np.empty(2,dtype=np.float32)})

    def test_raw_c_statuses_indices_and_output_preservation(self):
        center=np.array([1.,2.,3.]);matrix=np.eye(3,order='F');indices=np.array([3,1],dtype=np.intc)
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            for group in (api.multivariate_normal,api.multivariate_t):
                indices[:]=[3,1]
                raw=group._functions['marginal'];a=np.full(2,77.);b=np.full((2,2),88.,order='F');status=C.c_int(-1)
                raw(3,2,ptr(indices),ptr(center),ptr(matrix),ptr(a),ptr(b),C.byref(status))
                self.assertEqual(status.value,0);np.testing.assert_array_equal(a,[3.,1.])
                for n,m,selection,expected in ((0,2,[3,1],1),(3,0,[3,1],1),(3,2,[0,1],2),
                                              (3,2,[4,1],2),(3,2,[1,1],2),(1,2,[1,1],2)):
                    indices[:]=selection;a[:]=77.;b[:]=88.
                    raw(n,m,ptr(indices),ptr(center),ptr(matrix),ptr(a),ptr(b),C.byref(status))
                    self.assertEqual(status.value,expected)
                    np.testing.assert_array_equal(a,[77.,77.]);np.testing.assert_array_equal(b,np.full((2,2),88.))
                raw=group._functions['entropy'];value=C.c_double(99.)
                args=[0,ptr(center),ptr(matrix)]+([5.] if group is api.multivariate_t else [])
                raw(*args,C.byref(value),C.byref(status))
                self.assertEqual(status.value,1);self.assertEqual(value.value,99.)

    def test_native_and_c_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_multivariate_parameters_c_api.f90',path,'-o',exe])
                run([exe],env=env)
                if shutil.which('gcc'):
                    client=folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 double mean[3]={1,2,3},matrix[9]={1,0,0,0,1,0,0,0,1},value,a[2],b[4];
 int indices[2]={3,1},status;
 scifort_multivariate_normal_entropy_f64(3,mean,matrix,&value,&status);
 if(status || fabs(value-4.256815599614018)>1e-13) return 1;
 scifort_multivariate_t_entropy_f64(3,mean,matrix,5,&value,&status);
 if(status || !isfinite(value)) return 2;
 scifort_multivariate_normal_marginal_f64(3,2,indices,mean,matrix,a,b,&status);
 if(status || a[0]!=3 || a[1]!=1 || b[0]!=1 || b[1]!=0 || b[3]!=1) return 3;
 scifort_multivariate_t_marginal_f64(3,2,indices,mean,matrix,a,b,&status);
 return status || a[0]!=3 || a[1]!=1 || b[0]!=1 || b[1]!=0 || b[3]!=1;
}
''',encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',client])
