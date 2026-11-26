# SPDX-License-Identifier: MIT
"""Optional multivariate samplers preserve native streams and error contracts."""
import ctypes as C
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_library import build
from scifort import Scifort
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]


class VectorSamplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_optional_combinations_streams_and_normal_limit(self):
        mu=np.array([2.,-1.]);matrix=np.array([[4.,77.],[1.,3.]])
        for _,_,api in self.builds:
            for location in (None,mu):
                for covariance in (None,matrix):
                    for student in (False,True):
                        group=api.multivariate_t if student else api.multivariate_normal
                        kwargs=({'loc':location,'shape':covariance,'df':np.inf} if student else
                                {'mean':location,'cov':covariance})
                        state=api.rng(777);reference=api.rng(777)
                        r=group.rvs_array(state,n=2,n_samples=31,**kwargs)
                        self.assertEqual(r.status,0)
                        first=group.rvs_array(reference,n=2,n_samples=7,**kwargs)
                        rest=group.rvs_array(reference,n=2,n_samples=24,**kwargs)
                        np.testing.assert_array_equal(r.outputs['samples'],np.concatenate(
                            [first.outputs['samples'],rest.outputs['samples']],axis=1))
                        np.testing.assert_array_equal(state.state,reference.state)
                        consumed=api.rng(777);consumed.uniform(62)
                        np.testing.assert_array_equal(state.state,consumed.state)
                        normal_state=api.rng(777)
                        normal=api.multivariate_normal.rvs_array(normal_state,location,covariance,n=2,n_samples=31)
                        # Native t adds location after the zero-mean normal draw;
                        # normal batch accumulation includes location from the start.
                        np.testing.assert_allclose(r.outputs['samples'],normal.outputs['samples'],
                                                   rtol=2e-15,atol=2e-15)
                        np.testing.assert_array_equal(state.state,normal_state.state)
            self.assertEqual(api.multivariate_normal.rvs_array(api.rng()).outputs['samples'].shape,(1,1))
            state=api.rng(3);reference=api.rng(3)
            default=api.multivariate_t.rvs_array(state,n=2,n_samples=5)
            explicit=api.multivariate_t.rvs_array(reference,df=1,n=2,n_samples=5)
            np.testing.assert_array_equal(default.outputs['samples'],explicit.outputs['samples'])
            np.testing.assert_array_equal(state.state,reference.state)
            consumed=api.rng(3);consumed.uniform(15)
            np.testing.assert_array_equal(state.state,consumed.state)

    def test_sample_moments_singular_and_lower_triangle(self):
        mu=np.array([2.,-1.]);matrix=np.array([[4.,1.],[1.,3.]])
        for _,_,api in self.builds:
            for student in (False,True):
                state=api.rng(234);size=12000
                if student:r=api.multivariate_t.rvs_array(state,mu,matrix,df=10,n=2,n_samples=size)
                else:r=api.multivariate_normal.rvs_array(state,mu,matrix,n=2,n_samples=size)
                self.assertEqual(r.status,0);samples=r.outputs['samples']
                cov=matrix*(10/8 if student else 1)
                np.testing.assert_array_less(abs(samples.mean(axis=1)-mu),6*np.sqrt(np.diag(cov)/size))
                np.testing.assert_allclose(np.cov(samples,bias=True),cov,rtol=.12,atol=.08)
                kwargs={'loc':mu,'shape':np.diag([1.,0.]),'df':5} if student else {'mean':mu,'cov':np.diag([1.,0.])}
                group=api.multivariate_t if student else api.multivariate_normal
                singular=group.rvs_array(api.rng(4),n=2,n_samples=10,**kwargs)
                self.assertEqual(singular.status,0);np.testing.assert_array_equal(singular.outputs['samples'][1],mu[1])
                if student:kwargs['shape']=np.zeros((2,2))
                else:kwargs['cov']=np.zeros((2,2))
                zero=group.rvs_array(api.rng(4),n=2,n_samples=10,**kwargs)
                np.testing.assert_array_equal(zero.outputs['samples'],np.repeat(mu[:,None],10,axis=1))
                a=api.rng(5);b=api.rng(5)
                asymmetric=matrix.copy();asymmetric[0,1]=99
                if student:
                    x=group.rvs_array(a,mu,asymmetric,5,n=2,n_samples=10)
                    y=group.rvs_array(b,mu,matrix,5,n=2,n_samples=10)
                else:
                    x=group.rvs_array(a,mu,asymmetric,n=2,n_samples=10)
                    y=group.rvs_array(b,mu,matrix,n=2,n_samples=10)
                np.testing.assert_array_equal(x.outputs['samples'],y.outputs['samples'])

    def test_empty_batches_invalid_parameters_and_failure_advancement(self):
        for _,_,api in self.builds:
            state=api.rng(11);saved=state.state
            for group in (api.multivariate_normal,api.multivariate_t):
                self.assertEqual(group.rvs_array(state,n=2,n_samples=0).status,0)
                student=group is api.multivariate_t
                for bad in ([np.nan,0.],[np.inf,0.]):
                    kwargs={'loc':bad} if student else {'mean':bad}
                    result=group.rvs_array(state,n=2,n_samples=2,**kwargs)
                    self.assertEqual(result.status,2);self.assertTrue(np.all(np.isnan(result.outputs['samples'])))
                for bad,code in (([[1.,0.],[0.,-1.]],3),([[np.nan,0.],[0.,1.]],2)):
                    kwargs={'shape':bad} if student else {'cov':bad}
                    result=group.rvs_array(state,n=2,n_samples=2,**kwargs)
                    self.assertEqual(result.status,code);self.assertTrue(np.all(np.isnan(result.outputs['samples'])))
                    empty=group.rvs_array(state,n=2,n_samples=0,**kwargs)
                    self.assertEqual(empty.status,0 if student else code)
                np.testing.assert_array_equal(state.state,saved)
            for df in (0.,-1.,-np.inf,np.nan):
                self.assertEqual(api.multivariate_t.rvs_array(state,df=df,n=2,n_samples=0).status,0)
                result=api.multivariate_t.rvs_array(state,df=df,n=2,n_samples=2)
                self.assertEqual(result.status,2);self.assertTrue(np.all(np.isnan(result.outputs['samples'])))
                np.testing.assert_array_equal(state.state,saved)
            self.assertEqual(api.multivariate_t.rvs_array(state,n=0,n_samples=0).status,0)
            self.assertEqual(api.multivariate_t.rvs_array(state,n=0,n_samples=1).status,1)
            # Underflow occurs after drawing the normal vector and one chi-square uniform.
            expected=api.rng(11);expected.uniform(3)
            tiny=api.multivariate_t.rvs_array(state,df=1e-300,n=2,n_samples=2)
            self.assertEqual(tiny.status,2);self.assertTrue(np.all(np.isnan(tiny.outputs['samples'])))
            np.testing.assert_array_equal(state.state,expected.state)

    def test_shapes_types_ownership_and_caller_buffers(self):
        for _,_,api in self.builds:
            for group in (api.multivariate_normal,api.multivariate_t):
                state=api.rng(9);saved=state.state
                out=np.empty((2,3),order='F')
                result=group.rvs_array(state,n=2,n_samples=3,out={'samples':out})
                self.assertEqual(result.status,0);self.assertIs(result.outputs['samples'],out)
                saved=state.state
                for key in ('n','n_samples'):
                    for bad in (-1,2**31,True,1.5):
                        with self.assertRaises((TypeError,ValueError)):group.rvs_array(state,**{key:bad})
                for vector,matrix in (([1.],None),(None,np.eye(3)),([1.+1j,2.],None)):
                    kwargs=({'loc':vector,'shape':matrix} if group is api.multivariate_t else {'mean':vector,'cov':matrix})
                    with self.assertRaises((TypeError,ValueError)):group.rvs_array(state,n=2,**kwargs)
                readonly=np.empty((2,3),order='F');readonly.flags.writeable=False
                for bad in (np.empty((2,3)),np.empty((2,3),dtype=np.intc,order='F'),readonly,np.empty((3,2),order='F')):
                    with self.assertRaises((ValueError,TypeError)):group.rvs_array(state,n=2,n_samples=3,out={'samples':bad})
                vector=np.array([1.,2.]);kwargs={'loc':vector} if group is api.multivariate_t else {'mean':vector}
                with self.assertRaises(ValueError):group.rvs_array(state,n=2,out={'samples':vector.reshape(2,1)},**kwargs)
                with self.assertRaises(TypeError):group.rvs_array([1,2,3,4,5,6])
                np.testing.assert_array_equal(state.state,saved)
                # Explicit output dimension is authoritative, never inferred from parameters/out.
                with self.assertRaises(ValueError):group.rvs_array(state,vector)

    def test_raw_c_guard_preservation_and_python_agreement(self):
        ptr=lambda x:x.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            for group in (api.multivariate_normal,api.multivariate_t):
                student=group is api.multivariate_t
                rng=api.rng(8);state=rng.state.copy();out=np.full((2,3),77.,order='F')
                vector=np.array([2.,-1.]);matrix=np.array([[4.,1.],[1.,3.]],order='F');status=C.c_int(-1)
                raw=group._functions['rvs_array']
                def call(n=2,count=3,flag1=1,flag2=1):
                    args=[n,count,ptr(state),ptr(out),flag1,ptr(vector),flag2,ptr(matrix)]
                    if student:args.append(5.)
                    raw(*args,C.byref(status))
                call()
                result=(group.rvs_array(rng,vector,matrix,5,n=2,n_samples=3) if student else
                        group.rvs_array(rng,vector,matrix,n=2,n_samples=3))
                self.assertEqual(status.value,0);np.testing.assert_array_equal(out,result.outputs['samples'])
                np.testing.assert_array_equal(state,rng.state)
                saved=state.copy();out[:]=77
                for kwargs in ({'n':-1},{'count':-1},{'flag1':2},{'flag2':-1}):
                    call(**kwargs);self.assertEqual(status.value,1)
                    np.testing.assert_array_equal(out,77.);np.testing.assert_array_equal(state,saved)
                # Absent optional arrays ignore buffer values while pointers stay valid.
                vector[:]=np.nan;matrix[:]=np.nan
                call(flag1=0,flag2=0)
                self.assertEqual(status.value,0);self.assertTrue(np.all(np.isfinite(out)))
                out[:]=77;state[:]=0
                for count in (0,3):
                    call(count=count);self.assertEqual(status.value,4)
                    np.testing.assert_array_equal(out,77.);np.testing.assert_array_equal(state,0)

    def test_native_c_and_cpp_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_multivariate_sampling_c_api.f90',path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    client=folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 int64_t state[6]; double mu[2]={1,2},matrix[4]={0},out[4]; int status;
 scifort_rng_seed_i64(123,state,&status); if(status) return 1;
 scifort_multivariate_normal_rvs_array_f64(2,2,state,out,1,mu,1,matrix,&status);
 if(status || out[0]!=1 || out[1]!=2 || out[2]!=1 || out[3]!=2) return 2;
 scifort_multivariate_t_rvs_array_f64(2,2,state,out,1,mu,1,matrix,INFINITY,&status);
 return status || out[0]!=1 || out[1]!=2;
}
''',encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',client])
