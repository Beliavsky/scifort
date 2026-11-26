# SPDX-License-Identifier: MIT
"""Matrix samplers: optional parameters, streams, moments and defined failures."""
import ctypes as C
import itertools
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


class MatrixSamplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_all_optional_combinations_streams_and_defaults(self):
        mean=np.arange(6.).reshape(2,3);row=np.array([[2.,.3],[.3,1.]])
        col=np.array([[1.,.1,.2],[.1,2.,.3],[.2,.3,3.]])
        for _,_,api in self.builds:
            for student in (False,True):
                group=api.matrix_t if student else api.matrix_normal
                for flags in itertools.product((False,True),repeat=3):
                    args=[value if flag else None for value,flag in zip((mean,row,col),flags)]
                    kwargs={'df':5} if student else {}
                    state=api.rng(2**40+123);reference=api.rng(2**40+123)
                    result=group.rvs_array(state,*args,m=2,n=3,n_samples=9,**kwargs)
                    self.assertEqual(result.status,0);samples=result.outputs['samples']
                    first=group.rvs_array(reference,*args,m=2,n=3,n_samples=4,**kwargs)
                    rest=group.rvs_array(reference,*args,m=2,n=3,n_samples=5,**kwargs)
                    np.testing.assert_array_equal(samples,np.concatenate([first.outputs['samples'],rest.outputs['samples']],axis=2))
                    np.testing.assert_array_equal(state.state,reference.state)
                    consumed=api.rng(2**40+123);consumed.uniform((9 if student else 6)*9)
                    np.testing.assert_array_equal(state.state,consumed.state)
                    explicit=[value if flag else default for value,flag,default in
                              zip((mean,row,col),flags,(np.zeros((2,3)),np.eye(2),np.eye(3)))]
                    same=group.rvs_array(api.rng(2**40+123),*explicit,m=2,n=3,n_samples=9,**kwargs)
                    np.testing.assert_array_equal(samples,same.outputs['samples'])
                self.assertEqual(group.rvs_array(api.rng()).outputs['samples'].shape,(1,1,1))
            a=api.rng(12);b=api.rng(12)
            default=api.matrix_t.rvs_array(a,m=2,n=3,n_samples=3)
            explicit=api.matrix_t.rvs_array(b,df=1,m=2,n=3,n_samples=3)
            np.testing.assert_array_equal(default.outputs['samples'],explicit.outputs['samples'])
            np.testing.assert_array_equal(a.state,b.state)

    def test_moments_row_column_covariance_and_lower_triangles(self):
        mean=np.arange(6.).reshape(2,3);row=np.array([[2.,.3],[.3,1.]])
        col=np.array([[1.,.1,.2],[.1,2.,.3],[.2,.3,3.]])
        for _,_,api in self.builds:
            for student in (False,True):
                group=api.matrix_t if student else api.matrix_normal
                kwargs={'df':12} if student else {};size=14000
                result=group.rvs_array(api.rng(237),mean,row,col,m=2,n=3,n_samples=size,**kwargs)
                self.assertEqual(result.status,0);samples=result.outputs['samples']
                covariance=np.kron(col,row)/(10 if student else 1)
                variance=np.diag(covariance).reshape((2,3),order='F')
                np.testing.assert_array_less(abs(samples.mean(axis=2)-mean),6*np.sqrt(variance/size))
                actual=np.cov(samples.reshape((6,size),order='F'),bias=True)
                # Scale covariance errors by marginal SDs, including small cross terms.
                normalized=(actual-covariance)/np.sqrt(np.outer(np.diag(covariance),np.diag(covariance)))
                self.assertLess(np.max(abs(normalized)),.08)
                asym_row=row.copy();asym_row[0,1]=99
                asym_col=col.copy();asym_col[np.triu_indices(3,1)]=99
                first=group.rvs_array(api.rng(11),mean,row,col,m=2,n=3,n_samples=3,**kwargs)
                second=group.rvs_array(api.rng(11),mean,asym_row,asym_col,m=2,n=3,n_samples=3,**kwargs)
                np.testing.assert_array_equal(first.outputs['samples'],second.outputs['samples'])

    def test_empty_invalid_parameters_and_failure_state(self):
        for _,_,api in self.builds:
            for student in (False,True):
                group=api.matrix_t if student else api.matrix_normal
                state=api.rng(11);saved=state.state
                kwargs={'df':np.nan} if student else {}
                empty=group.rvs_array(state,np.full((2,3),np.nan),np.zeros((2,2)),np.zeros((3,3)),
                                      m=2,n=3,n_samples=0,**kwargs)
                self.assertEqual(empty.status,0)
                for m,n in ((0,3),(2,0),(0,0)):
                    self.assertEqual(group.rvs_array(state,m=m,n=n,n_samples=0).status,0)
                    self.assertEqual(group.rvs_array(state,m=m,n=n).status,1)
                # Nonfinite mean keeps the existing family-specific native status.
                result=group.rvs_array(state,np.full((2,3),np.nan),m=2,n=3,n_samples=3)
                self.assertEqual(result.status,2 if student else 1)
                self.assertTrue(np.all(np.isnan(result.outputs['samples'])))
                row_name='row_spread' if student else 'rowcov'
                col_name='col_spread' if student else 'colcov'
                for bad,code in ((np.full((2,2),np.nan),2),(np.zeros((2,2)),3)):
                    out=np.full((2,3,3),77.,order='F')
                    result=group.rvs_array(state,m=2,n=3,n_samples=3,out={'samples':out},**{row_name:bad})
                    self.assertEqual(result.status,code);self.assertTrue(np.all(np.isnan(out)))
                np.testing.assert_array_equal(state.state,saved)
                # Matrix-t validates column SPD after the inverse-Wishart draw.
                result=group.rvs_array(state,m=2,n=3,n_samples=3,**{col_name:np.zeros((3,3))})
                self.assertEqual(result.status,3);self.assertTrue(np.all(np.isnan(result.outputs['samples'])))
                reference=api.rng(11)
                if student:reference.uniform(3)
                np.testing.assert_array_equal(state.state,reference.state)
            for df in (0,-1,np.inf,-np.inf,np.nan):
                state=api.rng(9);saved=state.state
                self.assertEqual(api.matrix_t.rvs_array(state,df=df,n_samples=0).status,0)
                result=api.matrix_t.rvs_array(state,df=df,n_samples=3)
                self.assertEqual(result.status,2);self.assertTrue(np.all(np.isnan(result.outputs['samples'])))
                np.testing.assert_array_equal(state.state,saved)
            state=api.rng(11);reference=api.rng(11);reference.uniform(1)
            tiny=api.matrix_t.rvs_array(state,df=1e-300,n_samples=3)
            self.assertEqual(tiny.status,3);self.assertTrue(np.all(np.isnan(tiny.outputs['samples'])))
            np.testing.assert_array_equal(state.state,reference.state)

            # A fixed later mixture failure preserves a completed first matrix.
            state=api.rng(2**40+123);state.uniform(10)
            reference=api.rng(2**40+123);reference.uniform(10)
            first=api.matrix_t.rvs_array(reference,df=.0001)
            self.assertEqual(first.status,0)
            self.assertTrue(np.all(np.isfinite(first.outputs['samples'])))
            reference.uniform(1)
            out=np.full((1,1,3),77.,order='F')
            result=api.matrix_t.rvs_array(state,df=.0001,n_samples=3,out={'samples':out})
            self.assertEqual(result.status,3)
            np.testing.assert_array_equal(out[:,:,0],first.outputs['samples'][:,:,0])
            self.assertTrue(np.all(np.isnan(out[:,:,1:])))
            np.testing.assert_array_equal(state.state,reference.state)

    def test_shape_counts_buffers_and_ownership(self):
        for _,_,api in self.builds:
            for group in (api.matrix_normal,api.matrix_t):
                state=api.rng(9);saved=state.state
                for name in ('m','n','n_samples'):
                    for bad in (-1,True,1.2,2**31):
                        with self.assertRaises((ValueError,TypeError)):group.rvs_array(state,**{name:bad})
                for args in ((np.ones((3,2)),),(None,np.eye(3)),(None,None,np.eye(2)),(np.ones((2,3),dtype=complex),)):
                    with self.assertRaises((ValueError,TypeError)):group.rvs_array(state,*args,m=2,n=3,n_samples=0)
                readonly=np.empty((2,3,2),order='F');readonly.flags.writeable=False
                for bad in (np.empty((2,3,2)),np.empty((2,3,2),dtype=np.intc,order='F'),np.empty((3,2,2),order='F'),readonly):
                    with self.assertRaises((ValueError,TypeError)):group.rvs_array(state,m=2,n=3,n_samples=2,out={'samples':bad})
                mean=np.zeros((2,3),order='F')
                with self.assertRaises(ValueError):group.rvs_array(state,mean,m=2,n=3,out={'samples':mean.reshape((2,3,1),order='F')})
                with self.assertRaises(ValueError):group.rvs_array(state,mean)
                with self.assertRaises(TypeError):group.rvs_array([1,2,3,4,5,6])
                np.testing.assert_array_equal(state.state,saved)
                out=np.empty((2,3,2),order='F')
                result=group.rvs_array(state,m=2,n=3,n_samples=2,out={'samples':out})
                self.assertEqual(result.status,0);self.assertIs(result.outputs['samples'],out)

    def test_raw_c_flags_state_and_python_agreement(self):
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            for student in (False,True):
                group=api.matrix_t if student else api.matrix_normal
                rng=api.rng(12);state=rng.state.copy();out=np.full((2,3,2),77.,order='F')
                mean=np.arange(6.).reshape((2,3),order='F');row=np.eye(2,order='F');col=np.eye(3,order='F')
                status=C.c_int(-1);raw=group._functions['rvs_array']
                def call(m=2,n=3,count=2,flags=(1,1,1)):
                    args=[m,n,count,ptr(state),ptr(out),flags[0],ptr(mean),flags[1],ptr(row),flags[2],ptr(col)]
                    if student:args.append(5.)
                    raw(*args,C.byref(status))
                call();kwargs={'df':5} if student else {}
                expected=group.rvs_array(rng,mean,row,col,m=2,n=3,n_samples=2,**kwargs)
                self.assertEqual(status.value,0);np.testing.assert_array_equal(out,expected.outputs['samples'])
                np.testing.assert_array_equal(state,rng.state);saved=state.copy();out[:]=77
                for flags in ((2,1,1),(1,-1,1),(1,1,2)):
                    call(flags=flags);self.assertEqual(status.value,1)
                    np.testing.assert_array_equal(out,77.);np.testing.assert_array_equal(state,saved)
                for kwargs in ({'m':-1},{'n':-1},{'count':-1}):
                    call(**kwargs);self.assertEqual(status.value,1)
                    np.testing.assert_array_equal(out,77.);np.testing.assert_array_equal(state,saved)
                mean[:]=np.nan;row[:]=np.nan;col[:]=np.nan
                call(flags=(0,0,0));self.assertEqual(status.value,0);self.assertTrue(np.all(np.isfinite(out)))
                state[:]=0;out[:]=77
                for count in (0,2):
                    call(count=count);self.assertEqual(status.value,4)
                    np.testing.assert_array_equal(out,77.);np.testing.assert_array_equal(state,0)

    def test_native_c_and_cpp_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_matrix_sampling_c_api.f90',path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    client=folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 int64_t state[6]; double mean[6]={0},row[4]={1,0,0,1},col[9]={1,0,0,0,1,0,0,0,1},out[12]; int status;
 scifort_rng_seed_i64(123,state,&status); if(status) return 1;
 scifort_matrix_normal_rvs_array_f64(2,3,2,state,out,1,mean,1,row,1,col,&status);
 if(status || !isfinite(out[0]) || !isfinite(out[11])) return 2;
 scifort_matrix_t_rvs_array_f64(2,3,2,state,out,0,mean,0,row,0,col,5,&status);
 return status || !isfinite(out[0]) || !isfinite(out[11]);
}
''',encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',client])
