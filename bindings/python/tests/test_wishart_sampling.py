# SPDX-License-Identifier: MIT
"""Matrix batch layout, native streams and defined failure remainders."""
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


class WishartSamplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_streams_layout_split_and_scalar_batches(self):
        scale=np.array([[2.,.4],[.4,1.]])
        for _,_,api in self.builds:
            for group in (api.wishart,api.invwishart):
                state=api.rng(2**40);reference=api.rng(2**40)
                out=np.empty((2,2,13),order='F')
                r=group.rvs_array(state,7,scale,n_samples=13,out={'samples':out})
                self.assertEqual(r.status,0);self.assertIs(r.outputs['samples'],out)
                first=group.rvs_array(reference,7,scale,n_samples=4).outputs['samples']
                rest=group.rvs_array(reference,7,scale,n_samples=9).outputs['samples']
                np.testing.assert_array_equal(out,np.concatenate([first,rest],axis=2))
                np.testing.assert_array_equal(state.state,reference.state)
                consumed=api.rng(2**40);consumed.uniform(3*13)
                np.testing.assert_array_equal(state.state,consumed.state)
                reference=api.rng(2**40)
                individual=np.concatenate([group.rvs_array(reference,7,scale).outputs['samples'] for _ in range(13)],axis=2)
                np.testing.assert_array_equal(out,individual)
                asymmetric=scale.copy();asymmetric[0,1]=99
                a=group.rvs_array(api.rng(3),7,asymmetric,n_samples=4)
                b=group.rvs_array(api.rng(3),7,scale,n_samples=4)
                np.testing.assert_array_equal(a.outputs['samples'],b.outputs['samples'])

    def test_positive_definiteness_symmetry_and_moments(self):
        scale=np.array([[2.,.4],[.4,1.]]);size=10000;df=12
        for _,_,api in self.builds:
            for name in ('wishart','invwishart'):
                group=getattr(api,name)
                r=group.rvs_array(api.rng(543),df,scale,n_samples=size)
                self.assertEqual(r.status,0);samples=r.outputs['samples']
                np.testing.assert_allclose(samples,samples.swapaxes(0,1),rtol=2e-15,atol=2e-15)
                self.assertTrue(np.all(np.linalg.eigvalsh(np.moveaxis(samples,2,0))>0))
                mean=group.mean(df,scale);variance=group.var(df,scale)
                np.testing.assert_array_less(abs(samples.mean(axis=2)-mean),6*np.sqrt(variance/size))
                np.testing.assert_allclose(samples.var(axis=2),variance,rtol=.18,atol=.002)

    def test_empty_invalid_and_generation_failure_remainders(self):
        for _,_,api in self.builds:
            for group in (api.wishart,api.invwishart):
                state=api.rng(11);saved=state.state
                for scale in (np.full((2,2),np.nan),np.diag([1.,-1.]),np.empty((0,0))):
                    empty=group.rvs_array(state,np.nan,scale,n_samples=0)
                    self.assertEqual(empty.status,0);self.assertEqual(empty.outputs['samples'].shape,(*scale.shape,0))
                for df,scale,code in ((1.,np.eye(2),2),(np.nan,np.eye(2),2),
                                      (np.inf,np.eye(2),2),(7.,np.full((2,2),np.nan),2),
                                      (7.,np.diag([1.,-1.]),3),(7.,np.diag([1.,0.]),3)):
                    out=np.full((2,2,4),77.,order='F')
                    result=group.rvs_array(state,df,scale,n_samples=4,out={'samples':out})
                    self.assertEqual(result.status,code);self.assertTrue(np.all(np.isnan(out)))
                    np.testing.assert_array_equal(state.state,saved)
                self.assertEqual(group.rvs_array(state,7,np.empty((0,0))).status,1)
                np.testing.assert_array_equal(state.state,saved)
            state=api.rng(11);expected=api.rng(11);expected.uniform(1)
            out=np.full((1,1,3),77.,order='F')
            tiny=api.invwishart.rvs_array(state,1e-300,[[1.]],n_samples=3,out={'samples':out})
            self.assertEqual(tiny.status,3);self.assertTrue(np.all(np.isnan(out)))
            np.testing.assert_array_equal(state.state,expected.state)
            # Preserve the existing native Wishart underflow behavior explicitly.
            tiny=api.wishart.rvs_array(api.rng(11),1e-300,[[1.]])
            self.assertEqual(tiny.status,0);self.assertTrue(np.all(np.isnan(tiny.outputs['samples'])))

            # A deterministic later failure retains the completed first matrix.
            state=api.rng(2**40+123);state.uniform(10)
            reference=api.rng(2**40+123);reference.uniform(10)
            first=api.invwishart.rvs_array(reference,.0001,[[1.]])
            self.assertEqual(first.status,0)
            reference.uniform(1)
            out=np.full((1,1,3),77.,order='F')
            result=api.invwishart.rvs_array(state,.0001,[[1.]],n_samples=3,out={'samples':out})
            self.assertEqual(result.status,3)
            np.testing.assert_array_equal(out[:,:,0],first.outputs['samples'][:,:,0])
            self.assertTrue(np.all(np.isnan(out[:,:,1:])))
            np.testing.assert_array_equal(state.state,reference.state)

    def test_shapes_counts_and_buffers(self):
        for _,_,api in self.builds:
            for group in (api.wishart,api.invwishart):
                state=api.rng(12);saved=state.state
                for bad in (-1,True,1.5,2**31):
                    with self.assertRaises((TypeError,ValueError)):group.rvs_array(state,7,np.eye(2),n_samples=bad)
                for scale in (np.ones((2,3)),[1.,2.],np.ones((2,2,1)),np.eye(2,dtype=complex)):
                    with self.assertRaises((TypeError,ValueError)):group.rvs_array(state,7,scale,n_samples=0)
                readonly=np.empty((2,2,3),order='F');readonly.flags.writeable=False
                for out in (np.empty((2,2,3)),np.empty((2,2,3),dtype=np.intc,order='F'),
                            np.empty((2,3,2),order='F'),readonly):
                    with self.assertRaises((TypeError,ValueError)):group.rvs_array(state,7,np.eye(2),n_samples=3,out={'samples':out})
                scale=np.eye(2,order='F')
                with self.assertRaises(ValueError):group.rvs_array(state,7,scale,out={'samples':scale.reshape(2,2,1,order='F')})
                np.testing.assert_array_equal(state.state,saved)

    def test_raw_c_guards_and_python_agreement(self):
        ptr=lambda x:x.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            for group in (api.wishart,api.invwishart):
                rng=api.rng(8);state=rng.state.copy();scale=np.array([[2.,.4],[.4,1.]],order='F')
                out=np.full((2,2,3),77.,order='F');status=C.c_int(-1);raw=group._functions['rvs_array']
                def call(p=2,count=3):raw(p,count,ptr(state),7.,ptr(scale),ptr(out),C.byref(status))
                call();expected=group.rvs_array(rng,7,scale,n_samples=3)
                self.assertEqual(status.value,0);np.testing.assert_array_equal(out,expected.outputs['samples'])
                np.testing.assert_array_equal(state,rng.state);saved=state.copy();out[:]=77
                for args in ((-1,3),(2,-1)):
                    call(*args);self.assertEqual(status.value,1)
                    np.testing.assert_array_equal(state,saved);np.testing.assert_array_equal(out,77.)
                state[:]=0
                for count in (0,3):
                    call(count=count);self.assertEqual(status.value,4)
                    np.testing.assert_array_equal(state,0);np.testing.assert_array_equal(out,77.)

    def test_native_c_and_cpp_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_wishart_sampling_c_api.f90',path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    client=folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 int64_t state[6]; double scale[4]={2,.4,.4,1},out[8]; int status;
 scifort_rng_seed_i64(123,state,&status); if(status) return 1;
 scifort_wishart_rvs_array_f64(2,2,state,7,scale,out,&status);
 if(status || out[0]<=0 || out[3]<=0 || out[4]<=0 || out[7]<=0) return 2;
 scifort_invwishart_rvs_array_f64(2,2,state,7,scale,out,&status);
 return status || !isfinite(out[0]) || out[0]<=0;
}
''',encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',client])
