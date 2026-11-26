# SPDX-License-Identifier: MIT
"""Explicit-state streams, generated Dirichlet sampling and transport error contracts."""
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
from scifort import Scifort, RandomState
from test_multivariate_t import run
ROOT=Path(__file__).resolve().parents[3]


class RandomTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran'):raise unittest.SkipTest('gfortran required')
        cls.builds=[]
        for checked in (True,False):
            path,metadata=build(checked=checked)
            cls.builds.append((path,metadata,Scifort(path,dll_directories=metadata['dll_directories'])))

    def test_uniform_reproducibility_snapshots_and_seed_range(self):
        for _,_,api in self.builds:
            for seed in (0,1,-1,2**40,-2**63,2**63-1):
                first=api.rng(seed);second=RandomState(api,seed)
                self.assertEqual(first.state.dtype,np.dtype(np.int64))
                self.assertEqual(first.state.shape,(6,))
                expected=second.uniform(20)
                actual=np.concatenate([first.uniform(7),first.uniform(13)])
                np.testing.assert_array_equal(actual,expected)
                self.assertTrue(np.all((actual>0)&(actual<1)))
                np.testing.assert_array_equal(first.state,second.state)
                snapshot=first.state
                expected=first.uniform(12)
                first.set_state(snapshot)
                np.testing.assert_array_equal(first.uniform(12),expected)
                copy=first.state;copy[:]=0
                self.assertFalse(np.all(first.state==0))
                snapshot=first.state
                self.assertEqual(first.uniform(0).shape,(0,))
                np.testing.assert_array_equal(first.state,snapshot)

    def test_seed_state_uniform_shapes_and_buffer_errors(self):
        for _,_,api in self.builds:
            for seed in (True,1.,-2**63-1,2**63,np.uint64(2**64-1)):
                with self.assertRaises((ValueError,TypeError)):api.rng(seed)
            state=api.rng(7);saved=state.state
            invalid=[np.zeros(6,dtype=np.int64),[-1,1,1,1,1,1],
                     [4294967087,1,1,1,1,1],[1,1,1,4294944443,1,1],
                     [0,0,0,1,1,1],[1,1,1,0,0,0],[1,2],np.ones(6),
                     np.full(6,2**64-1,dtype=np.uint64)]
            for values in invalid:
                with self.assertRaises((ValueError,TypeError)):state.set_state(values)
                np.testing.assert_array_equal(state.state,saved)
            restored=saved.astype('>i8');state.set_state(restored)
            np.testing.assert_array_equal(state.state,saved)
            for n in (-1,2**31,True,1.2):
                with self.assertRaises((ValueError,TypeError)):state.uniform(n)
            readonly=np.empty(3);readonly.flags.writeable=False
            for out in (np.empty(4),np.empty(3,dtype=np.float32),np.empty((3,1)),readonly,np.empty(6)[::2]):
                with self.assertRaises((ValueError,TypeError)):state.uniform(3,out=out)
            np.testing.assert_array_equal(state.state,saved)
            out=np.empty(3)
            self.assertIs(state.uniform(3,out=out),out)

    def test_dirichlet_reproducibility_state_consumption_and_moments(self):
        alpha=np.array([2.,3.,4.]);n=2048
        for _,_,api in self.builds:
            state=api.rng(123);reference=api.rng(123)
            out=np.empty((3,n),order='F')
            result=api.dirichlet.rvs_array(state,alpha,n_samples=n,out={'samples':out})
            self.assertEqual(result.status,0);self.assertIs(result.outputs['samples'],out)
            reference.uniform(3*n)
            np.testing.assert_array_equal(state.state,reference.state)
            self.assertTrue(np.all(np.isfinite(out)) and np.all(out>=0))
            np.testing.assert_allclose(out.sum(axis=0),1,rtol=0,atol=4e-16)
            mean=alpha/alpha.sum();var=mean*(1-mean)/(alpha.sum()+1)
            np.testing.assert_array_less(abs(out.mean(axis=1)-mean),6*np.sqrt(var/n))
            np.testing.assert_allclose(out.var(axis=1),var,rtol=.2)
            snapshot=state.state
            first=api.dirichlet.rvs_array(state,alpha).outputs['samples']
            self.assertEqual(first.shape,(3,1))
            state.set_state(snapshot)
            np.testing.assert_array_equal(api.dirichlet.rvs_array(state,alpha).outputs['samples'],first)
            self.assertTrue(np.all(api.dirichlet.rvs_array(state,[2.],n_samples=4).outputs['samples']==1))

    def test_empty_invalid_parameters_and_numerical_failure_state(self):
        for _,_,api in self.builds:
            state=api.rng(17);saved=state.state
            for alpha in ([],[-1.,np.nan]):
                result=api.dirichlet.rvs_array(state,alpha,n_samples=0)
                self.assertEqual(result.status,0)
                self.assertEqual(result.outputs['samples'].shape,(len(alpha),0))
                np.testing.assert_array_equal(state.state,saved)
            for alpha in ([],[-1.,2.],[np.nan,2.],[np.inf,2.],[0.,2.]):
                result=api.dirichlet.rvs_array(state,alpha,n_samples=2)
                self.assertEqual(result.status,1)
                self.assertTrue(np.all(np.isnan(result.outputs['samples'])))
                np.testing.assert_array_equal(state.state,saved)
            # Positive tiny shapes underflow the inverse-Gamma draws to zero.
            reference=api.rng(17);reference.uniform(2)
            result=api.dirichlet.rvs_array(state,[1e-300,1e-300],n_samples=2)
            self.assertEqual(result.status,1)
            self.assertTrue(np.all(np.isnan(result.outputs['samples'])))
            np.testing.assert_array_equal(state.state,reference.state)

    def test_sampler_validation_and_state_ownership(self):
        for path,metadata,api in self.builds:
            state=api.rng(5);saved=state.state
            other=Scifort(path,dll_directories=metadata['dll_directories'])
            for bad in (state.state,other.rng(5)):
                with self.assertRaises(TypeError):api.dirichlet.rvs_array(bad,[2.,3.])
            for count in (-1,True,2**31,1.5):
                with self.assertRaises((ValueError,TypeError)):api.dirichlet.rvs_array(state,[2.,3.],n_samples=count)
            for alpha in ([[2.,3.]],np.ones((1,2,1)),[2.+1j,3.]):
                with self.assertRaises((ValueError,TypeError)):api.dirichlet.rvs_array(state,alpha)
            readonly=np.empty((2,2),order='F');readonly.flags.writeable=False
            for out in (np.empty((3,2),order='F'),np.empty((2,2)),np.empty((2,2),dtype=np.float32,order='F'),readonly):
                with self.assertRaises((ValueError,TypeError)):
                    api.dirichlet.rvs_array(state,[2.,3.],n_samples=2,out={'samples':out})
            alpha=np.array([2.,3.])
            with self.assertRaises(ValueError):api.dirichlet.rvs_array(state,alpha,out={'samples':alpha.reshape((2,1))})
            np.testing.assert_array_equal(state.state,saved)
            for variant in (alpha.astype('>f8'),np.array([2.,9.,3.])[::2]):
                result=api.dirichlet.rvs_array(state,variant)
                self.assertEqual(result.status,0)

    def test_raw_c_invalid_transport_preserves_outputs(self):
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            rng=api.rng(123);saved=rng.state;state=saved.copy()
            out=np.full((2,2),77.,order='F');alpha=np.array([2.,3.]);status=C.c_int(-1)
            raw=api.dirichlet._functions['rvs_array']
            for values in (np.zeros(6,dtype=np.int64),np.array([4294967087,1,1,1,1,1],dtype=np.int64)):
                before=values.copy()
                for count in (0,2):
                    raw(2,count,ptr(values),ptr(out),ptr(alpha),C.byref(status))
                    self.assertEqual(status.value,3)
                    np.testing.assert_array_equal(values,before);np.testing.assert_array_equal(out,np.full((2,2),77.))
            raw(-1,2,ptr(state),ptr(out),ptr(alpha),C.byref(status))
            self.assertEqual(status.value,2);np.testing.assert_array_equal(state,saved)
            np.testing.assert_array_equal(out,np.full((2,2),77.))
            raw_uniform=rng._uniform
            raw_uniform(-1,ptr(state),ptr(out),C.byref(status))
            self.assertEqual(status.value,1);np.testing.assert_array_equal(state,saved)
            state[:]=0
            raw_uniform(0,ptr(state),ptr(out),C.byref(status))
            self.assertEqual(status.value,3);np.testing.assert_array_equal(out,np.full((2,2),77.))

    def test_discrete_sampling_streams_moments_and_counts(self):
        count=4096
        for _,_,api in self.builds:
            for group,params,n,mean,var in (
                (api.multinomial,[.2,.3,.5],12,np.array([2.4,3.6,6.]),12*np.array([.16,.21,.25])),
                (api.multivariate_hypergeom,[10,20,30],12,np.array([2.,4.,6.]),
                 12*np.array([1/6,1/3,1/2])*np.array([5/6,2/3,1/2])*48/59)):
                state=api.rng(471);reference=api.rng(471)
                def draw(rng,size,out=None):
                    if group is api.multinomial:
                        return group.rvs_array(rng,n,params,n_samples=size,out=out)
                    return group.rvs_array(rng,params,n,n_samples=size,out=out)
                out=np.empty((3,count),dtype=np.intc,order='F')
                result=draw(state,count,{'samples':out})
                self.assertEqual(result.status,0);self.assertIs(result.outputs['samples'],out)
                first=draw(reference,17).outputs['samples']
                rest=draw(reference,count-17).outputs['samples']
                np.testing.assert_array_equal(out,np.concatenate([first,rest],axis=1))
                np.testing.assert_array_equal(state.state,reference.state)
                self.assertTrue(np.all(out>=0));np.testing.assert_array_equal(out.sum(axis=0),n)
                if group is api.multivariate_hypergeom:self.assertTrue(np.all(out<=np.array(params)[:,None]))
                np.testing.assert_array_less(abs(out.mean(axis=1)-mean),6*np.sqrt(var/count))
                np.testing.assert_allclose(out.var(axis=1),var,rtol=.15)
                self.assertEqual(draw(state,1).outputs['samples'].shape,(3,1))
                saved=state.state;self.assertEqual(draw(state,0).status,0)
                np.testing.assert_array_equal(state.state,saved)
                readonly=np.empty((3,2),dtype=np.intc,order='F');readonly.flags.writeable=False
                for bad in (np.empty((3,2)),np.empty((3,2),dtype=np.int64,order='F'),
                            np.empty((3,2),dtype=np.intc),readonly,np.empty((2,2),dtype=np.intc,order='F')):
                    with self.assertRaises((TypeError,ValueError)):draw(state,2,{'samples':bad})
                np.testing.assert_array_equal(state.state,saved)

    def test_discrete_parameters_empty_batches_and_overflow(self):
        for _,_,api in self.builds:
            state=api.rng(15);saved=state.state
            for p in ([],[np.nan,-1.]):
                self.assertEqual(api.multinomial.rvs_array(state,-1,p,n_samples=0).status,0)
            for p in ([],[np.nan,1.],[-1.,2.],[np.inf,0.]):
                r=api.multinomial.rvs_array(state,5,p,n_samples=2)
                self.assertEqual(r.status,1);self.assertTrue(np.all(r.outputs['samples']==0))
            self.assertEqual(api.multinomial.rvs_array(state,-1,[.5,.5]).status,1)
            for m,n in (([-1,2],1),([1,2],4),([1,2],-1)):
                for count in (0,2):
                    r=api.multivariate_hypergeom.rvs_array(state,m,n,n_samples=count)
                    self.assertEqual(r.status,1);self.assertTrue(np.all(r.outputs['samples']==0))
            for m in ([],[2**30,2**30],np.array([2**64-1],dtype=np.uint64),[1.,2.]):
                with self.assertRaises((ValueError,TypeError)):
                    api.multivariate_hypergeom.rvs_array(state,m,1)
            for group in (api.multinomial,api.multivariate_hypergeom):
                args=(state,2**31-1,[.5,.5]) if group is api.multinomial else (state,[1,2],2**31-1)
                with self.assertRaises(ValueError):group.rvs_array(*args)
            np.testing.assert_array_equal(state.state,saved)
            np.testing.assert_array_equal(api.multinomial.rvs_array(state,5,[9.],n_samples=3).outputs['samples'],5)
            np.testing.assert_array_equal(api.multivariate_hypergeom.rvs_array(state,[5],5,n_samples=3).outputs['samples'],5)
            np.testing.assert_array_equal(state.state,saved)
            # Degenerate zero draws do not consume uniforms.
            self.assertTrue(np.all(api.multinomial.rvs_array(state,0,[.5,.5]).outputs['samples']==0))
            self.assertTrue(np.all(api.multivariate_hypergeom.rvs_array(state,[0,0],0).outputs['samples']==0))
            np.testing.assert_array_equal(state.state,saved)

    def test_discrete_raw_c_guards_and_python_agreement(self):
        ptr=lambda a:a.ctypes.data_as(C.c_void_p)
        for _,_,api in self.builds:
            for group in (api.multinomial,api.multivariate_hypergeom):
                rng=api.rng(123);raw_state=rng.state.copy();saved=raw_state.copy()
                output=np.full((2,3),77,dtype=np.intc,order='F');status=C.c_int(-1)
                param=np.array([.4,.6]) if group is api.multinomial else np.array([4,6],dtype=np.intc)
                raw=group._functions['rvs_array']
                def call(k=2,columns=3,n=4):
                    if group is api.multinomial:raw(k,columns,ptr(raw_state),ptr(output),n,ptr(param),C.byref(status))
                    else:raw(k,columns,ptr(raw_state),ptr(param),n,ptr(output),C.byref(status))
                call()
                expected=(group.rvs_array(rng,4,param,n_samples=3) if group is api.multinomial
                          else group.rvs_array(rng,param,4,n_samples=3))
                self.assertEqual(status.value,0);np.testing.assert_array_equal(output,expected.outputs['samples'])
                np.testing.assert_array_equal(raw_state,rng.state)
                saved=raw_state.copy();output[:]=77
                for k,columns,n in ((-1,3,4),(2,-1,4),(2,3,2**31-1)):
                    call(k,columns,n);self.assertEqual(status.value,2)
                    np.testing.assert_array_equal(output,77);np.testing.assert_array_equal(raw_state,saved)
                if group is api.multivariate_hypergeom:
                    param[:]=2**30;call();self.assertEqual(status.value,2)
                    np.testing.assert_array_equal(output,77);np.testing.assert_array_equal(raw_state,saved)
                raw_state[:]=0
                for columns in (0,3):
                    if group is api.multivariate_hypergeom:param[:]=[4,6]
                    call(columns=columns);self.assertEqual(status.value,3)
                    np.testing.assert_array_equal(output,77);np.testing.assert_array_equal(raw_state,0)

    def test_native_and_c_clients(self):
        for path,metadata,_ in self.builds:
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp);env=os.environ.copy();env['PATH']=str(path.parent)+os.pathsep+env.get('PATH','')
                flags=['-O0','-g','-fcheck=all','-fbacktrace'] if '-fcheck=all' in metadata['flags'] else ['-O3']
                exe=folder/'native.exe'
                run([metadata['compiler_path'],*flags,'-I',path.parent,ROOT/'test/test_rng_c_api.f90',path,'-o',exe]);run([exe],env=env)
                if shutil.which('gcc'):
                    client=folder/'client.c'
                    client.write_text('''#include "scifort.h"
#include <math.h>
int main(void) {
 int counts[4],m[2]={4,6};
 int64_t state[6],saved[6]; double u[6],alpha[2]={2,3},samples[4]; int status;
 scifort_rng_seed_i64(INT64_MIN,state,&status); if(status) return 1;
 for(int i=0;i<6;i++) saved[i]=state[i];
 scifort_rng_uniform_f64(6,state,u,&status); if(status) return 2;
 for(int i=0;i<6;i++) if(!(u[i]>0 && u[i]<1)) return 3;
 scifort_dirichlet_rvs_array_f64(2,2,state,samples,alpha,&status);
 if(status || fabs(samples[0]+samples[1]-1)>1e-14 || fabs(samples[2]+samples[3]-1)>1e-14) return 4;
 scifort_rng_seed_i64(INT64_C(1099511627776),state,&status); if(status) return 5;
 scifort_multinomial_rvs_array_i32(2,2,state,counts,4,alpha,&status);
 /* alpha={2,3} is invalid after final-entry adjustment. */
 if(status!=1) return 6;
 scifort_multivariate_hypergeom_rvs_array_i32(2,2,state,m,4,counts,&status);
 if(status || counts[0]+counts[1]!=4 || counts[2]+counts[3]!=4) return 7;
 for(int i=0;i<6;i++) state[i]=0;
 scifort_rng_uniform_f64(0,state,u,&status);
 return status!=3 || saved[0]<=0;
}
''',encoding='ascii')
                    exe=folder/'c.exe'
                    run(['gcc','-std=c11','-Wall','-Wextra','-Werror','-I',ROOT/'include',client,path,'-o',exe]);run([exe],env=env)
                    if shutil.which('g++'):run(['g++','-std=c++17','-Wall','-Wextra','-Werror','-I',ROOT/'include','-fsyntax-only',client])
